%% beadAnalysisBatch.m
% Batch script to run estimateBeadDistance across multiple MovieLists
% (conditions) and compare bead statistics between conditions.
%
% Each movie image is divided into non-overlapping 500x500 px tiles.
% Dark/edge tiles (mean intensity < 30% of image mean) are skipped
% automatically, avoiding uneven-illumination border regions.
% Statistics are computed per tile and pooled per movie and per condition.
%
% Modeled after adhesionAnalysisBatch.m
% Sangyoon Han. May 2020 / Updated for tiled ROI analysis, 2026.
%
% Outputs (saved to AnalysisSummary_BeadAnalysis<specificName>/):
%   Figs/        - all comparison figures (boxplots, histograms)
%   Data/        - CSV tables and .mat summary
%   MovieImages/ - per-tile overlay images

%% Open necessary MLs
clear

% ---- PSF sigma override ---------------------------------------------------
% Set to [] to use the value derived from each channel's emissionWavelength_,
% NA, and pixelSize_. Set to a positive number (pixels) to override for all
% movies ? useful when metadata-derived sigma is wrong.
%   Theoretical: sigma = 0.21 * lambda_em(nm) / NA / pixSize(nm)
%   e.g. 560nm / 1.45 / 108nm = 0.75 px  (too small for actual PSF)
%   Override to ~1.5-1.6 px to match the observed bead size.
psfSigmaOverride = 1.6;   % px ? set to [] to use metadata-derived value

% ---- Tiling settings ------------------------------------------------------
tileSize    = 500;   % px ? side length of each square ROI tile
minMeanFrac = 0.3;   % skip tiles whose mean < this fraction of image mean

[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat, ...
    specificName, ~, MLdirect] = chooseSelectedFolders;
nameList = groupNames';

%% Output folders
rootAnalysis = fileparts(pathAnalysisAll{1});
summaryPath  = [rootAnalysis '/AnalysisSummary_BeadAnalysis' specificName];
ii = 0;
while exist(summaryPath, 'dir')
    ii = ii + 1;
    summaryPath = [rootAnalysis '/AnalysisSummary_BeadAnalysis' specificName num2str(ii)];
end
figPath  = [summaryPath '/Figs'];  mkdir(figPath)
dataPath = [summaryPath '/Data'];  mkdir(dataPath)

save([rootAnalysis filesep 'selectedFolders' specificName '.mat'], ...
    'rootAnalysis', 'pathAnalysisAll', 'MLNames', 'groupNames')

%% Load MovieLists for each condition
numConditions = numel(pathAnalysisAll);
MLAll = MovieList.empty(0, numConditions);
for k = 1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end

%% Per-condition storage
N                    = zeros(numConditions, 1);
beadDistAllGroup     = cell(numConditions, 1);
beadDensityGroup     = cell(numConditions, 1);
aggregationLevelGroup= cell(numConditions, 1);
meanSegAreaUm2Group  = cell(numConditions, 1);

%% Main loop: condition -> movie -> tile
for ii = 1:numConditions
    curML = MLAll(ii);
    N(ii) = numel(curML.movieDataFile_);

    beadDistAllCond      = {};
    beadDensityCond      = [];
    aggregationLevelCond = [];
    meanSegAreaUm2Cond   = [];

    for k = 1:N(ii)
        curMD = MovieData.load(curML.movieDataFile_{k});

        % Identify bead channel
        iTFM = curMD.getPackageIndex('TFMPackage');
        if isempty(iTFM)
            iBeadChan = 1;
        else
            tPack       = curMD.getPackage(iTFM);
            dispCalProc = tPack.getProcess(2);
            iBeadChan   = dispCalProc.checkChannelOutput;
        end
        beadChan = curMD.getChannel(iBeadChan);

        % Resolve psfSigma
        if ~isempty(psfSigmaOverride) && psfSigmaOverride > 0
            psfSigma = psfSigmaOverride;
        else
            try,   psfSigma = beadChan.psfSigma_; catch, psfSigma = []; end
            if isempty(psfSigma) || psfSigma <= 0
                psfSigma = 0.21 * beadChan.emissionWavelength_ / ...
                           curMD.numAperture_ / curMD.pixelSize_;
            end
        end

        % Load full image
        fullImg = beadChan.loadImage(1);

        % Tile the image; skip dark border regions automatically
        [tiles, tileOrigins, validMask] = tileBeadImage(fullImg, tileSize, minMeanFrac);
        nTiles      = numel(tiles);
        nValidTiles = sum(validMask);
        fprintf('  [%s] movie %d/%d ? %d/%d valid tiles\n', ...
                groupNames{ii}, k, N(ii), nValidTiles, nTiles)

        % Per-movie image output folder
        [~, mdName] = fileparts(fileparts(curMD.getFullPath));
        movieImgDir = fullfile(summaryPath, 'MovieImages', groupNames{ii}, mdName);
        if ~exist(movieImgDir, 'dir'), mkdir(movieImgDir); end

        % Sub-loop over valid tiles
        for t = 1:nTiles
            if ~validMask(t), continue, end

            tileImg   = tiles{t};
            tileRow   = tileOrigins(t, 1);
            tileCol   = tileOrigins(t, 2);
            tileLabel = sprintf('r%04d_c%04d', tileRow, tileCol);
            tileDir   = fullfile(movieImgDir, tileLabel);
            if ~exist(tileDir, 'dir'), mkdir(tileDir); end

            try
                [beadDist, beadDens, aggLevel, segArea] = ...
                    estimateBeadDistance(tileImg, curMD.pixelSize_, ...
                                         psfSigma, tileDir);
            catch ME
                warning('beadAnalysisBatch:tileFailed', ...
                    '[%s] movie %d tile %s: %s', ...
                    groupNames{ii}, k, tileLabel, ME.message)
                beadDist = NaN; beadDens = NaN;
                aggLevel = NaN; segArea  = NaN;
            end

            beadDistAllCond{end+1,1}      = beadDist(:);  %#ok<AGROW>
            beadDensityCond(end+1,1)      = beadDens;     %#ok<AGROW>
            aggregationLevelCond(end+1,1) = aggLevel;     %#ok<AGROW>
            meanSegAreaUm2Cond(end+1,1)   = segArea;      %#ok<AGROW>
        end
    end

    beadDistAllGroup{ii}      = cell2mat(beadDistAllCond);
    beadDensityGroup{ii}      = beadDensityCond;
    aggregationLevelGroup{ii} = aggregationLevelCond;
    meanSegAreaUm2Group{ii}   = meanSegAreaUm2Cond;

    fprintf('  => %s: %d valid tiles, %d distance measurements\n', ...
            groupNames{ii}, numel(beadDensityCond), numel(beadDistAllGroup{ii}))
end

%% =========================================================
%% PLOTTING
%% =========================================================

colors = lines(numConditions);

%% 1) Bead-to-bead distance ? histogram per condition (overlaid)
h1 = figure(1); clf; hold on
legendEntries = cell(numConditions, 1);
for ii = 1:numConditions
    d = beadDistAllGroup{ii};
    d = d(isfinite(d));
    histogram(d, 'FaceColor', colors(ii,:), 'FaceAlpha', 0.4, ...
              'EdgeColor', 'none', 'Normalization', 'probability')
    legendEntries{ii} = sprintf('%s (mean=%.2f µm, n=%d tiles)', ...
        nameList{ii}, mean(d), numel(beadDensityGroup{ii}));
end
legend(legendEntries, 'Location', 'northeast')
xlabel('Nearest-neighbour distance (µm)')
ylabel('Probability')
title('Bead-to-bead distance by condition')
box on
saveFigure(h1, figPath, 'beadDistHist')

%% 2) Bead-to-bead distance ? boxplot
h2 = figure(2); clf
boxPlotCellArray(beadDistAllGroup, nameList, 1, false, true)
ylabel('Nearest-neighbour distance (µm)')
title('Bead-to-bead distance')
saveFigure(h2, figPath, 'beadDistBox')

%% 3) Bead density ? boxplot (per tile)
h3 = figure(3); clf
boxPlotCellArray(beadDensityGroup, nameList, 1, false, true)
ylabel('Bead density (#/µm²)')
title('Bead density per tile')
saveFigure(h3, figPath, 'beadDensity')

%% 4) Aggregation level ? boxplot (per tile)
h4 = figure(4); clf
boxPlotCellArray(aggregationLevelGroup, nameList, 1, false, true)
ylabel('Aggregation level (beads/region)')
title('Bead aggregation level per tile')
hold on
yline(1, '--k', 'Non-aggregated', 'LabelHorizontalAlignment', 'left')
saveFigure(h4, figPath, 'beadAggregation')

%% 5) Mean segmentation area ? boxplot (per tile)
h5 = figure(5); clf
boxPlotCellArray(meanSegAreaUm2Group, nameList, 1, false, true)
ylabel('Mean segmentation area (µm²)')
title('Mean mask segmentation area per tile')
saveFigure(h5, figPath, 'beadSegArea')

%% 6) Density vs aggregation level scatter (one dot per tile)
h6 = figure(6); clf; hold on
for ii = 1:numConditions
    scatter(beadDensityGroup{ii}, aggregationLevelGroup{ii}, 30, ...
            colors(ii,:), 'filled', 'DisplayName', nameList{ii}, ...
            'MarkerFaceAlpha', 0.5)
end
xlabel('Bead density (#/µm²)')
ylabel('Aggregation level (beads/region)')
title('Density vs Aggregation level (per tile)')
legend('Location', 'best')
box on
saveFigure(h6, figPath, 'densityVsAggregation')

%% =========================================================
%% SAVE DATA
%% =========================================================

% Per-condition pooled distance CSVs
for ii = 1:numConditions
    d = beadDistAllGroup{ii};
    writetable(table(d, 'VariableNames', {'distanceUm'}), ...
        fullfile(dataPath, ['beadDist_' nameList{ii} '.csv']))
end

% Per-tile summary table ? collect tile labels in a second pass
% (lightweight: only loads image to count tiles, no detection)
allCondNames  = {};
allMovieNames = {};
allTileLabels = {};
allDensity    = [];
allAggLevel   = [];
allSegArea    = [];

for ii = 1:numConditions
    curML = MLAll(ii);
    for k = 1:numel(curML.movieDataFile_)
        curMD = MovieData.load(curML.movieDataFile_{k});
        [~, mdName] = fileparts(fileparts(curMD.getFullPath));
        fullImg = curMD.getChannel(1).loadImage(1);
        [~, tileOrigins, validMask] = tileBeadImage(fullImg, tileSize, minMeanFrac);
        validIdx = find(validMask);
        for t = 1:numel(validIdx)
            ti = validIdx(t);
            allCondNames{end+1,1}  = groupNames{ii};
            allMovieNames{end+1,1} = mdName;
            allTileLabels{end+1,1} = sprintf('r%04d_c%04d', ...
                tileOrigins(ti,1), tileOrigins(ti,2));
        end
    end
    allDensity  = [allDensity;  beadDensityGroup{ii}];       %#ok<AGROW>
    allAggLevel = [allAggLevel; aggregationLevelGroup{ii}];  %#ok<AGROW>
    allSegArea  = [allSegArea;  meanSegAreaUm2Group{ii}];     %#ok<AGROW>
end

summaryTable = table(allCondNames, allMovieNames, allTileLabels, ...
                     allDensity, allAggLevel, allSegArea, ...
    'VariableNames', {'Condition','Movie','Tile', ...
                      'BeadDensity_per_um2','AggregationLevel','MeanSegArea_um2'});
writetable(summaryTable, fullfile(dataPath, 'beadSummaryAllTiles.csv'))

% Save workspace
save(fullfile(dataPath, 'beadAnalysisData.mat'), '-v7.3')

disp('=== beadAnalysisBatch complete. Results saved to:')
disp(summaryPath)
close all

%% =========================================================
%% Local helper: save figure
%% =========================================================
function saveFigure(hFig, figPath, baseName)
    hgexport(hFig, fullfile(figPath, baseName), hgexport('factorystyle'), 'Format', 'eps')
    hgsave(hFig,   fullfile(figPath, baseName), '-v7.3')
    print(hFig,    fullfile(figPath, [baseName '.tif']), '-dtiff')
end