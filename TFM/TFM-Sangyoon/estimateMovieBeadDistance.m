function [] = estimateMovieBeadDistance(ML)
%function [] = estimateMovieBeadDistance(ML) reads ML and MDs, pulls up one
% bead image, then runs estimateBeadDistance and quantifies statistics about
% bead distance, density, aggregation level, and segmentation area.
%
%   input:
%       ML:         MovieList
%
%   output:
%       Graphs and data stored where ML is located.
%
% Sangyoon Han. May 2020
% Updated: aggregationLevel, meanSegAreaUm2 boxplots; per-movie image saving.

%% Output folder
folderPath = [ML.getPath filesep 'BeadDistance_' ML.getFilename];
folderPath = folderPath(1:end-4);
if ~exist(folderPath, 'dir')
    mkdir(folderPath)
end

%% Read ML
numMovies = numel(ML.movieDataFile_);

%% Per-movie storage
beadDistAllCell      = cell(numMovies, 1);
beadDensityAll       = zeros(numMovies, 1);
aggregationLevelAll  = zeros(numMovies, 1);
meanSegAreaUm2All    = zeros(numMovies, 1);

for ii = 1:numMovies
    % Load MovieData
    curMD = MovieData.load(ML.movieDataFile_{ii});

    % Identify bead channel
    iTFM = curMD.getPackageIndex('TFMPackage');
    if isempty(iTFM)
        disp(['No TFM was done in this movie: ' curMD.getFilename '.'])
        disp('Using the first channel for bead image.');
        iBeadChan = 1;
    else
        tPack        = curMD.getPackage(iTFM);
        dispCalProc  = tPack.getProcess(2);
        iBeadChan    = dispCalProc.checkChannelOutput;
    end
    beadChan = curMD.getChannel(iBeadChan);

    % psfSigma_ is read-only/derived in some u-Track versions.
    % Try reading it directly; fall back to computing from first principles.
    try
        psfSigma = beadChan.psfSigma_;
    catch
        psfSigma = [];
    end
    if isempty(psfSigma) || psfSigma <= 0
        psfSigma = 0.21 * beadChan.emissionWavelength_ / ...
                   curMD.numAperture_ / curMD.pixelSize_;
    end

    % Per-movie output subfolder for overlay images
    [~, mdName] = fileparts(curMD.getFilename);
    movieImgFolder = fullfile(folderPath, mdName);
    if ~exist(movieImgFolder, 'dir')
        mkdir(movieImgFolder)
    end

    % Run detection + quantification + save overlay images
    [beadDistAll, beadDensity, aggregationLevel, meanSegAreaUm2] = ...
        estimateBeadDistance(beadChan.loadImage(1), ...
                             curMD.pixelSize_, ...
                             psfSigma, ...
                             movieImgFolder);

    % Store results
    beadDistAllCell{ii}     = beadDistAll;
    beadDensityAll(ii)      = beadDensity;
    aggregationLevelAll(ii) = aggregationLevel;
    meanSegAreaUm2All(ii)   = meanSegAreaUm2;
end

%% ---- Bead distance histogram ------------------------------------------
beadDistGroup = cell2mat(beadDistAllCell);

h1    = figure(1); clf
hHist = histogram(beadDistGroup);
xlabel('Distance (\mum)')
ylabel('Number of occurrence')
[maxCount, iMaxCount] = max(hHist.BinCounts);
text(hHist.BinEdges(iMaxCount+1) + 0.5*hHist.BinWidth, maxCount*0.95, ...
     ['mean: ' num2str(mean(beadDistGroup), 2) ' \mum']);
text(hHist.BinEdges(iMaxCount+1) + 0.5*hHist.BinWidth, maxCount*0.85, ...
     ['std: '  num2str(std(beadDistGroup),  2) ' \mum']);
title('Bead-to-bead distance')

hgexport(h1, strcat(folderPath, '/beadDist'), hgexport('factorystyle'), 'Format', 'eps')
hgsave(h1,   strcat(folderPath, '/beadDist'), '-v7.3')
print(h1,    strcat(folderPath, '/beadDist.tif'), '-dtiff')
writetable(table(beadDistGroup), strcat(folderPath, '/beadDist.csv'))

%% ---- Bead density boxplot ---------------------------------------------
beadDensityAllCell = {beadDensityAll};

h2 = figure(2); clf
boxPlotCellArray(beadDensityAllCell, {'All movies'}, 1, 0, 1);
ylabel('Bead density (#/\mum^2)')
title('Bead density')

hgexport(h2, strcat(folderPath, '/beadDensity'), hgexport('factorystyle'), 'Format', 'eps')
hgsave(h2,   strcat(folderPath, '/beadDensity'), '-v7.3')
print(h2,    strcat(folderPath, '/beadDensity.tif'), '-dtiff')
writetable(table(beadDensityAll), strcat(folderPath, '/beadDensity.csv'))

%% ---- Aggregation level boxplot ----------------------------------------
aggregationLevelAllCell = {aggregationLevelAll};

h3 = figure(3); clf
boxPlotCellArray(aggregationLevelAllCell, {'All movies'}, 1, 0, 1);
ylabel('Aggregation level (beads / region)')
title('Bead aggregation level')
% Reference line at 1 = perfectly non-aggregated
hold on
yline(1, '--k', 'Non-aggregated', 'LabelHorizontalAlignment', 'left')

hgexport(h3, strcat(folderPath, '/beadAggregation'), hgexport('factorystyle'), 'Format', 'eps')
hgsave(h3,   strcat(folderPath, '/beadAggregation'), '-v7.3')
print(h3,    strcat(folderPath, '/beadAggregation.tif'), '-dtiff')
writetable(table(aggregationLevelAll), strcat(folderPath, '/beadAggregation.csv'))

%% ---- Mean segmentation area boxplot -----------------------------------
meanSegAreaUm2AllCell = {meanSegAreaUm2All};

h4 = figure(4); clf
boxPlotCellArray(meanSegAreaUm2AllCell, {'All movies'}, 1, 0, 1);
ylabel('Mean segmentation area (\mum^2)')
title('Mean mask segmentation area')

hgexport(h4, strcat(folderPath, '/beadSegArea'), hgexport('factorystyle'), 'Format', 'eps')
hgsave(h4,   strcat(folderPath, '/beadSegArea'), '-v7.3')
print(h4,    strcat(folderPath, '/beadSegArea.tif'), '-dtiff')
writetable(table(meanSegAreaUm2All), strcat(folderPath, '/beadSegArea.csv'))

%% ---- Combined summary table -------------------------------------------
movieNames = cellfun(@(f) fileparts(f), ML.movieDataFile_, 'UniformOutput', false);
summaryTable = table(movieNames(:), beadDensityAll, aggregationLevelAll, meanSegAreaUm2All, ...
    'VariableNames', {'Movie', 'BeadDensity_per_um2', 'AggregationLevel', 'MeanSegArea_um2'});
writetable(summaryTable, strcat(folderPath, '/beadSummary.csv'))
disp('Done. Results saved to:')
disp(folderPath)
end