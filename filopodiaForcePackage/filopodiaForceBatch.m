%% filopodiaForceBatch
% Batch comparison of FilopodiaForcePackage results across conditions
% (e.g. WT, ArhGAP39KO, 1YA).
%
% Loads FilopodiaStatisticsProcess (P6) output (filoStats.mat) from every
% movie in every condition's MovieList, then compares between conditions:
%   - filopodia count per cell (per-movie scalar)
%   - length, lifetime, protrusion/retraction speed, persistence (per-filopodium pooled)
%   - traction at tip / shaft / base, talin at tip / base (per-filopodium pooled)
%   - tip axial traction, length-normalized shaft profiles (traction/axial/talin)
%
% Mirrors strainEnergyBatch layout: boxPlotCellArray + CSV + EPS/fig/tif.
% Stats output is raw data + boxplots only (no significance test).
% Sangyoon J. Han / 2026

%% open necessary MLs
[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat, specificName] = chooseSelectedFolders;
numConditions = numel(pathAnalysisAll);
for k = 1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]); %#ok<SAGROW>
end

%% Output
rootAnalysis = pwd;
figPath  = [rootAnalysis '/AnalysisSummary_Filopodia' specificName '/Figs'];  mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary_Filopodia' specificName '/Data'];  mkdir(dataPath)

%% group containers (one cell per condition)
% per-movie scalars
count_Group        = cell(numConditions,1);   % # filopodia per movie

% per-filopodium pooled (cell per condition -> {Nmovies x 1} -> vector per movie)
L_Group            = cell(numConditions,1);    % mean length (nm)
Lmax_Group         = cell(numConditions,1);
lifetime_Group     = cell(numConditions,1);    % s
vProt_Group        = cell(numConditions,1);    % nm/s
vRet_Group         = cell(numConditions,1);    % nm/s (negative)
persistence_Group  = cell(numConditions,1);
netElong_Group     = cell(numConditions,1);

tipForce_Group     = cell(numConditions,1);    % Pa
tipAxial_Group     = cell(numConditions,1);    % Pa
shaftForce_Group   = cell(numConditions,1);    % Pa
baseForce_Group    = cell(numConditions,1);    % Pa
tipTalin_Group     = cell(numConditions,1);
baseTalin_Group    = cell(numConditions,1);
shaftTalin_Group   = cell(numConditions,1);

% length-normalized shaft profiles: per condition, accumulate columns (one per filopodium)
normForce_Group     = cell(numConditions,1);
normAxial_Group     = cell(numConditions,1);
normTalin_Group     = cell(numConditions,1);
sN_ref = [];   % common normalized arc-length grid

N = zeros(numConditions,1);

for ii = 1:numConditions
    curML = MLAll(ii);
    curMovies = curML.movies_;
    N(ii) = numel(curMovies);

    curCount = nan(N(ii),1);
    [curL,curLmax,curLife,curVp,curVr,curPers,curNet] = deal(cell(N(ii),1));
    [curTipF,curTipA,curShaftF,curBaseF,curTipI,curBaseI,curShaftI] = deal(cell(N(ii),1));
    curNF = []; curNA = []; curNT = [];

    for k = 1:N(ii)
        curMovie = curMovies{k};
        iProc = curMovie.getProcessIndex('FilopodiaStatisticsProcess',1,0);
        if isempty(iProc), warning('Movie %d in cond %d has no FilopodiaStatisticsProcess.',k,ii); continue; end
        proc = curMovie.getProcess(iProc);

        sFile = proc.outFilePaths_{1,1};
        if isempty(sFile)||exist(sFile,'file')~=2
            sFile = fullfile(proc.funParams_.OutputDirectory,'filoStats.mat');
        end
        if exist(sFile,'file')~=2, warning('No filoStats.mat for movie %d cond %d.',k,ii); continue; end
        S = load(sFile,'stats'); st = S.stats;

        curCount(k)   = st.nFilopodia;
        curL{k}       = st.pooled.Lmean_nm(:);
        curLife{k}    = st.pooled.lifetime_s(:);
        curVp{k}      = st.pooled.vProt_nmps(:);
        curVr{k}      = st.pooled.vRet_nmps(:);
        curPers{k}    = st.pooled.persistence(:);
        curNet{k}     = st.pooled.netElong_nm(:);
        curTipF{k}    = st.pooled.tipForceMean(:);
        curTipA{k}    = st.pooled.tipAxialMean(:);
        curTipI{k}    = st.pooled.tipTalinMean(:);

        % per-filopodium fields not in pooled: pull from perFilo
        pf = st.perFilo;
        curLmax{k}   = [pf.Lmax_nm]';
        curShaftF{k} = [pf.shaftForceMean]';
        curBaseF{k}  = [pf.baseForceMean]';
        curBaseI{k}  = [pf.baseTalinMean]';
        curShaftI{k} = [pf.shaftTalinMean]';

        % normalized profiles
        if isfield(st,'normProfile') && st.normProfile.nFilo > 0
            np = st.normProfile;
            if isempty(sN_ref), sN_ref = np.sN(:); end
            curNF = [curNF, np.force];       %#ok<AGROW>
            curNA = [curNA, np.forceAxial];  %#ok<AGROW>
            curNT = [curNT, np.talin];       %#ok<AGROW>
        end
    end

    count_Group{ii}       = curCount;
    L_Group{ii}           = curL;
    Lmax_Group{ii}        = curLmax;
    lifetime_Group{ii}    = curLife;
    vProt_Group{ii}       = curVp;
    vRet_Group{ii}        = curVr;
    persistence_Group{ii} = curPers;
    netElong_Group{ii}    = curNet;
    tipForce_Group{ii}    = curTipF;
    tipAxial_Group{ii}    = curTipA;
    shaftForce_Group{ii}  = curShaftF;
    baseForce_Group{ii}   = curBaseF;
    tipTalin_Group{ii}    = curTipI;
    baseTalin_Group{ii}   = curBaseI;
    shaftTalin_Group{ii}  = curShaftI;
    normForce_Group{ii}   = curNF;
    normAxial_Group{ii}   = curNA;
    normTalin_Group{ii}   = curNT;
end
disp('Loading done')

%% group names (same logic as strainEnergyBatch)
if usedSelectedFoldersMat
    groupNames2 = groupNames;
    for ii = 1:numConditions
        [~, finalFolder] = fileparts(pathAnalysisAll{ii});
        groupNames{ii} = finalFolder;
    end
    nameList = groupNames';
    if any(cellfun(@isempty,nameList))
        nameList = MLNames;
        if numel(nameList)>1 && strcmp(nameList{1},nameList{2})
            for ii = 1:numConditions
                curPath = fileparts(pathAnalysisAll{ii});
                [~,nameList{ii}] = fileparts(curPath);
            end
        end
    end
    [uniqueNames, ~, idx] = unique(nameList,'stable'); %#ok<ASGLU>
    counts = accumarray(idx(:),1);
    if any(counts>1), nameList = groupNames2'; end
else
    nameList = groupNames';
end

%% ---- helper for pooled per-filopodium box plots ----
plotPooled = @(G,ylab,ttl,fname) localBoxPooled(G,nameList,ylab,ttl,figPath,dataPath,fname);

%% counts per cell (per-movie scalar)
count_CellArray = cellfun(@(x) x(:), count_Group, 'unif', false);
h1 = figure; boxPlotCellArray(count_CellArray,nameList,1,false,true);
ylabel('# filopodia per cell'); title('Filopodia count per cell');
saveAll(h1,figPath,'filopodiaCount');
writetable(table(count_CellArray,'RowNames',nameList), [dataPath '/filopodiaCount.csv'],'WriteRowNames',true);

%% geometry / dynamics
plotPooled(L_Group,           'Mean length (nm)',            'Filopodia length',               'length');
plotPooled(Lmax_Group,        'Max length (nm)',             'Filopodia max length',           'lengthMax');
plotPooled(lifetime_Group,    'Lifetime (s)',                'Filopodia lifetime',             'lifetime');
plotPooled(vProt_Group,       'Protrusion speed (nm/s)',     'Protrusion speed',               'vProt');
plotPooled(vRet_Group,        'Retraction speed (nm/s)',     'Retraction speed',               'vRet');
plotPooled(persistence_Group, 'Persistence (net/total)',     'Length-change persistence',      'persistence');
plotPooled(netElong_Group,    'Net elongation (nm)',         'Net elongation',                 'netElong');

%% mechanics: traction
plotPooled(tipForce_Group,   'Tip traction (Pa)',          'Traction at filopodium tip',     'tipForce');
plotPooled(shaftForce_Group, 'Shaft traction (Pa)',        'Traction along filopodium shaft','shaftForce');
plotPooled(baseForce_Group,  'Base traction (Pa)',         'Traction at filopodium base',    'baseForce');
plotPooled(tipAxial_Group,   'Tip axial traction (Pa)',    'Axial traction at tip (+outward)','tipAxial');

%% mechanics: talin
plotPooled(tipTalin_Group,   'Tip talin (a.u.)',           'Talin at filopodium tip',        'tipTalin');
plotPooled(shaftTalin_Group, 'Shaft talin (a.u.)',         'Talin along filopodium shaft',   'shaftTalin');
plotPooled(baseTalin_Group,  'Base talin (a.u.)',          'Talin at filopodium base',       'baseTalin');

%% length-normalized shaft profile comparison (condition means overlaid)
compareNormProfile(sN_ref, normForce_Group, nameList, 'traction (Pa)', ...
    'Normalized shaft traction profile', figPath, dataPath, 'normProfile_traction');
compareNormProfile(sN_ref, normAxial_Group, nameList, 'axial traction (Pa, + outward)', ...
    'Normalized shaft axial-traction profile', figPath, dataPath, 'normProfile_axialTraction');
compareNormProfile(sN_ref, normTalin_Group, nameList, 'talin (a.u.)', ...
    'Normalized shaft talin profile', figPath, dataPath, 'normProfile_talin');

%% save workspace
save([dataPath filesep 'allData.mat'])
disp('filopodiaForceBatch done')

% =====================================================================
function localBoxPooled(G, nameList, ylab, ttl, figPath, dataPath, fname)
% G: cell{numCond} -> {Nmovies x 1} cells (or numeric per-movie) of per-filo vectors
CellArray = cellfun(@(x) poolCond(x), G, 'unif', false);
h1 = figure;
try
    boxPlotCellArray(CellArray,nameList,1,false,true);
catch
    boxPlotCellArray(CellArray,nameList',1,false,true);
end
ylabel(ylab); title(ttl);
saveAll(h1,figPath,fname);
writetable(table(CellArray,'RowNames',nameList), [dataPath '/' fname '.csv'],'WriteRowNames',true);
end

function v = poolCond(x)
% concatenate per-movie vectors of one condition into a single pooled vector
if iscell(x)
    v = [];
    for k = 1:numel(x), v = [v; x{k}(:)]; end %#ok<AGROW>
else
    v = x(:);
end
end

% =====================================================================
function compareNormProfile(sN, Gmat, nameList, ylab, ttl, figPath, dataPath, fname)
% Gmat: cell{numCond} -> nNorm x nFilo matrix of normalized profiles
if isempty(sN), warning('No normalized profiles found; skipping %s.',fname); return; end
nc = numel(Gmat);
cmap = lines(nc);
h1 = figure('Color','w'); hold on;
legendEntries = {};
meanTable = sN(:);  varNames = {'normS'};
for ii = 1:nc
    M = Gmat{ii};
    if isempty(M), continue; end
    mu  = mean(M,2,'omitnan');
    sem = std(M,0,2,'omitnan')./sqrt(max(sum(isfinite(M),2),1));
    % shaded SEM band
    xb = [sN(:); flipud(sN(:))];
    yb = [mu-sem; flipud(mu+sem)];
    ok = isfinite(xb)&isfinite(yb);
    if any(ok)
        fill(xb(ok), yb(ok), cmap(ii,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    end
    plot(sN, mu, '-', 'Color',cmap(ii,:), 'LineWidth',2.5);
    legendEntries{end+1} = sprintf('%s (n=%d)', nameList{ii}, size(M,2)); %#ok<AGROW>
    meanTable = [meanTable, mu, sem]; %#ok<AGROW>
    varNames = [varNames, {[matlab.lang.makeValidName(nameList{ii}) '_mean'], ...
                           [matlab.lang.makeValidName(nameList{ii}) '_sem']}]; %#ok<AGROW>
end
hold off; box on;
xlabel('normalized arc length  (0 = tip \rightarrow 1 = base)');
ylabel(ylab); title(ttl); xlim([0 1]);
legend(legendEntries,'Location','best');
saveAll(h1,figPath,fname);
T = array2table(meanTable,'VariableNames',varNames);
writetable(T, [dataPath '/' fname '.csv']);
end

% =====================================================================
function saveAll(h, figPath, fname)
try
    hgexport(h, [figPath '/' fname], hgexport('factorystyle'), 'Format','eps');
    hgsave(h, [figPath '/' fname], '-v7.3');
catch
end
print(h, [figPath '/' fname '.tif'], '-dtiff');
end
