function computeMovieFilopodiaStats(movieData)
%COMPUTEMOVIEFILOPODIASTATS  Process 6. Aggregate filopodium dynamics (P4)
%and force/intensity (P5) into per-movie statistics and correlations.
%
% Produces, per filopodium and pooled over the movie:
%   - geometry/dynamics: length L, net & absolute elongation, protrusion vs
%     retraction speed, lifetime, persistence,
%   - mechanics: tip traction (magnitude & axial), base traction, talin at
%     tip/base, mean shaft traction & talin,
%   - correlations: tip talin vs tip axial traction; protrusion velocity vs
%     tip axial traction; length vs tip traction.
% Output is a stats struct (per-filopodium table + pooled distributions +
% correlation coefficients) for downstream WT/KO group comparison.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaStatisticsProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaStatisticsProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

iCls = getfielddef(p,'ClassProcessIndex',[]);
if isempty(iCls), iCls = movieData.getProcessIndex('FilopodiaClassificationProcess',1,0); end
iSmp = getfielddef(p,'SampleProcessIndex',[]);
if isempty(iSmp), iSmp = movieData.getProcessIndex('FilopodiaSamplingProcess',1,0); end
assert(~isempty(iCls), 'Run FilopodiaClassificationProcess (P4) first.');
clsProc = movieData.processes_{iCls};
clsFile = clsProc.outFilePaths_{1,iChan};
if isempty(clsFile)||exist(clsFile,'file')~=2, clsFile = fullfile(clsProc.funParams_.OutputDirectory,'filoClassification.mat'); end
Scls = load(clsFile, 'filopodia');
filopodia = Scls.filopodia;

haveSmp = ~isempty(iSmp);
filoSamples = [];
if haveSmp
    smpProc = movieData.processes_{iSmp};
    smpFile = smpProc.outFilePaths_{1,iChan};
    if isempty(smpFile)||exist(smpFile,'file')~=2, smpFile = fullfile(smpProc.funParams_.OutputDirectory,'filoSamples.mat'); end
    if exist(smpFile,'file')==2
        Ssmp = load(smpFile, 'filoSamples'); filoSamples = Ssmp.filoSamples;
    else
        haveSmp = false;
    end
end

%% I/O
inFilePaths = cell(1, numel(movieData.channels_)); inFilePaths{1,iChan} = clsFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoStats.mat');
outFilePaths = cell(1, numel(movieData.channels_)); outFilePaths{1,iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

dt  = movieData.timeInterval_; if isempty(dt), dt = 1; end
minLifeStat = getfielddef(p,'MinLifetimeForStats',3);

%% per-filopodium stats
nFil = numel(filopodia);
F = struct('tipTrackId',{},'lifetime',{},'lifetime_s',{}, ...
    'Lmean_nm',{},'Lmax_nm',{},'netElong_nm',{},'absElong_nm',{}, ...
    'vProt_nmps',{},'vRet_nmps',{},'vMeanAbs_nmps',{},'persistence',{}, ...
    'tipForceMean',{},'tipForceMax',{},'tipAxialMean',{},'tipLateralMean',{}, ...
    'baseForceMean',{},'tipTalinMean',{},'baseTalinMean',{}, ...
    'shaftForceMean',{},'shaftTalinMean',{});

for f = 1:nFil
    fil = filopodia(f);
    if fil.lifetime < minLifeStat, continue; end
    L = fil.L_nm(:);
    v = fil.velocity_nmps(:);            % dL/dt (nm/s)
    k = numel(F)+1;
    F(k).tipTrackId   = fil.tipTrackId;
    F(k).lifetime     = fil.lifetime;
    F(k).lifetime_s   = fil.lifetime * dt;
    F(k).Lmean_nm     = mean(L,'omitnan');
    F(k).Lmax_nm      = max(L,[],'omitnan');
    F(k).netElong_nm  = L(end) - L(1);
    F(k).absElong_nm  = sum(abs(diff(L)),'omitnan');
    F(k).vProt_nmps   = mean(v(v>0),'omitnan');           % elongation speed
    F(k).vRet_nmps    = mean(v(v<0),'omitnan');           % retraction speed (negative)
    F(k).vMeanAbs_nmps= mean(abs(v),'omitnan');
    % persistence: net / total path of the length trajectory
    tot = sum(abs(diff(L)),'omitnan');
    F(k).persistence  = abs(L(end)-L(1)) / max(tot, eps);

    % mechanics from P5 (match by tipTrackId)
    [tF,tFmax,tA,tLat,bF,tI,bI,sF,sI] = deal(NaN);
    if haveSmp
        s = matchSample(filoSamples, fil.tipTrackId);
        if ~isempty(s)
            tF  = mean(s.tipForce,'omitnan');
            tFmax = max(s.tipForce,[],'omitnan');
            tA  = mean(s.tipForceAxial,'omitnan');
            tLat= mean(s.tipForceLateral,'omitnan');
            bF  = mean(s.baseForce,'omitnan');
            tI  = mean(s.tipTalin,'omitnan');
            bI  = mean(s.baseTalin,'omitnan');
            [sF, sI] = shaftMeans(s);
        end
    end
    F(k).tipForceMean   = tF;   F(k).tipForceMax = tFmax;
    F(k).tipAxialMean   = tA;   F(k).tipLateralMean = tLat;
    F(k).baseForceMean  = bF;
    F(k).tipTalinMean   = tI;   F(k).baseTalinMean = bI;
    F(k).shaftForceMean = sF;   F(k).shaftTalinMean = sI;
end

%% pooled distributions + correlations
pooled.n            = numel(F);
pooled.Lmean_nm     = [F.Lmean_nm];
pooled.lifetime_s   = [F.lifetime_s];
pooled.vProt_nmps   = [F.vProt_nmps];
pooled.vRet_nmps    = [F.vRet_nmps];
pooled.netElong_nm  = [F.netElong_nm];
pooled.persistence  = [F.persistence];
pooled.tipForceMean = [F.tipForceMean];
pooled.tipAxialMean = [F.tipAxialMean];
pooled.tipTalinMean = [F.tipTalinMean];

corrStats = struct();
corrStats.talin_vs_axial   = safeCorr([F.tipTalinMean], [F.tipAxialMean]);
corrStats.vProt_vs_axial   = safeCorr([F.vProt_nmps],   [F.tipAxialMean]);
corrStats.L_vs_tipForce    = safeCorr([F.Lmean_nm],     [F.tipForceMean]);
corrStats.talin_vs_force   = safeCorr([F.tipTalinMean], [F.tipForceMean]);

stats.perFilo   = F;
stats.pooled    = pooled;
stats.corr      = corrStats;
stats.nFilopodia= numel(F);
stats.movieName = movieData.movieDataFileName_;

%% length-normalized shaft profiles (per-filopodium time-average), aggregated over the movie
normProf = struct('sN',[],'force',[],'forceAxial',[],'talin',[], ...
    'forceMean',[],'forceAxialMean',[],'talinMean',[], ...
    'forceSEM',[],'forceAxialSEM',[],'talinSEM',[],'nFilo',0);
if haveSmp && ~isempty(filoSamples) && isfield(filoSamples,'normForceMean')
    keepIds = [F.tipTrackId];                 % only filopodia passing MinLifetimeForStats
    sN = filoSamples(1).normS;
    Pf=[]; Pa=[]; Pt=[];
    for q = 1:numel(filoSamples)
        if ~ismember(filoSamples(q).tipTrackId, keepIds), continue; end
        Pf = [Pf, filoSamples(q).normForceMean(:)];      %#ok<AGROW>
        Pa = [Pa, filoSamples(q).normForceAxialMean(:)]; %#ok<AGROW>
        Pt = [Pt, filoSamples(q).normTalinMean(:)];      %#ok<AGROW>
    end
    normProf.sN             = sN;
    normProf.force          = Pf;     % nNorm x nFilo
    normProf.forceAxial     = Pa;
    normProf.talin          = Pt;
    normProf.nFilo          = size(Pf,2);
    normProf.forceMean      = mean(Pf,2,'omitnan');
    normProf.forceAxialMean = mean(Pa,2,'omitnan');
    normProf.talinMean      = mean(Pt,2,'omitnan');
    normProf.forceSEM       = std(Pf,0,2,'omitnan')./sqrt(max(sum(isfinite(Pf),2),1));
    normProf.forceAxialSEM  = std(Pa,0,2,'omitnan')./sqrt(max(sum(isfinite(Pa),2),1));
    normProf.talinSEM       = std(Pt,0,2,'omitnan')./sqrt(max(sum(isfinite(Pt),2),1));
end
stats.normProfile = normProf;

save(outFile, 'stats', '-v7.3');
fprintf(['Stats: %d filopodia | median L %.0f nm | median lifetime %.1f s | ' ...
    'median tip traction %.1f Pa\n'], numel(F), median(pooled.Lmean_nm,'omitnan'), ...
    median(pooled.lifetime_s,'omitnan'), median(pooled.tipForceMean,'omitnan'));
fprintf('corr: talin~axial r=%.2f | vProt~axial r=%.2f | talin~force r=%.2f\n', ...
    corrStats.talin_vs_axial, corrStats.vProt_vs_axial, corrStats.talin_vs_force);

%% figure: length-normalized shaft profiles (individual thin + bold mean)
makeFig = getfielddef(p,'MakeFigures',true);
if makeFig && normProf.nFilo > 0
    figDir = fullfile(outDir,'Figs'); if exist(figDir,'dir')~=7, mkdir(figDir); end
    drawNormProfile(normProf.sN, normProf.force,      normProf.forceMean, ...
        'traction (Pa)', 'Normalized shaft profile: traction', fullfile(figDir,'normProfile_traction'));
    drawNormProfile(normProf.sN, normProf.forceAxial, normProf.forceAxialMean, ...
        'axial traction (Pa, + outward)', 'Normalized shaft profile: axial traction', fullfile(figDir,'normProfile_axialTraction'));
    drawNormProfile(normProf.sN, normProf.talin,      normProf.talinMean, ...
        'talin intensity (a.u.)', 'Normalized shaft profile: talin', fullfile(figDir,'normProfile_talin'));
end
end

% ===================================================================
function drawNormProfile(sN, M, mu, ylab, ttl, fbase)
% M: nNorm x nFilo individual profiles; mu: nNorm x 1 mean
h = figure('Color','w','Visible','off');
hold on;
plot(sN, M, '-', 'Color',[0.7 0.7 0.85], 'LineWidth',0.5);
plot(sN, mu, '-', 'Color',[0.10 0.30 0.80], 'LineWidth',2.5);
hold off; box on;
xlabel('normalized arc length  (0 = tip \rightarrow 1 = base)');
ylabel(ylab); title(sprintf('%s  (n=%d)', ttl, size(M,2)));
xlim([0 1]);
try
    hgexport(h, fbase, hgexport('factorystyle'), 'Format','eps');
    hgsave(h, fbase, '-v7.3');
catch
end
print(h, [fbase '.tif'], '-dtiff');
close(h);
end

% ===================================================================
function s = matchSample(filoSamples, tipTrackId)
s = [];
ids = [filoSamples.tipTrackId];
k = find(ids == tipTrackId, 1);
if ~isempty(k), s = filoSamples(k); end
end

% ===================================================================
function [mF, mI] = shaftMeans(s)
mF = NaN; mI = NaN;
if ~isfield(s,'shaftProfile') || isempty(s.shaftProfile), return; end
allF = []; allI = [];
for j = 1:numel(s.shaftProfile)
    sp = s.shaftProfile{j};
    if isempty(sp), continue; end
    allF = [allF; sp.force(:)];  allI = [allI; sp.talin(:)]; %#ok<AGROW>
end
mF = mean(allF,'omitnan'); mI = mean(allI,'omitnan');
end

% ===================================================================
function r = safeCorr(x, y)
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 3, r = NaN; return; end
c = corrcoef(x(ok), y(ok));
r = c(1,2);
end

% ===================================================================
function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
