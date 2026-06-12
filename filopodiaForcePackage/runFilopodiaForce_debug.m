%% runFilopodiaForce_debug.m
% Interactive driver to run + debug Process 1 (segmentation) and Process 2
% (detection: tip/base/shaft) on one MovieData, with QC overlays for tuning.
% Run cell by cell (Ctrl+Enter). Edit the CONFIG block first.
% Sangyoon J. Han / 2026

%% ===================== CONFIG =====================
% Provide either a path to movieData.mat or an MD already in the workspace.
mdPath   = '/mnt/nas/Collaborations/Celine_San_Diego/20260411/Soft-KO/KO_HELA_TALIN_G_BEAD_R_soft_40_nm_P_03_Airyscan_Processing_Processed_2_P3/KO HELA TALIN-G BEAD-R-soft 40 nm_P-03-Airyscan Processing-Processed 2.mat';            % e.g. '/path/to/Cell04/movieData.mat'; '' -> use MD var
iChanTal = 1;             % talin-GFP channel index
tFrame   = 1;             % frame to inspect in QC

if isempty(mdPath)
    assert(exist('MD','var')==1, 'Load MD or set mdPath.');
else
    % S = load(mdPath); fn = fieldnames(S); MD = S.(fn{1});
    MD = MovieData.load(mdPath);
    % MD.sanityCheck;       % refresh paths if relocated
end
fprintf('Movie: %d channels, %d frames, pixelSize=%g nm, dt=%g s\n', ...
    numel(MD.channels_), MD.nFrames_, MD.pixelSize_, MD.timeInterval_);

%% ===================== PROCESS 1: SEGMENTATION =====================
iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
if isempty(iP1)
    p1 = FilopodiaSegmentationProcess(MD);
    MD.addProcess(p1);
    iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
end
segProc = MD.processes_{iP1};

pp = segProc.funParams_;
pp.ChannelIndex   = iChanTal;
pp.SteerableOrder = 4;            % ridge detector
pp.SigmaArray     = [1 2 4];      % px; lower these if shafts are very thin
pp.BodyThreshold  = 'rosin';      % 'rosin'|'otsu'|numeric
pp.ProcessFrames  = [];           % set to tFrame to test a single frame fast

% pp.SigmaArray        = [1 2];    % sigma=4 deletion -> central blob supression
pp.BodyClosingRadius = 4;       % px, boundary smoothing


segProc.setPara(pp);

tic; segProc.run(); toc

%% ---- QC P1: body mask + ridge response on tFrame ----
img = double(MD.channels_(iChanTal).loadImage(tFrame));
bodyMask = segProc.loadChannelOutput(iChanTal, tFrame, 'output','bodyMask');
res      = segProc.loadChannelOutput(iChanTal, tFrame, 'output','res');

figure('Name','P1 QC','Color','w');
subplot(1,3,1); imagesc(img); axis image off; colormap(gca,gray); title('talin');
hold on; visboundaries(bodyMask,'Color','c','LineWidth',0.5); hold off;
subplot(1,3,2); imagesc(res); axis image off; colormap(gca,hot); title('steerable res');
subplot(1,3,3); imagesc(res.*~bodyMask); axis image off; colormap(gca,hot);
title('res outside body (shaft signal)');

%% ===================== PROCESS 2: DETECTION =====================
iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
if isempty(iP2)
    p2 = FilopodiaDetectionProcess(MD);
    MD.addProcess(p2);
    iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
end
detProc = MD.processes_{iP2};

pp2 = detProc.funParams_;
pp2.ChannelIndex     = iChanTal;
% PSF sigma in px: use channel value if present, else a sensible Airyscan guess
if ~isempty(MD.channels_(iChanTal).psfSigma_)
    pp2.PSFsigma = MD.channels_(iChanTal).psfSigma_;
else
    pp2.PSFsigma = 2.3;          % TUNE: ~ (PSF FWHM in px)/2.355
end
pp2.Alpha             = 0.05;
pp2.TipMaxDistFromBody = 50;     % TUNE: max filopodium reach (px) from body
% pp2.BaseSearchBand     = 5;      % TUNE: px band at body edge counted as base
pp2.BaseSearchBand   = 12;       % px, base expansion
pp2.MaxTipBaseDist     = 100;    % TUNE: max plausible filo length (px)
pp2.OrientTolerance    = 30;     % TUNE: deg; raise if traces miss bends
pp2.OrientLambda       = 2;      % TUNE: raise to force ridge-following
pp2.MinFiloLength      = 5;      % px
pp2.ProcessFrames      = [];     % set to tFrame for a fast single-frame test
detProc.setPara(pp2);

tic; detProc.run(); toc

%% ===================== QC P2: tip/base/shaft overlay =====================
filoInfo = detProc.loadChannelOutput(iChanTal, 'output','filoInfo');
fi = filoInfo{tFrame};
fprintf('Frame %d: %d filopodia detected\n', tFrame, numel(fi));

img = double(MD.channels_(iChanTal).loadImage(tFrame));
figure('Name','P2 QC','Color','w');
imagesc(img); axis image off; colormap(gray); hold on;
clim([100 600])

for k = 1:numel(fi)
    cl = fi(k).centerline;
    plot(cl(:,1), cl(:,2), 'y-', 'LineWidth', 1.2);      % shaft trace
    plot(fi(k).basePos(1), fi(k).basePos(2), 'co', 'MarkerSize',7,'LineWidth',1.5); % base
    plot(fi(k).tipPos(1),  fi(k).tipPos(2),  'r+', 'MarkerSize',9,'LineWidth',1.5); % tip
end
title(sprintf('tip(+) base(o) shaft(-) | n=%d', numel(fi)));
hold off;

% length distribution (px -> nm)
if ~isempty(fi)
    L = [fi.length] * MD.pixelSize_;
    figure('Name','P2 lengths','Color','w');
    histogram(L); xlabel('filopodium length (nm)'); ylabel('count');
    title(sprintf('median = %.0f nm', median(L)));
end

%% ===================== diagnostics for tuning =====================
% If too few/many tips: adjust pp2.PSFsigma, pp2.Alpha, pp2.TipMaxDistFromBody.
% If shafts cut corners or wander: raise OrientLambda / lower OrientTolerance.
% If bases land mid-shaft: lower BaseSearchBand. Re-run the PROCESS 2 cell.
