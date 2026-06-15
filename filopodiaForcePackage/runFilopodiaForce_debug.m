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
pp2.TipMaxDistFromBody = 60;     % TUNE: max filopodium reach (px) from body
% pp2.BaseSearchBand     = 5;      % TUNE: px band at body edge counted as base
pp2.BaseSearchBand   = 12;       % px, base expansion
pp2.MaxTipBaseDist     = 80;    % TUNE: max plausible filo length (px)
pp2.OrientTolerance    = 30;     % TUNE: deg; raise if traces miss bends
pp2.OrientLambda       = 3;      % TUNE: raise to force ridge-following
pp2.MinFiloLength      = 5;      % px
pp2.ProcessFrames      = [];     % set to tFrame for a fast single-frame test
pp2.CarveDistalFrac  = 0.8;   % fraction to block distal portion. 1.0 will block even base

pp2.ShaftAbsorbRadius  = 3;     % shaft middle absortiopn radius (if increased, it will be more aggressive)

pp2.TipScoreWeightDist = 0.5;   % dist vs amplitude even consideration
pp2.TipScoreWeightAmp  = 0.5;
pp2.TipScoreMinPrctile = 50;    % low 50% removed; If increased, it will be more strict
pp2.TipMinAmplitude    = 0;     % absolute amplitude floor (0=off); kills isolated dim noise
pp2.TipRidgeBand = 1;     % If inreased(3~4), more generous. If decreased (to 1), more strict.

detProc.setPara(pp2);

tic; detProc.run(); toc

%% ===================== QC P2: tip/base/shaft overlay =====================
filoInfo = detProc.loadChannelOutput(iChanTal, 'output','filoInfo');
fi = filoInfo{tFrame};
fprintf('Frame %d: %d filopodia detected\n', tFrame, numel(fi));

img = double(MD.channels_(iChanTal).loadImage(tFrame));
figure('Name','P2 QC','Color','w');

subplot(1,2,2); imshow(res.*~bodyMask,[]); axis image off; colormap(gca,hot);

subplot(1,2,1);
imshow(img); axis image off; colormap(gray); hold on;
clim([100 500])

for k = 1:numel(fi)
    cl = fi(k).centerline;
    plot(cl(:,1), cl(:,2), 'y-', 'LineWidth', 1.2);      % shaft trace
    plot(fi(k).basePos(1), fi(k).basePos(2), 'co', 'MarkerSize',7,'LineWidth',1.5); % base
    plot(fi(k).tipPos(1),  fi(k).tipPos(2),  'r+', 'MarkerSize',9,'LineWidth',1.5); % tip
end
title(sprintf('tip(+) base(o) shaft(-) | n=%d', numel(fi)));
hold off;

% length distribution (px -> nm)
% if ~isempty(fi)
%     L = [fi.length] * MD.pixelSize_;
%     figure('Name','P2 lengths','Color','w');
%     histogram(L); xlabel('filopodium length (nm)'); ylabel('count');
%     title(sprintf('median = %.0f nm', median(L)));
% end

%% ===================== diagnostics for tuning =====================
% If too few/many tips: adjust pp2.PSFsigma, pp2.Alpha, pp2.TipMaxDistFromBody.
% If shafts cut corners or wander: raise OrientLambda / lower OrientTolerance.
% If bases land mid-shaft: lower BaseSearchBand. Re-run the PROCESS 2 cell.

%% ===================== QC: WHY are tips detected? (ridge gate audit) =====================
% ? ?? tip? shaftMask(strong ridge) ?? ??? ???? ????.
tFrame = 1;
iChanTal = 1;

% P1, P2 ?? ??
segDir = MD.processes_{MD.getProcessIndex('FilopodiaSegmentationProcess',1,0)}.funParams_.OutputDirectory;
seg = load(fullfile(segDir, sprintf('filoSeg_frame_%04d.mat', tFrame)));
detProc = MD.processes_{MD.getProcessIndex('FilopodiaDetectionProcess',1,0)};
filoInfo = detProc.loadChannelOutput(iChanTal, 'output','filoInfo');
fi = filoInfo{tFrame};

img   = double(MD.channels_(iChanTal).loadImage(tFrame));
res   = seg.res;  shaft = seg.shaftMask>0;  body = seg.bodyMask>0;
[H,W] = size(res);

% ??? ?? ridge ???? ???? ???
p = detProc.funParams_;
ridge = imdilate(shaft, strel('disk', round(p.TipRidgeBand)));

% ? tip??: ridge ???, ? ?? res ?, ? ?? raw ??
nF = numel(fi);
tipRes = zeros(1,nF); tipOnRidge = false(1,nF); tipInt = zeros(1,nF);
for k = 1:nF
    c = round(fi(k).tipPos);                 % [x y]
    cx = min(max(c(1),1),W); cy = min(max(c(2),1),H);
    tipRes(k)     = res(cy,cx);
    tipOnRidge(k) = ridge(cy,cx);
    tipInt(k)     = img(cy,cx);
end

%% Figure 1: tip? ridge ?/??? ? ?? + shaftMask ??
figure('Name','QC tip audit','Color','w','Position',[50 50 1400 700]);

subplot(1,2,1);
imagesc(img); axis image off; colormap(gca,gray); hold on;
% shaftMask? ??? ????
sm = imdilate(shaft, strel('disk',round(p.TipRidgeBand)));
h = imagesc(cat(3, zeros(H,W), 0.4*ones(H,W), ones(H,W)));
set(h,'AlphaData', 0.25*double(sm));
visboundaries(body,'Color','c','LineWidth',0.5);
for k = 1:nF
    cl = fi(k).centerline;
    plot(cl(:,1),cl(:,2),'-','Color',[1 1 0 0.5],'LineWidth',0.8);
    col = 'rg'; col = col(tipOnRidge(k)+1);  % off-ridge=red, on-ridge=green
    plot(fi(k).tipPos(1),fi(k).tipPos(2),'+','Color',col,'MarkerSize',9,'LineWidth',1.5);
end
title(sprintf('green=on ridge (%d), red=off ridge (%d) | blue=shaftMask(+band)', ...
    sum(tipOnRidge), sum(~tipOnRidge)));

% raw ?? shaftMask ??? (ridge? ???? ???)
subplot(1,2,2);
imagesc(img); axis image off; colormap(gca,gray); hold on;
visboundaries(shaft,'Color',[0.2 0.6 1],'LineWidth',0.4);
visboundaries(body,'Color','c','LineWidth',0.5);
plot([fi.tipPos]*0+nan, nan); % noop
for k = 1:nF
    plot(fi(k).tipPos(1),fi(k).tipPos(2),'r+','MarkerSize',8,'LineWidth',1.2);
end
title('shaftMask outline (blue) on raw ? are these real ridges?');

%% Figure 2: tip ?? ?? ?? (? ????)
figure('Name','QC tip properties','Color','w','Position',[60 60 1200 380]);
subplot(1,3,1); histogram(tipRes); xlabel('res at tip'); ylabel('count');
title(sprintf('res at tip (median %.1f)', median(tipRes)));
subplot(1,3,2); histogram(tipInt); xlabel('raw intensity at tip');
title('raw talin intensity at tip');
subplot(1,3,3);
bar([sum(tipOnRidge) sum(~tipOnRidge)]); set(gca,'XTickLabel',{'on ridge','off ridge'});
title('ridge gate outcome'); ylabel('count');

%% ?? ??
fprintf('\n=== tip audit (frame %d) ===\n', tFrame);
fprintf('total tips: %d | on-ridge: %d | off-ridge: %d\n', nF, sum(tipOnRidge), sum(~tipOnRidge));
fprintf('res at tip: min %.1f  median %.1f  max %.1f\n', min(tipRes), median(tipRes), max(tipRes));
fprintf('(shaftMask covers %.1f%% of image; +band %.1f%%)\n', ...
    100*mean(shaft(:)), 100*mean(sm(:)));
%% P3 running
iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
if isempty(iP3), MD.addProcess(FilopodiaTrackingProcess(MD)); 
   iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0); end
trkProc = MD.processes_{iP3};
pp3 = trkProc.funParams_;
pp3.MaxLinkDist   = 10;   % tip? ???? ???? ?? px (protrusion ??? ??)
pp3.MaxGapFrames  = 2;    % ?? ??? ? ??
pp3.MinTrackLength = 3;   % 3??? ?? ?? (false ??)
trkProc.setPara(pp3);
trkProc.run();

%% P3 QC
trkProc = MD.processes_{MD.getProcessIndex('FilopodiaTrackingProcess',1,0)};
ft = trkProc.loadChannelOutput(1,'output','filoTracks');

%% 1) track ??(lifetime) ?? ? 3? ????, ?? ?? ? ??
lifeF = [ft.lifetime];
figure; histogram(lifeF); xlabel('track lifetime (frames)'); ylabel('count');
title(sprintf('median %d, max %d frames | n=%d tracks', median(lifeF), max(lifeF), numel(ft)));

%% 2) ?? overlay ? ???? track? ?? filopodia ????
tF = 1; img = double(MD.channels_(1).loadImage(tF));
figure; imshow(img,[100 500]); hold on;
for i=1:numel(ft)
    j = find(ft(i).frames==tF,1);
    if isempty(j), continue; end
    plot(ft(i).tipTrack(j,1), ft(i).tipTrack(j,2),'r+','MarkerSize',8,'LineWidth',1.2);
end
title(sprintf('tips of tracks present at frame %d', tF));

%% 3) velocity ?? sanity ? ?? ? ?? ?? linking ??
allV = abs([ft.velocity_nmps]);
fprintf('|velocity|: median %.0f nm/s, 95pct %.0f nm/s, max %.0f nm/s\n', ...
    median(allV), prctile(allV,95), max(allV));