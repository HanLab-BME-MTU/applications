%% runFilopodiaForce_debug.m
% Interactive driver to run and debug Process 1 (segmentation) and Process 2
% (detection: tip/base/shaft) on one MovieData, with QC overlays for tuning.
% Run cell by cell (Ctrl+Enter). Edit the CONFIG block first.
% Sangyoon J. Han / 2026

%% ===================== CONFIG =====================
mdPath   = '/mnt/nas/Collaborations/Celine_San_Diego/20260411/Soft-KO/KO_HELA_TALIN_G_BEAD_R_soft_40_nm_P_03_Airyscan_Processing_Processed_2_P3/KO HELA TALIN-G BEAD-R-soft 40 nm_P-03-Airyscan Processing-Processed 2.mat';            % e.g. '/path/to/Cell04/movieData.mat'; '' -> use MD in workspace
iChanTal = 1;             % talin-GFP channel index
tFrame   = 1;             % frame to inspect in QC

if isempty(mdPath)
    assert(exist('MD','var')==1, 'Load MD first or set mdPath.');
else
    MD = MovieData.load(mdPath);
    MD.sanityCheck;
end
fprintf('Movie: %d channels, %d frames, pixelSize=%.4g nm, dt=%.4g s\n', ...
    numel(MD.channels_), MD.nFrames_, MD.pixelSize_, MD.timeInterval_);

%% ===================== PROCESS 1: SEGMENTATION =====================
iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
if isempty(iP1)
    MD.addProcess(FilopodiaSegmentationProcess(MD));
    iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
end
segProc = MD.processes_{iP1};

pp = segProc.funParams_;
pp.ChannelIndex      = iChanTal;
pp.SteerableOrder    = 4;          % even order -> ridge detector
pp.SigmaArray        = [1 2];      % px; tune to talin PSF
pp.BodyThreshold     = 'otsu';    % 'rosin' | 'otsu' | numeric
pp.BodyClosingRadius = 12;         % px, smooths body boundary
pp.ProcessFrames     = tFrame;     % set [] for all frames
pp.GaussianBlurSigma = 2; 
segProc.setPara(pp);

tic; segProc.run(); toc

%% --- P1 QC: body mask + ridge response ---
img      = double(MD.channels_(iChanTal).loadImage(tFrame));
bodyMask = segProc.loadChannelOutput(iChanTal, tFrame, 'output','bodyMask');
res      = segProc.loadChannelOutput(iChanTal, tFrame, 'output','res');

figure('Name','P1 QC','Color','w');
subplot(1,3,1); imshow(img,[100 500]); axis image off; colormap(gca,gray);
hold on; visboundaries(bodyMask,'Color','c','LineWidth',0.5); hold off;
title('talin + body boundary');
subplot(1,3,2); imshow(res,[]); colormap(gca,hot);
title('steerable response (res)');
subplot(1,3,3); imagesc(res .* ~bodyMask); axis image off; colormap(gca,hot);
title('res outside body (shaft signal)');

%% ===================== PROCESS 2: DETECTION =====================
iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
if isempty(iP2)
    MD.addProcess(FilopodiaDetectionProcess(MD));
    iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
end
detProc = MD.processes_{iP2};

pp2 = detProc.funParams_;
pp2.ChannelIndex = iChanTal;
% PSF sigma: use channel value if set, otherwise a sensible Airyscan estimate
if ~isempty(MD.channels_(iChanTal).psfSigma_)
    pp2.PSFsigma = MD.channels_(iChanTal).psfSigma_;
else
    pp2.PSFsigma = 1.6;     % TUNE: ~(PSF FWHM in px)/2.355
end
pp2.Alpha             = 0.05;  % TUNE: lower -> fewer detections
pp2.TipMaxDistFromBody = 70;  % TUNE: max filopodium reach (px)
pp2.MaxTipBaseDist     = 90;  % keep ~1.2 x TipMaxDistFromBody
pp2.OrientLambda       = 3;    % TUNE: higher -> straighter shafts
pp2.OrientTolerance    = 30;
pp2.TipRidgeBand       = 2;    % TUNE: px around shaftMask for ridge gate
pp2.ShaftAbsorbRadius  = 3;
pp2.CarveDistalFrac    = 0.8;
pp2.ProcessFrames      = tFrame;
detProc.setPara(pp2);

tic; detProc.run(); toc

%% --- P2 QC: tip / base / shaft overlay ---
filoInfo = detProc.loadChannelOutput(iChanTal, 'output','filoInfo');
fi = filoInfo{tFrame};
fprintf('Frame %d: %d filopodia detected\n', tFrame, numel(fi));

img = double(MD.channels_(iChanTal).loadImage(tFrame));
figure('Name','P2 QC','Color','w');
imagesc(img); axis image off; colormap(gray); hold on;
for k = 1:numel(fi)
    cl = fi(k).centerline;
    plot(cl(:,1),cl(:,2),'y-','LineWidth',1.2);
    plot(fi(k).basePos(1),fi(k).basePos(2),'co','MarkerSize',7,'LineWidth',1.5);
    plot(fi(k).tipPos(1), fi(k).tipPos(2), 'r+','MarkerSize',9,'LineWidth',1.5);
end
title(sprintf('tip(+) base(o) shaft(-) | n=%d', numel(fi)));

% Length distribution (px -> nm)
if ~isempty(fi)
    L = [fi.length] * MD.pixelSize_;
    figure('Name','P2 lengths','Color','w');
    histogram(L); xlabel('filopodium length (nm)'); ylabel('count');
    title(sprintf('median = %.0f nm', median(L)));
end

%% ===================== QC: WHY are tips detected? (ridge gate audit) =====================
% Diagnose whether false tips pass the ridge gate because shaftMask is too
% broad (background texture included) or because the gate itself is ineffective.
segDir = segProc.funParams_.OutputDirectory;
seg    = load(fullfile(segDir, sprintf('filoSeg_frame_%04d.mat', tFrame)));

res   = seg.res;  shaft = seg.shaftMask > 0;  body = seg.bodyMask > 0;
[H,W] = size(res);
p     = detProc.funParams_;

% Replicate the ridge gate used in detectOneFrame
ridgeDil = imdilate(shaft, strel('disk', round(p.TipRidgeBand)));

% Evaluate each detected tip
nFi = numel(fi);
tipRes     = zeros(1, nFi);
tipOnRidge = false(1, nFi);
tipInt     = zeros(1, nFi);
for k = 1:nFi
    cx = min(max(round(fi(k).tipPos(1)),1),W);
    cy = min(max(round(fi(k).tipPos(2)),1),H);
    tipRes(k)     = res(cy, cx);
    tipOnRidge(k) = ridgeDil(cy, cx);
    tipInt(k)     = img(cy, cx);
end

% Figure A: color tips by on-ridge (green) vs off-ridge (red) over shaftMask
figure('Name','QC: ridge gate audit','Color','w','Position',[50 50 1400 700]);
subplot(1,2,1);
imagesc(img); axis image off; colormap(gca,gray); hold on;
hS = imagesc(cat(3, zeros(H,W), 0.4*ones(H,W), ones(H,W)));   % blue tint
set(hS, 'AlphaData', 0.25 * double(ridgeDil));
visboundaries(body,'Color','c','LineWidth',0.5);
for k = 1:nFi
    cl = fi(k).centerline;
    plot(cl(:,1),cl(:,2),'-','Color',[1 1 0 0.5],'LineWidth',0.8);
    col = 'rg'; col = col(tipOnRidge(k)+1);
    plot(fi(k).tipPos(1),fi(k).tipPos(2),'+','Color',col,'MarkerSize',9,'LineWidth',1.5);
end
title(sprintf('green = on ridge (%d)   red = off ridge (%d)   blue = shaftMask+band', ...
    sum(tipOnRidge), sum(~tipOnRidge)));

% Figure A right: shaftMask outline on raw -- is the mask capturing real ridges?
subplot(1,2,2);
imagesc(img); axis image off; colormap(gca,gray); hold on;
visboundaries(shaft, 'Color',[0.2 0.6 1],'LineWidth',0.4);
visboundaries(body,  'Color','c','LineWidth',0.5);
for k = 1:nFi
    plot(fi(k).tipPos(1),fi(k).tipPos(2),'r+','MarkerSize',8,'LineWidth',1.2);
end
title('shaftMask outline (blue) on raw image -- are these real ridges?');

% Figure B: tip-property distributions to understand what slipped through
figure('Name','QC: tip properties','Color','w','Position',[60 60 1200 380]);
subplot(1,3,1); histogram(tipRes); xlabel('res at tip'); ylabel('count');
title(sprintf('steerable res at tip (median %.1f)', median(tipRes)));
subplot(1,3,2); histogram(tipInt); xlabel('raw talin intensity at tip');
title('raw talin intensity at tip');
subplot(1,3,3);
bar([sum(tipOnRidge) sum(~tipOnRidge)]);
set(gca,'XTickLabel',{'on ridge','off ridge'}); ylabel('count');
title('ridge gate outcome');

% Console summary
fprintf('\n=== tip audit frame %d ===\n', tFrame);
fprintf('total tips: %d  |  on-ridge: %d  |  off-ridge: %d\n', ...
    nFi, sum(tipOnRidge), sum(~tipOnRidge));
fprintf('res at tip:  min %.1f   median %.1f   max %.1f\n', ...
    min(tipRes), median(tipRes), max(tipRes));
fprintf('shaftMask coverage: %.1f%%   after +%dpx band: %.1f%%\n', ...
    100*mean(shaft(:)), round(p.TipRidgeBand), 100*mean(ridgeDil(:)));
fprintf('\nInterpretation:\n');
fprintf('  If most false tips are GREEN  -> shaftMask too broad (P1 hysteresis too loose)\n');
fprintf('  If most false tips are RED    -> ridge gate not applied correctly (bug)\n');
fprintf('  If shaftMask coverage > 15%%  -> background texture captured; tighten HysteresisHigh\n');
fprintf('  If res at false tips is low   -> add a per-tip res floor (TipMinRes param)\n');

%% ===================== TUNING NOTES =====================
% --- P1 (shaftMask is too broad) ---
%   pp.HysteresisHigh = 25;    % raise to tighten ridge mask (default auto)
%   pp.HysteresisLow  = 10;
%   segProc.setPara(pp); segProc.run();
%
% --- P2 (tip sensitivity) ---
%   pp2.Alpha         = 0.01;  % fewer detections
%   pp2.PSFsigma      = 2.0;   % larger PSF -> only bigger spots
%   pp2.TipRidgeBand  = 1;     % stricter: must be almost exactly on ridge
%   pp2.OrientLambda  = 4;     % straighter shafts
%   detProc.setPara(pp2); detProc.run();
