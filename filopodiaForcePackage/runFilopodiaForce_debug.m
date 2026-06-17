%% runFilopodiaForce_debug.m
% Interactive driver to run and debug Process 1 (segmentation) and Process 2
% (detection) on one MovieData, with QC overlays for tuning.
% Run cell by cell (Ctrl+Enter). Edit the CONFIG block first.
% Sangyoon J. Han / 2026
clear 
clear class
%% ===================== CONFIG =====================
mdPath   = '/mnt/nas/Collaborations/Celine_San_Diego/20260411/Soft-KO/KO_HELA_TALIN_G_BEAD_R_soft_40_nm_P_03_Airyscan_Processing_Processed_2_P3/KO HELA TALIN-G BEAD-R-soft 40 nm_P-03-Airyscan Processing-Processed 2.mat';
iChanTal = 1;             % talin-GFP channel index
tFrame   = [];             % frame to inspect in QC

if isempty(mdPath)
    assert(exist('MD','var')==1, 'Load MD first or set mdPath.');
else
    MD = MovieData.load(mdPath);
    MD.sanityCheck;
end
fprintf('Movie: %d channels, %d frames, pixelSize=%.4g nm, dt=%.4g s\n', ...
    numel(MD.channels_), MD.nFrames_, MD.pixelSize_, MD.timeInterval_);

%% ===================== PACKAGE REGISTRATION =====================
% Register FilopodiaForcePackage on MD if not already present.
% Only needs to run once; persists after MD.save().
pkg = setupFilopodiaForcePackage(MD);

%% ===================== PROCESS 1: SEGMENTATION =====================
iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
if isempty(iP1)
    MD.addProcess(FilopodiaSegmentationProcess(MD));
    iP1 = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
end
segProc = MD.processes_{iP1};

pp = segProc.funParams_;
pp.ChannelIndex      = iChanTal;
pp.SteerableOrder    = 4;      % even order -> ridge detector
pp.SigmaArray        = [1 2];  % px; tune to talin PSF
pp.BodyThreshold     = 'otsu'; % 'rosin' | 'otsu' | numeric
pp.GaussianBlurSigma = 2;      % px; blur before body threshold (reduces noise)
pp.BodyOpenRadius    = 8;      % px; opening removes filopodia roots (despike)
pp.BodyClosingRadius = 8;      % px; closing rounds/smooths the body edge
pp.ProcessFrames     = tFrame; % set [] for all frames
segProc.setPara(pp);

tic; segProc.run(); toc

%% --- P1 QC: body mask + ridge response ---
tFrameToInspect = 3;
img  = double(MD.channels_(iChanTal).loadImage(tFrameToInspect));
segDir = segProc.funParams_.OutputDirectory;
seg  = load(fullfile(segDir, sprintf('filoSeg_frame_%04d.mat', tFrameToInspect)));
bodyMask = seg.bodyMask;
res      = seg.res;

figure('Name','P1 QC','Color','w');
subplot(1,3,1); imshow(img,[100 500]); colormap(gca,gray); hold on;
visboundaries(bodyMask,'Color','c','LineWidth',0.5); hold off;
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

% merge any newly added default fields into existing params (e.g. UseRidgeTips)
defP2 = FilopodiaDetectionProcess.getDefaultParams(MD);
pp2 = detProc.funParams_;
fn = fieldnames(defP2);
for ii = 1:numel(fn)
    if ~isfield(pp2, fn{ii}), pp2.(fn{ii}) = defP2.(fn{ii}); end
end
pp2.ChannelIndex = iChanTal;
if ~isempty(MD.channels_(iChanTal).psfSigma_)
    pp2.PSFsigma = MD.channels_(iChanTal).psfSigma_;
else
    pp2.PSFsigma = 1.6;        % TUNE: ~(PSF FWHM in px)/2.355
end
pp2.Alpha              = 0.05; % TUNE: lower -> fewer detections
pp2.TipMaxDistFromBody = 70;   % TUNE: max filopodium reach (px)
pp2.MaxTipBaseDist     = 90;   % keep ~1.2 x TipMaxDistFromBody
pp2.BaseInsideBand     = 4;    % px inside body edge to include base adhesions
pp2.OrientLambda       = 3;    % TUNE: higher -> straighter shafts (tip mode only)
pp2.OrientTolerance    = 30;
pp2.TipRidgeBand       = -1;   % ridge gate (<0 = off)
pp2.ShaftAbsorbRadius  = 3;
pp2.CarveDistalFrac    = 0.8;
pp2.DetectMode         = 'auto'; % 'auto' -> 'all' for multi-frame, 'tip' for 1-frame
pp2.ProcessFrames      = tFrame;
detProc.setPara(pp2);

tic; detProc.run(); toc

%% --- P2 QC: overlay detections on image ---
detFile = fullfile(detProc.funParams_.OutputDirectory, 'filoDetection.mat');
S = load(detFile);
detectMode = S.detectMode;
fprintf('Frame %d: DetectMode = %s\n', tFrame, detectMode);

img = double(MD.channels_(iChanTal).loadImage(tFrameToInspect));
figure('Name','P2 QC','Color','w');

if strcmp(detectMode, 'all')
    % 'all' mode: show all detected adhesion puncta, colored by signed distance
    % to body edge (positive = outside, negative = inside).
    adh = S.adhesionInfo{tFrameToInspect};
    fprintf('Frame %d: %d adhesions detected (all mode)\n', tFrame, numel(adh));

    if ~isempty(adh)
        apos  = cat(1, adh.pos);          % [x y]
        adist = [adh.dist];               % signed dist: + outside, - inside
        imagesc(img); axis image off; colormap(gray); hold on;
        visboundaries(bodyMask,'Color','c','LineWidth',0.5);
        % color by distance: warm = distal (tip candidates), cool = proximal/inside
        sc = scatter(apos(:,1), apos(:,2), 20, adist, 'filled');
        colormap(gca, jet); cb = colorbar; ylabel(cb,'signed dist to edge (px)');
        caxis([-pp2.BaseInsideBand, pp2.TipMaxDistFromBody]);
        title(sprintf('all adhesions | n=%d | warm=distal, cool=inside body', numel(adh)));
    else
        imagesc(img); axis image off; colormap(gray);
        title('no adhesions detected');
    end

else
    % 'tip' mode: show filopodia (tip/base/shaft) -- same as before.
    fi = S.filoInfo{tFrame};
    fprintf('Frame %d: %d filopodia detected (tip mode)\n', tFrame, numel(fi));

    imagesc(img); axis image off; colormap(gray); hold on;
    for k = 1:numel(fi)
        cl = fi(k).centerline;
        plot(cl(:,1), cl(:,2), 'y-', 'LineWidth', 1.2);
        plot(fi(k).basePos(1), fi(k).basePos(2), 'co', 'MarkerSize',7,'LineWidth',1.5);
        plot(fi(k).tipPos(1),  fi(k).tipPos(2),  'r+', 'MarkerSize',9,'LineWidth',1.5);
    end
    title(sprintf('tip(+) base(o) shaft(-) | n=%d', numel(fi)));

    % Length distribution (px -> nm)
    if ~isempty(fi)
        L = [fi.length] * MD.pixelSize_;
        figure('Name','P2 lengths','Color','w');
        histogram(L); xlabel('filopodium length (nm)'); ylabel('count');
        title(sprintf('median = %.0f nm', median(L)));
    end
end

%% --- P2 QC: adhesion distance distribution ('all' mode) ---
if strcmp(detectMode, 'all') && ~isempty(adh)
    adist = [adh.dist];
    ares  = [adh.res];
    figure('Name','P2 adhesion distributions','Color','w','Position',[60 60 1200 420]);
    subplot(1,3,1); histogram(adist, 30); xlabel('signed dist to body edge (px)');
    ylabel('count'); title(sprintf('distance distribution (n=%d)', numel(adh)));
    xline(0,'r--','body edge','LabelVerticalAlignment','bottom');
    subplot(1,3,2); histogram(ares, 30); xlabel('steerable res at adhesion');
    title('ridge response at adhesions');
    subplot(1,3,3); scatter(adist, ares, 10, 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel('dist to edge (px)'); ylabel('res'); title('dist vs res (distal+bright = real tip?)');
    fprintf('dist: min %.1f  median %.1f  max %.1f px\n', min(adist), median(adist), max(adist));
    fprintf('res:  min %.1f  median %.1f  max %.1f\n', min(ares), median(ares), max(ares));
end

%% ===================== PROCESS 3: TRACKING =====================
% (only meaningful with full multi-frame detection; skip for single-frame QC)
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    if isempty(iP3)
        MD.addProcess(FilopodiaTrackingProcess(MD));
        iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    end
    trkProc = MD.processes_{iP3};
    pp3 = trkProc.funParams_;
    pp3.MaxLinkDist    = 10;   % px/frame; raise if filopodia move fast
    pp3.MaxGapFrames   = 2;    % frames; allow 2-frame gaps
    pp3.MinTrackLength = 3;    % discard tracks shorter than this
    % ensure output lands under FilopodiaForcePackage (may be wrong if
    % the process was created before the folder convention was set)
    pp3.OutputDirectory = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaTracking');
    trkProc.setPara(pp3);
    tic; trkProc.run(); toc

    trkFile = fullfile(trkProc.funParams_.OutputDirectory, 'filoTracks.mat');
    T = load(trkFile, 'adhesionTracks');
    ft = T.adhesionTracks;
    lifeF = [ft.lifetime];
    fprintf('Tracking: %d tracks | median lifetime %d frames | max %d frames\n', ...
        numel(ft), median(lifeF), max(lifeF));
    figure('Name','P3 lifetime','Color','w');
    histogram(lifeF); xlabel('track lifetime (frames)'); ylabel('count');
    title(sprintf('adhesion track lifetimes | n=%d', numel(ft)));
end

%% ===================== PROCESS 4: CLASSIFICATION (tip/base/shaft) =====================
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP4 = MD.getProcessIndex('FilopodiaClassificationProcess',1,0);
    if isempty(iP4)
        MD.addProcess(FilopodiaClassificationProcess(MD));
        iP4 = MD.getProcessIndex('FilopodiaClassificationProcess',1,0);
    end
    clsProc = MD.processes_{iP4};
    pp4 = clsProc.funParams_;
    pp4.MaxShaftLen     = 160;   % px; max shaft trace length (tip reach)
    pp4.BaseReachTol    = 2;     % px; base must land within this of body
    pp4.MinReachFrac    = 0.5;   % shaft must reach body in >= this fraction of frames
    pp4.MinTipLifetime  = 5;     % frames; tip track must persist this long
    pp4.OutputDirectory = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaClassification');
    clsProc.setPara(pp4);
    tic; clsProc.run(); toc

    clsFile = fullfile(clsProc.funParams_.OutputDirectory, 'filoClassification.mat');
    C = load(clsFile, 'filopodia', 'roleByTrack');
    nTip = numel(C.filopodia);
    fprintf('Filopodia: %d tip tracks (of %d) reached body and persisted\n', ...
        nTip, numel(C.roleByTrack));

    if nTip > 0
        figure('Name','P4 filopodium length L(t)','Color','w'); hold on;
        nShow = min(nTip, 8);
        for i = 1:nShow
            f = C.filopodia(i);
            plot(f.frames*MD.timeInterval_, f.L_nm, '-','LineWidth',1.2);
        end
        xlabel('time (s)'); ylabel('filopodium length L (nm)');
        title(sprintf('L(t) for %d of %d filopodia', nShow, nTip));
    end
end

%% ===================== VIEW IN movieViewer =====================
% After P1/P2/P3 have run on the full movie, open the framework viewer.
% Each process contributes toggleable overlays:
%   P1 -> 'Body mask', 'Ridge response (res)', 'Shaft mask'
%   P2 -> 'Adhesions' (all detected puncta, markers)
%   P3 -> 'Adhesion tracks' (u-track trajectories with dragtails)
% Use the checkboxes in the viewer's Overlay panel to turn them on/off,
% and the time slider to scrub frames.
movieViewer(MD);

%% ===================== TUNING NOTES =====================
% --- P1 body segmentation ---
%   pp.GaussianBlurSigma = 3;   % more blur -> smoother threshold
%   pp.BodyOpenRadius    = 12;  % larger -> cut deeper filopodia roots
%   pp.BodyClosingRadius = 12;  % larger -> rounder body edge
%   segProc.setPara(pp); segProc.run();
%
% --- P2 adhesion detection ('all' mode) ---
%   pp2.TipMaxDistFromBody = 100; % extend reach
%   pp2.BaseInsideBand     = 6;   % capture deeper base FAs
%   pp2.Alpha              = 0.01; % stricter -> fewer adhesions
%   detProc.setPara(pp2); detProc.run();
%
% --- P3 tracking ---
%   pp3.MaxLinkDist    = 15;   % allow faster-moving tips
%   pp3.MaxGapFrames   = 3;    % bridge longer disappearances
%   pp3.MinTrackLength = 5;    % stricter false-track filtering
%   trkProc.setPara(pp3); trkProc.run();