%% runFilopodiaForce_debug.m
% Interactive driver to run and debug the FilopodiaForcePackage pipeline.
% Run cell by cell (Ctrl+Enter). Edit the CONFIG block first.
% Sangyoon J. Han / 2026

%% ===================== CONFIG =====================
mdPath   = '/mnt/nas/Collaborations/Celine_San_Diego/20260411/Soft-KO/KO_HELA_TALIN_G_BEAD_R_soft_40_nm_P_03_Airyscan_Processing_Processed_2_P3/KO HELA TALIN-G BEAD-R-soft 40 nm_P-03-Airyscan Processing-Processed 2.mat';
iChanTal = 1;             % talin-GFP channel index
tFrame   = [];             % [] for all frames; set to e.g. 5 for single-frame QC

if isempty(mdPath)
    assert(exist('MD','var')==1, 'Load MD first or set mdPath.');
else
    MD = MovieData.load(mdPath);
    MD.sanityCheck;
end
fprintf('Movie: %d channels, %d frames, pixelSize=%.4g nm, dt=%.4g s\n', ...
    numel(MD.channels_), MD.nFrames_, MD.pixelSize_, MD.timeInterval_);

%% ===================== PACKAGE REGISTRATION =====================
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
pp.SteerableOrder    = 4;
pp.SigmaArray        = [1 2];
pp.BodyThreshold     = 'otsu';
pp.GaussianBlurSigma = 2;
pp.BodyOpenRadius    = 8;
pp.BodyClosingRadius = 8;
pp.ProcessFrames     = tFrame;
segProc.setPara(pp);

tic; segProc.run(); toc

%% --- P1 QC: body mask + ridge response ---
tQC = 5;   % frame to inspect (change as needed)
img  = double(MD.channels_(iChanTal).loadImage(tQC));
segDir = segProc.funParams_.OutputDirectory;
seg  = load(fullfile(segDir, sprintf('filoSeg_frame_%04d.mat', tQC)));

figure('Name','P1 QC','Color','w');
subplot(1,3,1); imshow(img,[prctile(img(:),1) prctile(img(:),99)]); colormap(gca,gray); hold on;
visboundaries(seg.bodyMask,'Color','c','LineWidth',0.5); hold off;
title('talin + body boundary');
subplot(1,3,2); imshow(seg.res,[]); colormap(gca,hot);
title('steerable response (res)');
subplot(1,3,3); imagesc(seg.res .* ~seg.bodyMask); axis image off; colormap(gca,hot);
title('res outside body (shaft signal)');

%% --- P1 QC: shaftMask + body outline ---
figure('Name','P1 QC: shaft vs body','Color','w');
imshow(img,[prctile(img(:),1) prctile(img(:),99)]); colormap(gca,gray); hold on;
visboundaries(seg.bodyMask,'Color','c','LineWidth',0.5);
visboundaries(seg.shaftMask & ~seg.bodyMask,'Color',[0.2 0.6 1],'LineWidth',0.3);
title('body (cyan) + shaftMask outside body (blue)');

%% ===================== PROCESS 2: DETECTION =====================
iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
if isempty(iP2)
    MD.addProcess(FilopodiaDetectionProcess(MD));
    iP2 = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
end
detProc = MD.processes_{iP2};

defP2 = FilopodiaDetectionProcess.getDefaultParams(MD);
pp2 = detProc.funParams_;
fn = fieldnames(defP2);
for ii = 1:numel(fn)
    if ~isfield(pp2, fn{ii}), pp2.(fn{ii}) = defP2.(fn{ii}); end
end
pp2.ChannelIndex       = iChanTal;
pp2.PSFsigma           = 2.6;  % TIP scale (tips are larger than the PSF)
pp2.Alpha              = 0.05;
pp2.TipMaxDistFromBody = 100;
pp2.MaxTipBaseDist     = 120;
pp2.BaseInsideBand     = 4;
pp2.DetectMode         = 'auto';
pp2.ProcessFrames      = tFrame;
detProc.setPara(pp2);

tic; detProc.run(); toc

%% --- P2 QC: adhesion overlay on image ---
detFile = fullfile(detProc.funParams_.OutputDirectory, 'filoDetection.mat');
S = load(detFile);
detectMode = S.detectMode;
fprintf('DetectMode = %s\n', detectMode);

img = double(MD.channels_(iChanTal).loadImage(tQC));
figure('Name','P2 QC: adhesions','Color','w');

if strcmp(detectMode, 'all')
    adh = S.adhesionInfo{tQC};
    fprintf('Frame %d: %d adhesions detected\n', tQC, numel(adh));
    if ~isempty(adh)
        apos = cat(1, adh.pos); adist = [adh.dist];
        imagesc(img); axis image off; colormap(gray); hold on;
        visboundaries(seg.bodyMask,'Color','c','LineWidth',0.5);
        scatter(apos(:,1), apos(:,2), 20, adist, 'filled');
        colormap(gca, jet); cb = colorbar; ylabel(cb,'signed dist to edge (px)');
        caxis([-pp2.BaseInsideBand, pp2.TipMaxDistFromBody]);
        title(sprintf('all adhesions | n=%d | warm=distal', numel(adh)));
    end
end

%% --- P2 QC: distance & res distributions ---
if strcmp(detectMode,'all') && ~isempty(adh)
    adist = [adh.dist]; ares = [adh.res];
    figure('Name','P2 distributions','Color','w','Position',[60 60 1200 420]);
    subplot(1,3,1); histogram(adist,30); xlabel('signed dist to edge (px)'); ylabel('count');
    title(sprintf('dist distribution (n=%d)', numel(adh)));
    xline(0,'r--','body edge');
    subplot(1,3,2); histogram(ares,30); xlabel('steerable res'); title('res at adhesions');
    subplot(1,3,3); scatter(adist,ares,10,'filled','MarkerFaceAlpha',0.4);
    xlabel('dist (px)'); ylabel('res'); title('dist vs res');
end

%% ===================== PROCESS 3: TRACKING =====================
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    if isempty(iP3)
        MD.addProcess(FilopodiaTrackingProcess(MD));
        iP3 = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    end
    trkProc = MD.processes_{iP3};
    pp3 = trkProc.funParams_;
    pp3.MaxLinkDist    = 10;
    pp3.MaxGapFrames   = 2;
    pp3.MinTrackLength = 3;
    pp3.OutputDirectory = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaTracking');
    trkProc.setPara(pp3);
    tic; trkProc.run(); toc
end

%% --- P3 QC: track overlay on image ---
trkFile = fullfile(trkProc.funParams_.OutputDirectory, 'filoTracks.mat');
T = load(trkFile, 'adhesionTracks');
ft = T.adhesionTracks;
lifeF = [ft.lifetime];
fprintf('Tracking: %d tracks | median lifetime %d frames | max %d\n', ...
    numel(ft), median(lifeF), max(lifeF));

% overlay tracks at tQC: show each track as its trajectory up to tQC
img = double(MD.channels_(iChanTal).loadImage(tQC));
figure('Name',sprintf('P3 QC: tracks at frame %d',tQC),'Color','w');
imagesc(img); axis image off; colormap(gray); hold on;
visboundaries(seg.bodyMask,'Color','c','LineWidth',0.3);
colors = lines(min(numel(ft),256));
nShown = 0;
for i = 1:numel(ft)
    tr = ft(i);
    mask = tr.frames <= tQC & tr.frames >= max(1,tQC-20);  % show last 20 frames of trail
    if ~any(mask), continue; end
    pos = tr.pos(mask,:);
    ci = mod(i-1,size(colors,1))+1;
    plot(pos(:,1), pos(:,2), '-', 'Color', [colors(ci,:) 0.6], 'LineWidth', 0.8);
    if tr.frames(end) >= tQC && tr.frames(1) <= tQC
        % adhesion is present at tQC: plot current position
        idx = find(tr.frames == tQC, 1);
        if ~isempty(idx)
            plot(tr.pos(idx,1), tr.pos(idx,2), 'o', 'Color', colors(ci,:), 'MarkerSize', 4, 'LineWidth', 1);
            nShown = nShown + 1;
        end
    end
end
title(sprintf('tracks at frame %d: %d visible (trail=last 20 frames)', tQC, nShown));

% lifetime histogram
figure('Name','P3 QC: lifetimes','Color','w');
histogram(lifeF,30); xlabel('track lifetime (frames)'); ylabel('count');
title(sprintf('adhesion track lifetimes | n=%d | median=%d', numel(ft), median(lifeF)));

%% ===================== PROCESS 4: CLASSIFICATION =====================
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP4 = MD.getProcessIndex('FilopodiaClassificationProcess',1,0);
    if isempty(iP4)
        MD.addProcess(FilopodiaClassificationProcess(MD));
        iP4 = MD.getProcessIndex('FilopodiaClassificationProcess',1,0);
    end
    clsProc = MD.processes_{iP4};
    pp4 = clsProc.funParams_;
    pp4.MinTipLifetime  = 5;
    pp4.MinLinearFrac   = 0.85;
    pp4.MinTipDist      = 6;
    pp4.ShaftBand       = 4;
    pp4.BodyMaxAngle    = 60;
    pp4.MinReachFrac    = 0.5;
    pp4.LenPenalty      = 0.6;
    pp4.OutputDirectory = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaClassification');
    clsProc.setPara(pp4);
    tic; clsProc.run(); toc
end

%% --- P4 QC: tip/shaft/base overlay ---
clsFile = fullfile(clsProc.funParams_.OutputDirectory, 'filoClassification.mat');
C = load(clsFile);
nTip = numel(C.filopodia);
fprintf('Filopodia: %d well-tracked tips | weak-tip %d | background %d\n', ...
    nTip, sum(strcmp(C.roleByTrack,'weak-tip')), sum(strcmp(C.roleByTrack,'background')));

img = double(MD.channels_(iChanTal).loadImage(tQC));
figure('Name',sprintf('P4 QC: classification frame %d',tQC),'Color','w');
imshow(img,[100 500]); axis image off; colormap(gray); hold on;
visboundaries(seg.bodyMask,'Color','w','LineWidth',0.3);
% shafts (yellow lines)
sh = C.shaftByFrame{tQC};
if ~isempty(sh)
    plot(sh(:,1), sh(:,2), '-', 'Color', [1 1 0 0.7], 'LineWidth', 1.2);
end
% shaft adhesions (yellow dots)
sa = C.posByFrame.shaft{tQC};
if ~isempty(sa)
    plot(sa(:,1), sa(:,2), '.', 'Color', [1 1 0], 'MarkerSize', 10);
end
% base (cyan)
ba = C.posByFrame.base{tQC};
if ~isempty(ba)
    plot(ba(:,1), ba(:,2), 'co', 'MarkerSize', 6, 'LineWidth', 1.5);
end
% tips (red) last so they're on top
tp = C.posByFrame.tip{tQC};
if ~isempty(tp)
    plot(tp(:,1), tp(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
end
title(sprintf('frame %d: %d tips (red) | %d shaft adh (yellow) | %d base (cyan)', ...
    tQC, size(tp,1), size(sa,1), size(ba,1)));

%% --- P4 QC: L(t) plot ---
if nTip > 0
    figure('Name','P4 QC: filopodium length L(t)','Color','w'); hold on;
    nShow = min(nTip, 10);
    for i = 1:nShow
        f = C.filopodia(i);
        plot(f.frames*MD.timeInterval_, f.L_nm, '-', 'LineWidth', 1.2);
    end
    xlabel('time (s)'); ylabel('filopodium length L (nm)');
    title(sprintf('L(t) for %d of %d filopodia', nShow, nTip));
end

%% --- P4 diagnostic: tip eligibility ---
trkFile = fullfile(trkProc.funParams_.OutputDirectory, 'filoTracks.mat');
T = load(trkFile, 'adhesionTracks'); at = T.adhesionTracks;
nTr = numel(at); life = [at.lifetime]; maxDist = arrayfun(@(t) max(t.dist), at);
linFrac = zeros(1,nTr);
for i = 1:nTr
    P = at(i).pos; P = P-mean(P,1); Cv = (P'*P)/max(1,size(P,1)-1);
    ev = sort(eig(Cv),'descend'); linFrac(i) = ev(1)/max(sum(ev),eps);
end
pp4 = clsProc.funParams_;
passLife=life>=pp4.MinTipLifetime; passDist=maxDist>=pp4.MinTipDist;
passLin=linFrac>=pp4.MinLinearFrac | arrayfun(@(i) sqrt(sum(eig((at(i).pos-mean(at(i).pos))'*(at(i).pos-mean(at(i).pos))/max(1,size(at(i).pos,1)-1))))<3, 1:nTr);
fprintf('\n=== tip eligibility ===\n');
fprintf('pass lifetime>=%d: %d/%d\n', pp4.MinTipLifetime, sum(passLife), nTr);
fprintf('pass maxDist>=%.0f: %d/%d\n', pp4.MinTipDist, sum(passDist), nTr);
fprintf('pass linearity>=%.2f: %d/%d\n', pp4.MinLinearFrac, sum(passLin), nTr);
fprintf('pass ALL: %d/%d\n', sum(passLife&passDist&passLin), nTr);
fprintf('actual tips: %d\n', sum(strcmp(C.roleByTrack,'tip')));

%% ===================== VIEW IN movieViewer =====================
movieViewer(MD);

%% ===================== TUNING NOTES =====================
% --- P1 body segmentation ---
%   pp.GaussianBlurSigma = 3;   % more blur -> smoother threshold
%   pp.BodyOpenRadius    = 12;  % larger -> cut deeper filopodia roots
%   pp.BodyClosingRadius = 12;  % larger -> rounder body edge
%
% --- P2 adhesion detection ('all' mode) ---
%   pp2.PSFsigma           = 2.6; % TIP scale (tips > PSF)
%   pp2.TipMaxDistFromBody = 100; % extend reach
%   pp2.Alpha              = 0.01; % stricter -> fewer adhesions
%
% --- P3 tracking ---
%   pp3.MaxLinkDist    = 15;   % allow faster-moving tips
%   pp3.MaxGapFrames   = 3;    % bridge longer disappearances
%   pp3.MinTrackLength = 5;    % stricter false-track filtering
%
% --- P4 classification ---
%   pp4.MinTipLifetime = 8;    % longer-lived tracks only
%   pp4.MinLinearFrac  = 0.9;  % stricter linearity
%   pp4.BodyMaxAngle   = 45;   % tighter radial constraint
%   pp4.ShaftBand      = 5;    % wider shaft adhesion capture