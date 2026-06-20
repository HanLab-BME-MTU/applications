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
pp2.PSFsigma           = 2.1;  % TUNE: lower = more sensitive; 2.1 balances tip size vs noise
pp2.Alpha              = 0.05;
pp2.TipMaxDistFromBody = 110;
pp2.MaxTipBaseDist     = 130;
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
    pp3.MaxLinkDist    = 8;    % px/frame; small: generous detection keeps detections close frame-to-frame
    pp3.MaxGapFrames   = 3;    % frames; bridge short disappearances
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
imagesc(img); axis image off; colormap(gray); hold on;
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

%% ===================== PROCESS 5: FORCE / INTENSITY SAMPLING =====================
% Requires the TFMPackage ForceFieldCalculationProcess to have been run on MD
% (cross-package: traction is read from forceField.mat).
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP5 = MD.getProcessIndex('FilopodiaSamplingProcess',1,0);
    if isempty(iP5)
        MD.addProcess(FilopodiaSamplingProcess(MD));
        iP5 = MD.getProcessIndex('FilopodiaSamplingProcess',1,0);
    end
    smpProc = MD.processes_{iP5};
    pp5 = smpProc.funParams_;
    pp5.ChannelIndex     = iChanTal;
    pp5.ShaftSampleStep  = 3;     % px; arc-length step along shaft
    pp5.SampleRadius     = 1;     % px; local averaging radius for talin
    pp5.OutputDirectory  = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaSampling');
    smpProc.setPara(pp5);
    tic; smpProc.run(); toc
end

%% --- P5 QC: tip traction & talin distributions, and one shaft profile ---
smpFile = fullfile(smpProc.funParams_.OutputDirectory, 'filoSamples.mat');
Q = load(smpFile, 'filoSamples');
fs = Q.filoSamples;
fprintf('Sampling: %d filopodia\n', numel(fs));

% pool tip traction & talin across all filopodia/frames
tipF = []; tipFa = []; tipI = [];
for f = 1:numel(fs)
    tipF  = [tipF,  fs(f).tipForce];
    tipFa = [tipFa, fs(f).tipForceAxial];
    tipI  = [tipI,  fs(f).tipTalin];
end
goodF = isfinite(tipF);
fprintf('tip traction: %d/%d sampled | median %.1f Pa\n', nnz(goodF), numel(tipF), median(tipF(goodF)));

figure('Name','P5 QC: tip force & talin','Color','w','Position',[60 60 1200 380]);
subplot(1,3,1); histogram(tipF(goodF),30); xlabel('tip traction magnitude (Pa)'); ylabel('count');
title(sprintf('tip traction | median %.1f Pa', median(tipF(goodF))));
subplot(1,3,2); histogram(tipFa(isfinite(tipFa)),30); xlabel('tip axial traction (Pa, +=outward)');
title('axial traction at tip'); xline(0,'r--');
subplot(1,3,3); scatter(tipI(goodF), tipF(goodF), 12, 'filled', 'MarkerFaceAlpha',0.4);
xlabel('tip talin intensity'); ylabel('tip traction (Pa)'); title('talin vs traction at tip');

% one example shaft profile (force + talin vs arc length)
if ~isempty(fs)
    f = 1; j = 1;   % first filopodium, first frame
    sp = fs(f).shaftProfile{j};
    if ~isempty(sp) && ~isempty(sp.s_nm)
        figure('Name','P5 QC: shaft profile','Color','w');
        yyaxis left;  plot(sp.s_nm, sp.force, '-o','LineWidth',1.2); ylabel('traction (Pa)');
        yyaxis right; plot(sp.s_nm, sp.talin, '-s','LineWidth',1.2); ylabel('talin intensity');
        xlabel('arc length tip->base (nm)');
        title(sprintf('shaft profile: filopodium %d, frame %d', f, fs(f).frames(j)));
    end
end

%% ===================== PROCESS 6: STATISTICS =====================
if strcmp(detectMode,'all') && isempty(pp2.ProcessFrames)
    iP6 = MD.getProcessIndex('FilopodiaStatisticsProcess',1,0);
    if isempty(iP6)
        MD.addProcess(FilopodiaStatisticsProcess(MD));
        iP6 = MD.getProcessIndex('FilopodiaStatisticsProcess',1,0);
    end
    statProc = MD.processes_{iP6};
    pp6 = statProc.funParams_;
    pp6.ChannelIndex        = iChanTal;
    pp6.MinLifetimeForStats = 3;
    pp6.OutputDirectory     = fullfile(MD.outputDirectory_, ...
        'FilopodiaForcePackage', 'FilopodiaStatistics');
    statProc.setPara(pp6);
    tic; statProc.run(); toc
end

%% --- P6 QC: pooled distributions + correlations ---
statFile = fullfile(statProc.funParams_.OutputDirectory, 'filoStats.mat');
R = load(statFile, 'stats');
st = R.stats;
fprintf('Stats: %d filopodia\n', st.nFilopodia);
fprintf('corr talin~axial r=%.2f | vProt~axial r=%.2f | L~tipForce r=%.2f | talin~force r=%.2f\n', ...
    st.corr.talin_vs_axial, st.corr.vProt_vs_axial, st.corr.L_vs_tipForce, st.corr.talin_vs_force);

figure('Name','P6 QC: filopodium statistics','Color','w','Position',[40 40 1300 700]);
subplot(2,3,1); histogram(st.pooled.Lmean_nm,20); xlabel('mean length (nm)'); ylabel('count'); title('length');
subplot(2,3,2); histogram(st.pooled.lifetime_s,20); xlabel('lifetime (s)'); title('lifetime');
subplot(2,3,3); histogram(st.pooled.vProt_nmps(isfinite(st.pooled.vProt_nmps)),20);
xlabel('protrusion speed (nm/s)'); title('elongation speed');
subplot(2,3,4); histogram(st.pooled.tipForceMean(isfinite(st.pooled.tipForceMean)),20);
xlabel('tip traction (Pa)'); title('tip traction');
subplot(2,3,5);
x = st.pooled.tipTalinMean; y = st.pooled.tipAxialMean; ok = isfinite(x)&isfinite(y);
scatter(x(ok),y(ok),14,'filled','MarkerFaceAlpha',0.4);
xlabel('tip talin'); ylabel('tip axial traction (Pa)');
title(sprintf('talin vs axial (r=%.2f)', st.corr.talin_vs_axial));
subplot(2,3,6);
x = st.pooled.vProt_nmps; y = st.pooled.tipAxialMean; ok = isfinite(x)&isfinite(y);
scatter(x(ok),y(ok),14,'filled','MarkerFaceAlpha',0.4);
xlabel('protrusion speed (nm/s)'); ylabel('tip axial traction (Pa)');
title(sprintf('speed vs axial (r=%.2f)', st.corr.vProt_vs_axial));

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

%% --- P1 QC: steerable theta colormap + adhesion overlay ---
% Run this cell after P1 and P2 have completed to diagnose orientation issues.
tQC_theta = 5;   % change as needed
img_qc = double(MD.channels_(iChanTal).loadImage(tQC_theta));
seg_qc = load(fullfile(segProc.funParams_.OutputDirectory, ...
    sprintf('filoSeg_frame_%04d.mat', tQC_theta)), 'bodyMask','theta','res');
det_qc = S.adhesionInfo{tQC_theta};

% theta from steerableDetector is in [-pi/2, pi/2] (ridge orientation).
% Normalize to [0, pi) for display: mod(theta, pi).
% BUT: theta=-45 and theta=135 are the same ridge (mod pi ambiguity).
% Use mod(theta, pi) which correctly maps both to the same hue.
th  = seg_qc.theta;          % [-pi/2, pi/2] ridge orientation
r   = seg_qc.res;            % steerable response strength
bm  = seg_qc.bodyMask;

% mask: show only outside body where res is meaningful
outside = ~bm;
r_norm  = r ./ (max(r(outside & r>0)) + eps);   % [0,1]

% HSV colormap: hue encodes angle.
% theta in [-pi/2,pi/2]: map to [0,1] via (theta + pi/2) / pi
% so -90 deg (=90 deg ridge) -> hue=0 (red), 0 deg (horiz) -> hue=0.5 (cyan),
% +90 deg (vert) -> hue=1 (red again). Use mod for clean [0,1].
H_ch = mod(th + pi/2, pi) / pi;   % hue: -90->0(red), 0->0.5(cyan), 90->1(red)
S_ch = ones(size(th));
V_ch = r_norm .* outside;   % dark inside body
hsv_img = cat(3, H_ch, S_ch, V_ch);
rgb_theta = hsv2rgb(hsv_img);

figure('Name',sprintf('P1 QC: steerable theta (frame %d)',tQC_theta), ...
    'Color','w','Position',[50 50 1300 500]);

% panel 1: theta colormap (hue=direction, brightness=strength)
ax1 = subplot(1,3,1); imshow(rgb_theta); hold on;
visboundaries(bm,'Color','w','LineWidth',0.4);
if ~isempty(det_qc)
    apos = cat(1,det_qc.pos);
    plot(apos(:,1),apos(:,2),'ko','MarkerSize',7,'LineWidth',1.2);
end
title('theta (hue=direction, bright=strong ridge)');
% colorbar showing angle mapping
colormap(ax1, hsv(256));
cb = colorbar(ax1); cb.Ticks=[0 0.25 0.5 0.75 1];
cb.TickLabels={'-90°','-45°','0°','45°','90°'};
ylabel(cb,'ridge orientation theta (horiz=cyan, vert=red/both ends)');

% panel 2: res map outside body
ax2 = subplot(1,3,2);
imagesc(r .* outside); axis image off; colormap(ax2, hot); colorbar(ax2);
hold on; visboundaries(bm,'Color','c','LineWidth',0.4);
title('steerable response (res) outside body');

% panel 3: theta quiver: show local orientation as short lines on image
ax3 = subplot(1,3,3);
imshow(img_qc,[prctile(img_qc(:),1) prctile(img_qc(:),99)]); colormap(ax3,gray); hold on;
visboundaries(bm,'Color','c','LineWidth',0.4);
% subsample grid for quiver
step = 8;
[yy,xx] = meshgrid(1:step:size(th,1), 1:step:size(th,2));
yy=yy'; xx=xx';
mask_q = outside(sub2ind(size(th),min(yy,size(th,1)),min(xx,size(th,2))));
r_q    = r(sub2ind(size(r),min(yy(:),size(r,1)),min(xx(:),size(r,2))));
th_q   = th(sub2ind(size(th),min(yy(:),size(th,1)),min(xx(:),size(th,2))));
valid  = mask_q(:) & r_q(:) > prctile(r_q(r_q>0),50);   % only strong ridges
sc = step*0.45;
for qi = find(valid)'
    a = th_q(qi);
    x0_ = xx(qi); y0_ = yy(qi);
    dx = sc*cos(a); dy = sc*sin(a);
    plot([x0_-dx, x0_+dx],[y0_-dy, y0_+dy],'-','Color',[0.2 0.9 0.2],'LineWidth',0.8);
end
% overlay adhesions with their theta as a colored dot
if ~isempty(det_qc)
    apos=cat(1,det_qc.pos);
    ath=arrayfun(@(a) th(min(max(round(a.pos(2)),1),size(th,1)), ...
        min(max(round(a.pos(1)),1),size(th,2))), det_qc);
    scatter(apos(:,1),apos(:,2),50, mod(ath,pi)/pi,'filled','MarkerEdgeColor','k','LineWidth',0.5);
    colormap(ax3,hsv(256));
end
title(sprintf('theta quiver (green) + adhesion theta (colored dots) | frame %d', tQC_theta));

%% ===================== DETECTION STABILITY QC =====================
% Run after P2 to diagnose missing tips per frame.
detFile2 = fullfile(detProc.funParams_.OutputDirectory,'filoDetection.mat');
Sd = load(detFile2,'adhesionInfo');
nPerFrame = cellfun(@numel, Sd.adhesionInfo);
fprintf('\n=== P2 Detection stability ===\n');
fprintf('PSFsigma used: %.2f px\n', pp2.PSFsigma);
fprintf('Alpha: %.3f\n', pp2.Alpha);
fprintf('detections/frame: median=%d  min=%d  max=%d  std=%.1f\n', ...
    median(nPerFrame), min(nPerFrame), max(nPerFrame), std(double(nPerFrame)));
fprintf('frames with 0 detections: %d\n', sum(nPerFrame==0));
fprintf('frames with < half median: %d\n', sum(nPerFrame < median(nPerFrame)/2));

figure('Name','P2 QC: detections per frame','Color','w');
plot(nPerFrame,'b-o','MarkerSize',3); xlabel('frame'); ylabel('# detections');
title(sprintf('detections/frame | median=%d min=%d max=%d', ...
    median(nPerFrame),min(nPerFrame),max(nPerFrame)));
yline(median(nPerFrame),'r--','median');

%% ===================== TRACKING PARAMETER + STABILITY QC =====================
trkFile2 = fullfile(trkProc.funParams_.OutputDirectory,'filoTracks.mat');
T2 = load(trkFile2,'adhesionTracks');
at = T2.adhesionTracks;
nTr = numel(at);

fprintf('\n=== P3 Tracking parameters ===\n');
fprintf('MaxLinkDist:    %.1f px\n', pp3.MaxLinkDist);
fprintf('MaxGapFrames:   %d frames\n', pp3.MaxGapFrames);
fprintf('MinTrackLength: %d frames\n', pp3.MinTrackLength);

lifeF = [at.lifetime];
fprintf('\n=== P3 Track statistics (n=%d) ===\n', nTr);
fprintf('lifetime: median=%d  p25=%d  p75=%d  max=%d frames\n', ...
    median(lifeF), prctile(lifeF,25), prctile(lifeF,75), max(lifeF));
fprintf('  short (<=3 frames): %d (%.0f%%)\n', sum(lifeF<=3), 100*mean(lifeF<=3));
fprintf('  long  (>=10 frames): %d (%.0f%%)\n', sum(lifeF>=10), 100*mean(lifeF>=10));
fprintf('  very long (>=20): %d (%.0f%%)\n', sum(lifeF>=20), 100*mean(lifeF>=20));

% per-frame: how many tracks are ACTIVE (present) at each frame?
nActive = zeros(1,MD.nFrames_);
nNew    = zeros(1,MD.nFrames_);   % tracks starting this frame
nEnd    = zeros(1,MD.nFrames_);   % tracks ending this frame
for i = 1:nTr
    fr = at(i).frames;
    nActive(fr) = nActive(fr)+1;
    nNew(fr(1))   = nNew(fr(1))+1;
    nEnd(fr(end)) = nEnd(fr(end))+1;
end
fprintf('active tracks/frame: median=%d  min=%d  max=%d\n', ...
    median(nActive),min(nActive),max(nActive));
fprintf('new tracks/frame:    median=%d  max=%d\n', median(nNew),max(nNew));
fprintf('ending tracks/frame: median=%d  max=%d\n', median(nEnd),max(nEnd));

figure('Name','P3 QC: tracking stability','Color','w','Position',[50 50 1300 500]);
subplot(1,3,1);
histogram(lifeF,30,'FaceColor','b');
xlabel('lifetime (frames)'); ylabel('count'); title('track lifetime distribution');
xline(pp3.MinTrackLength,'r--',sprintf('min=%d',pp3.MinTrackLength));
xline(5,'g--','tip threshold=5');

subplot(1,3,2);
plot(nActive,'b-','LineWidth',1); hold on;
plot(nNew,'r-','LineWidth',0.8);
plot(nEnd,'g-','LineWidth',0.8);
xlabel('frame'); ylabel('count');
legend('active','new','ending','Location','best');
title('active / new / ending tracks per frame');

subplot(1,3,3);
% gap analysis: how many tracks have gaps (disappeared then reappeared)?
hasGap = arrayfun(@(t) numel(t.frames)<t.lifetime, at);
fprintf('tracks with gaps: %d/%d (%.0f%%)\n', sum(hasGap),nTr,100*mean(hasGap));
gapLen = arrayfun(@(t) ...
    sum(diff(sort(t.frames))>1), at);  % # of gap events
histogram(gapLen(gapLen>0),0:10,'FaceColor','r');
xlabel('# gap events per track'); ylabel('count');
title(sprintf('gap analysis | %d/%d tracks have gaps (MaxGapFrames=%d)', ...
    sum(hasGap),nTr,pp3.MaxGapFrames));

%% --- trajectory overlay: tracks present at tQC, showing last 30-frame trail ---
img_trk = double(MD.channels_(iChanTal).loadImage(tQC));
figure('Name',sprintf('P3 QC: trajectory trails at frame %d',tQC),'Color','w');
imagesc(img_trk); axis image off; colormap(gray); hold on;
visboundaries(seg.bodyMask,'Color','c','LineWidth',0.3);
nShown=0; cmap=jet(nTr);
for i=1:nTr
    tr=at(i);
    if tr.frames(end)<tQC || tr.frames(1)>tQC, continue; end
    mask=tr.frames>=max(1,tQC-29) & tr.frames<=tQC;
    pos=tr.pos(mask,:);
    if size(pos,1)<2, continue; end
    ci=mod(i-1,size(cmap,1))+1;
    plot(pos(:,1),pos(:,2),'-','Color',[cmap(ci,:) 0.6],'LineWidth',0.8);
    % current position (at tQC)
    j=find(tr.frames==tQC,1);
    if ~isempty(j)
        plot(tr.pos(j,1),tr.pos(j,2),'o','Color',cmap(ci,:),'MarkerSize',4,'LineWidth',1);
        nShown=nShown+1;
    end
end
title(sprintf('trajectory trails (last 30 frames) at frame %d | %d active tracks',tQC,nShown));
fprintf('\nTUNING SUGGESTIONS:\n');
fprintf('  - Too many short tracks / new tracks each frame -> lower MaxLinkDist or raise MinTrackLength\n');
fprintf('  - Tracks disappear then reappear often -> raise MaxGapFrames\n');
fprintf('  - Missing tip detections -> lower Alpha (e.g. 0.01) or raise PSFsigma\n');
