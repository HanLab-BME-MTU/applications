function classifyMovieFilopodia(movieData)
%CLASSIFYMOVIEFILOPODIA  Process 4. Assemble filopodia from WELL-TRACKED tips.
%
% Logic (per the tip-centric design):
%   1. P2 detects all talin adhesions; P3 tracks them.
%   2. A track is tip-eligible only if it is WELL TRACKED: long lifetime,
%      linear trajectory, and distal (far from body). These are the adhesions
%      most likely to be true filopodium tips. Only these drive geometry.
%   3. For each tip (per frame), one STRAIGHT shaft direction to the body is
%      found (steerable orientation + collinear consensus); base = where that
%      line meets the body.
%   4. Other, less-well-tracked adhesions lying ON that tip->base line
%      (perpendicular distance < ShaftBand) are assigned as SHAFT adhesions of
%      that filopodium. Shaft adhesions do NOT get their own direction/base.
%      The most distal adhesion on a line stays the tip; nearer ones become
%      its shaft adhesions.
% Outputs per-frame tip/base/shaft-adhesion positions and shaft lines, plus a
% per-filopodium struct with length L(t) and velocity.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaClassificationProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaClassificationProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

iSeg = p.SegProcessIndex; if isempty(iSeg), iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess',1,0); end
iDet = p.DetProcessIndex; if isempty(iDet), iDet = movieData.getProcessIndex('FilopodiaDetectionProcess',1,0); end
iTrk = p.TrkProcessIndex; if isempty(iTrk), iTrk = movieData.getProcessIndex('FilopodiaTrackingProcess',1,0); end
assert(~isempty(iSeg)&&~isempty(iDet)&&~isempty(iTrk), 'Processes 1-3 must be run first.');
segProc = movieData.processes_{iSeg};
detProc = movieData.processes_{iDet};
trkProc = movieData.processes_{iTrk};
segOutDir = segProc.funParams_.OutputDirectory;

detFile = detProc.outFilePaths_{1,iChan};
if isempty(detFile)||exist(detFile,'file')~=2, detFile = fullfile(detProc.funParams_.OutputDirectory,'filoDetection.mat'); end
Sdet = load(detFile, 'adhesionInfo');
adhesionInfo = Sdet.adhesionInfo;

trkFile = trkProc.outFilePaths_{1,iChan};
if isempty(trkFile)||exist(trkFile,'file')~=2, trkFile = fullfile(trkProc.funParams_.OutputDirectory,'filoTracks.mat'); end
Strk = load(trkFile, 'adhesionTracks', 'tracksFinal');
adhesionTracks = Strk.adhesionTracks;
tracksFinalAll = Strk.tracksFinal;

%% I/O
inFilePaths = cell(1, numel(movieData.channels_)); inFilePaths{1,iChan} = trkFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoClassification.mat');
outFilePaths = cell(1, numel(movieData.channels_)); outFilePaths{1,iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

%% preload P1 maps
nF = movieData.nFrames_;
bodyByFrame = cell(1,nF); thetaByFrame = cell(1,nF); distByFrame = cell(1,nF);
for t = 1:nF
    fn = fullfile(segOutDir, sprintf('filoSeg_frame_%04d.mat', t));
    if exist(fn,'file')~=2, continue; end
    s = load(fn, 'bodyMask','theta');
    bodyByFrame{t}=s.bodyMask; thetaByFrame{t}=s.theta; distByFrame{t}=bwdist(s.bodyMask);
end

%% (A) per-track features -> tip eligibility (WELL TRACKED)
nTr = numel(adhesionTracks);
minLife   = getfielddef(p,'MinTipLifetime',5);
linFrac   = getfielddef(p,'MinLinearFrac',0.85);   % fraction of trajectory variance on 1 axis
minSpread = getfielddef(p,'MinTrajSpread',3);       % px; below this a track is "stationary" (linearity n/a)
minTipDist= getfielddef(p,'MinTipDist',6);          % px; tip must reach at least this far from body

tipEligible = false(1,nTr);
for i = 1:nTr
    tr = adhesionTracks(i);
    if numel(tr.frames) < minLife, continue; end
    if max(tr.dist) < minTipDist, continue; end
    % linearity of the trajectory (PCA): elongated-along-a-line OR ~stationary
    P = tr.pos; P = P - mean(P,1);
    C = (P'*P)/max(1,size(P,1)-1);
    ev = sort(eig(C),'descend');
    spread = sqrt(sum(ev));
    frac = ev(1)/max(sum(ev),eps);
    if spread < minSpread || frac >= linFrac
        tipEligible(i) = true;
    end
end

% map detected adhesions per frame -> owning track id and eligibility
trkOfAdh = cell(1,nF);      % trkOfAdh{t}(k) = track index owning adhesionInfo{t}(k), or 0
for t=1:nF
    if isempty(adhesionInfo{t}), trkOfAdh{t}=zeros(1,0); else, trkOfAdh{t}=zeros(1,numel(adhesionInfo{t})); end
end
for i = 1:nTr
    tr = adhesionTracks(i);
    for j = 1:numel(tr.frames)
        t = tr.frames(j); k = tr.featIdx(j);
        if t>=1 && t<=nF && k>=1 && k<=numel(trkOfAdh{t}), trkOfAdh{t}(k) = i; end
    end
end

%% (B) per-frame assembly: tips drive geometry; nearby adhesions -> shaft
shaftBand = getfielddef(p,'ShaftBand',4);    % px; perpendicular distance to tip->base line
posByFrame.tip   = cell(1,nF);
posByFrame.base  = cell(1,nF);
posByFrame.shaft = cell(1,nF);    % shaft adhesions
shaftByFrame     = cell(1,nF);    % shaft lines (polylines, NaN separated)
roleByAdh        = cell(1,nF);    % per detected adhesion role at frame t
% per-track accumulation for L(t)
tipFrames = cell(1,nTr); tipLen = cell(1,nTr); tipPos = cell(1,nTr); tipBase = cell(1,nTr);

progressText(0,'Filopodia assembly','Filopodia assembly');
for t = 1:nF
    posByFrame.tip{t}=zeros(0,2); posByFrame.base{t}=zeros(0,2); posByFrame.shaft{t}=zeros(0,2);
    shaftByFrame{t}=zeros(0,2);
    a = adhesionInfo{t};
    if isempty(a) || isempty(bodyByFrame{t}), progressText(t/nF,'Filopodia assembly'); continue; end
    A = cat(1, a.pos); distA = [a.dist];
    n = size(A,1);
    roleByAdh{t} = repmat({'none'}, 1, n);

    elig = find(arrayfun(@(k) trkOfAdh{t}(k)>0 && tipEligible(trkOfAdh{t}(k)), 1:n));
    [~, ord] = sort(distA(elig), 'descend');   % most distal first
    elig = elig(ord);
    claimed = false(1,n);

    for kk = elig
        if claimed(kk), continue; end          % already a shaft of a more distal tip
        [shaft, base, ~, ~, reached] = straightShaftToBody(A(kk,1), A(kk,2), A, ...
            thetaByFrame{t}, bodyByFrame{t}, distByFrame{t}, p);
        if ~reached, continue; end
        % assign other adhesions on the tip->base segment as shaft adhesions
        seg = base - A(kk,:); Lseg = hypot(seg(1),seg(2));
        if Lseg < 1, continue; end
        u = seg / Lseg;                         % unit tip->base
        rel = A - A(kk,:);
        proj = rel(:,1)*u(1) + rel(:,2)*u(2);
        perp = abs(-rel(:,1)*u(2) + rel(:,2)*u(1));
        onLine = (perp < shaftBand) & (proj > 1) & (proj < Lseg) & (distA(:) < distA(kk)) & ~claimed(:);
        shIdx = find(onLine);
        claimed(shIdx) = true; claimed(kk) = true;

        roleByAdh{t}{kk} = 'tip';
        for q = shIdx.', roleByAdh{t}{q} = 'shaft'; end
        posByFrame.tip{t}(end+1,:)  = A(kk,:);
        posByFrame.base{t}(end+1,:) = base;
        if ~isempty(shIdx), posByFrame.shaft{t} = [posByFrame.shaft{t}; A(shIdx,:)]; end
        shaftByFrame{t} = [shaftByFrame{t}; shaft; nan(1,2)];

        % accumulate for this tip's track
        it = trkOfAdh{t}(kk);
        tipFrames{it}(end+1) = t;
        tipLen{it}(end+1)    = Lseg;
        tipPos{it}(end+1,:)  = A(kk,:);
        tipBase{it}(end+1,:) = base;
    end
    progressText(t/nF,'Filopodia assembly');
end

%% (C) per-filopodium struct (accepted tip tracks) + roleByTrack
dt  = movieData.timeInterval_;  if isempty(dt),  dt  = 1; end
pix = movieData.pixelSize_;     if isempty(pix), pix = 1; end
sw  = max(1, round(getfielddef(p,'VelSmoothWin',3)));
reachFrac = getfielddef(p,'MinReachFrac',0.5);

roleByTrack = repmat({'background'},1,nTr);
filopodia = struct('tipTrackId',{},'frames',{},'tipPos',{},'basePos',{}, ...
    'L',{},'L_nm',{},'velocity_nmps',{},'lifetime',{});
for i = 1:nTr
    if ~tipEligible(i) || isempty(tipFrames{i}), continue; end
    % accepted if it acted as a tip (reached body) in enough of its frames
    if numel(tipFrames{i}) < reachFrac * numel(adhesionTracks(i).frames), 
        roleByTrack{i} = 'weak-tip'; continue;
    end
    roleByTrack{i} = 'tip';
    L = tipLen{i}; Lsm = smoothLocal(L, sw);
    n = numel(filopodia)+1;
    filopodia(n).tipTrackId    = adhesionTracks(i).trackId;
    filopodia(n).frames        = tipFrames{i};
    filopodia(n).tipPos        = tipPos{i};
    filopodia(n).basePos       = tipBase{i};
    filopodia(n).L             = L;
    filopodia(n).L_nm          = L * pix;
    filopodia(n).velocity_nmps = gradient(Lsm)/dt * pix;
    filopodia(n).lifetime      = numel(tipFrames{i});
end

% tip tracks (u-track format) for TracksDisplay
tipTrkIdx = find(strcmp(roleByTrack,'tip'));
if ~isempty(tipTrkIdx) && ~isempty(tracksFinalAll)
    safe = tipTrkIdx(tipTrkIdx <= numel(tracksFinalAll));
    tipTracks = tracksFinalAll(safe);
else
    tipTracks = tracksFinalAll([]);
end

save(outFile, 'filopodia','tipTracks','roleByTrack','posByFrame','shaftByFrame','roleByAdh','-v7.3');
progressText(1,'Filopodia assembly');
fprintf('Assembly: %d filopodia (well-tracked tips) of %d tracks; %d tip-eligible.\n', ...
    numel(filopodia), nTr, nnz(tipEligible));
end

% ===================================================================
function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

% ===================================================================
function y = smoothLocal(x, w)
if w <= 1 || numel(x) < 3, y = x; return; end
k = ones(1, w) / w; y = conv(x, k, 'same');
h = floor(w/2);
for i = 1:min(h, numel(x))
    y(i) = mean(x(1:min(numel(x), i+h)));
    y(end-i+1) = mean(x(max(1,end-i+1-h):end));
end
end
