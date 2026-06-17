function classifyMovieFilopodia(movieData)
%CLASSIFYMOVIEFILOPODIA  Process 4 wrapper. Assemble filopodia by tracing
%shafts from tracked tip adhesions to the body.
%
% Approach (tip-centric, intensity-driven):
%   1. P2 detects bright talin tip adhesions (blobs) outside the body.
%   2. P3 tracks them.
%   3. Here, for each tracked tip at each frame, the shaft is traced from the
%      tip toward the body along the local steerable ridge orientation
%      (straightShaftToBody). The line endpoint on the body = base.
%   4. A track is accepted as a filopodium only if its shaft reaches the body
%      in a sufficient fraction of frames and it persists long enough; this
%      rejects background blobs (no coherent shaft to the body) and transient
%      noise. Filopodium length L(t) = traced shaft length; velocity = dL/dt.
% No ridge-mask topology / skeleton endpoints are used, so per-frame results
% are stable across frames (blobs track well).
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaClassificationProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaClassificationProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

iSeg = p.SegProcessIndex; if isempty(iSeg), iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess',1,0); end
iTrk = p.TrkProcessIndex; if isempty(iTrk), iTrk = movieData.getProcessIndex('FilopodiaTrackingProcess',1,0); end
assert(~isempty(iSeg) && ~isempty(iTrk), 'Processes 1-3 must be run before classification.');
segProc = movieData.processes_{iSeg};
trkProc = movieData.processes_{iTrk};
segOutDir = segProc.funParams_.OutputDirectory;

trkFile = trkProc.outFilePaths_{1, iChan};
if isempty(trkFile) || exist(trkFile,'file')~=2
    trkFile = fullfile(trkProc.funParams_.OutputDirectory, 'filoTracks.mat');
end
Strk = load(trkFile, 'adhesionTracks', 'tracksFinal');
adhesionTracks = Strk.adhesionTracks;
tracksFinalAll = Strk.tracksFinal;

%% I/O
inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = trkFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoClassification.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

%% preload per-frame body + theta from P1
nF = movieData.nFrames_;
bodyByFrame  = cell(1, nF);
thetaByFrame = cell(1, nF);
distByFrame  = cell(1, nF);
for t = 1:nF
    fn = fullfile(segOutDir, sprintf('filoSeg_frame_%04d.mat', t));
    if exist(fn,'file')~=2, continue; end
    s = load(fn, 'bodyMask', 'theta');
    bodyByFrame{t}  = s.bodyMask;
    thetaByFrame{t} = s.theta;
    distByFrame{t}  = bwdist(s.bodyMask);     % distance outside body
end

% all tracked tip positions per frame (collinear support for shaft direction)
nTr = numel(adhesionTracks);
allPtsByFrame = cell(1, nF);
for t = 1:nF, allPtsByFrame{t} = zeros(0,2); end
for i = 1:nTr
    tr = adhesionTracks(i);
    for j = 1:numel(tr.frames)
        t = tr.frames(j);
        if t>=1 && t<=nF
            allPtsByFrame{t}(end+1,:) = tr.pos(j,:);
        end
    end
end

%% straight shaft for every tracked tip, every frame
maxLen  = round(getfielddef(p,'MaxShaftLen',160));
reachTol= getfielddef(p,'BaseReachTol',2);    % px; base counts as on body if within this

posByFrame.tip  = cell(1,nF);
posByFrame.base = cell(1,nF);
shaftByFrame    = cell(1,nF);
for t = 1:nF, posByFrame.tip{t}=zeros(0,2); posByFrame.base{t}=zeros(0,2); shaftByFrame{t}=zeros(0,2); end

trackShaft = cell(1,nTr);     % per track: struct array over frames (path, base, len, reached)

progressText(0, 'Filopodia assembly', 'Filopodia assembly');
for i = 1:nTr
    tr = adhesionTracks(i);
    sh = struct('frame',{},'path',{},'base',{},'len',{},'reached',{});
    for j = 1:numel(tr.frames)
        t = tr.frames(j);
        if t<1 || t>nF || isempty(bodyByFrame{t}), continue; end
        xy = tr.pos(j,:);
        [pth, base, ~, ~, reached] = straightShaftToBody(xy(1), xy(2), allPtsByFrame{t}, ...
            thetaByFrame{t}, bodyByFrame{t}, distByFrame{t}, p);
        L = sum(sqrt(sum(diff(pth,1,1).^2,2)));   % straight length (px)
        k = numel(sh)+1;
        sh(k).frame=t; sh(k).path=pth; sh(k).base=base; sh(k).len=L; sh(k).reached=reached;
    end
    trackShaft{i} = sh;
    progressText(i/nTr, 'Filopodia assembly');
end

%% accept filopodium tracks: shaft reaches body often enough + persistent
minLife   = getfielddef(p,'MinTipLifetime',5);
reachFrac = getfielddef(p,'MinReachFrac',0.5);
dt  = movieData.timeInterval_;  if isempty(dt),  dt  = 1; end
pix = movieData.pixelSize_;     if isempty(pix), pix = 1; end
sw  = max(1, round(getfielddef(p,'VelSmoothWin',3)));

roleByTrack = repmat({'background'}, 1, nTr);
filopodia = struct('tipTrackId',{},'frames',{},'tipPos',{},'basePos',{}, ...
    'L',{},'L_nm',{},'velocity_nmps',{},'lifetime',{});
isTip = false(1,nTr);
for i = 1:nTr
    sh = trackShaft{i};
    if isempty(sh), continue; end
    reached = [sh.reached];
    if numel(sh) >= minLife && mean(reached) >= reachFrac
        isTip(i) = true;
        roleByTrack{i} = 'tip';
        fr = [sh.frame]; L = [sh.len];
        Lsm = smoothLocal(L, sw);
        n = numel(filopodia)+1;
        filopodia(n).tipTrackId    = adhesionTracks(i).trackId;
        filopodia(n).frames        = fr;
        filopodia(n).tipPos        = cell2mat(arrayfun(@(s) s.path(1,:), sh, 'unif',0)');
        filopodia(n).basePos       = cat(1, sh.base);
        filopodia(n).L             = L;
        filopodia(n).L_nm          = L * pix;
        filopodia(n).velocity_nmps = gradient(Lsm)/dt * pix;
        filopodia(n).lifetime      = numel(fr);
        % accumulate display positions per frame
        for j = 1:numel(sh)
            t = sh(j).frame;
            posByFrame.tip{t}(end+1,:)  = sh(j).path(1,:);
            posByFrame.base{t}(end+1,:) = sh(j).base;
            shaftByFrame{t} = [shaftByFrame{t}; sh(j).path; nan(1,2)]; % NaN separates filopodia
        end
    end
end

% u-track-format subset for TracksDisplay (accepted tip tracks only)
tipIdx = find(isTip);
if ~isempty(tipIdx) && ~isempty(tracksFinalAll)
    safeIdx = tipIdx(tipIdx <= numel(tracksFinalAll));
    tipTracks = tracksFinalAll(safeIdx);
else
    tipTracks = tracksFinalAll([]);
end

save(outFile, 'filopodia', 'tipTracks', 'roleByTrack', 'posByFrame', 'shaftByFrame', '-v7.3');
progressText(1, 'Filopodia assembly');
fprintf('Assembly done: %d filopodia (tip tracks) of %d adhesion tracks.\n', numel(tipIdx), nTr);
end

% ===================================================================
function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

% ===================================================================
function y = smoothLocal(x, w)
if w <= 1 || numel(x) < 3, y = x; return; end
k = ones(1, w) / w;
y = conv(x, k, 'same');
h = floor(w/2);
for i = 1:min(h, numel(x))
    y(i) = mean(x(1:min(numel(x), i+h)));
    y(end-i+1) = mean(x(max(1,end-i+1-h):end));
end
end