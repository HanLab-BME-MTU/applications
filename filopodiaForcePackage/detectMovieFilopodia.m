function detectMovieFilopodia(movieData)
%DETECTMOVIEFILOPODIA  Process 2 wrapper. Detect talin adhesions.
%
% Two modes (DetectMode):
%   'all' (default for multi-frame movies): detect EVERY talin adhesion
%       punctum near the cell edge as a point feature, without deciding which
%       is a tip. All adhesions (tip + shaft + base) go into a u-track-ready
%       movieInfo so that nothing is pre-selected and tip adhesions cannot be
%       underestimated. Tip/base/shaft roles are assigned AFTER tracking
%       (Process 4) from the tracked adhesions.
%   'tip' (default for single-frame movies): for a standalone frame with no
%       tracking, pick one tip adhesion per filopodium and trace its shaft to
%       the body (legacy behavior), so filopodia can be compared directly.
%
% Frames run in parallel (parfor) with a progress bar.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaDetectionProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaDetectionProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

% resolve detect mode
mode = lower(p.DetectMode);
if strcmp(mode,'auto'), if movieData.nFrames_ > 1, mode = 'all'; else, mode = 'tip'; end; end

% locate P1 output
iSeg = p.SegProcessIndex;
if isempty(iSeg), iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess', 1, 0); end
assert(~isempty(iSeg), 'FilopodiaSegmentationProcess must be run first.');
segProc   = movieData.processes_{iSeg};
segOutDir = segProc.funParams_.OutputDirectory;
assert(exist(segOutDir, 'dir') == 7, 'P1 output dir not found: %s', segOutDir);

segFiles  = dir(fullfile(segOutDir, 'filoSeg_frame_*.mat'));
segFrames = cellfun(@(n) sscanf(n, 'filoSeg_frame_%d.mat'), {segFiles.name});
frames = p.ProcessFrames; if isempty(frames), frames = 1:movieData.nFrames_; end
frames = intersect(frames, segFrames);
assert(~isempty(frames), 'No frames have P1 output. Run P1 first.');

sigma = p.PSFsigma;
psArgs = {'Alpha', p.Alpha};
if ~isempty(p.ConfRadius), psArgs = [psArgs, {'ConfRadius', p.ConfRadius}]; end
if ~isempty(p.WindowSize), psArgs = [psArgs, {'WindowSize', p.WindowSize}]; end

%% I/O
inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = segOutDir;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoDetection.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

%% preload images (avoid BioFormatsReader contention in parfor)
nF = movieData.nFrames_;  nProc = numel(frames);
chan = movieData.channels_(iChan);
fprintf('Detection mode: %s.  Loading %d frames...\n', mode, nProc);
imgs = cell(1, nProc);
for k = 1:nProc, imgs{k} = double(chan.loadImage(frames(k))); end

progressText(0, 'Filopodia detection', 'Filopodia detection');
useDQ = ~isempty(ver('parallel'));
if useDQ, dq = parallel.pool.DataQueue; afterEach(dq, @(~) progressTick(nProc)); end

adhCell  = cell(1, nProc);     % all-adhesion mode output
filoCell = cell(1, nProc);     % tip mode output
miCell   = cell(1, nProc);

parfor k = 1:nProc
    t = frames(k);
    img = imgs{k};
    s = load(fullfile(segOutDir, sprintf('filoSeg_frame_%04d.mat', t)));
    if strcmp(mode,'all')
        [adh, mi] = detectAdhesionsAll(img, s.bodyMask, s.res, sigma, psArgs, p, t);
        adhCell{k} = adh;  miCell{k} = mi;
    else
        [fi, mi] = detectTips(img, s.bodyMask, s.res, s.theta, sigma, psArgs, p, t);
        filoCell{k} = fi;  miCell{k} = mi;
    end
    if useDQ, send(dq, 1); end
end

%% assemble outputs indexed by absolute frame
movieInfo(1:nF) = struct('xCoord', [], 'yCoord', [], 'amp', []);
for k = 1:nProc
    if ~isempty(miCell{k}), movieInfo(frames(k)) = miCell{k}; end
end

if strcmp(mode,'all')
    adhesionInfo = cell(1, nF);
    for k = 1:nProc, adhesionInfo{frames(k)} = adhCell{k}; end
    detectMode = 'all';
    save(outFile, 'movieInfo', 'adhesionInfo', 'detectMode', '-v7.3');
else
    filoInfo = cell(1, nF);
    for k = 1:nProc, filoInfo{frames(k)} = filoCell{k}; end
    detectMode = 'tip';
    save(outFile, 'movieInfo', 'filoInfo', 'detectMode', '-v7.3');
end

progressText(1, 'Filopodia detection');
fprintf('Filopodia detection done (mode=%s). Frames: %s\n', mode, mat2str(frames));
end

% ===================================================================
function progressTick(N)
persistent c
if isempty(c), c = 0; end
c = c + 1;
progressText(c / N, 'Filopodia detection');
if c >= N, c = []; end
end

% ===================================================================
function [adh, mi] = detectAdhesionsAll(img, bodyMask, res, sigma, psArgs, p, t)
% Detect ALL talin adhesion puncta near the cell edge as point features.
adh = struct('pos',{}, 'amp',{}, 'dist',{}, 'res',{}, 'frame',{});
mi  = struct('xCoord', [], 'yCoord', [], 'amp', []);
[H, W] = size(bodyMask);

pstruct = pointSourceDetection(img, sigma, psArgs{:});
if isempty(pstruct) || isempty(pstruct.x), return; end
px = pstruct.x(:); py = pstruct.y(:); pa = pstruct.A(:);

% signed distance to body edge: + outside, - inside
sd = bwdist(bodyMask) - bwdist(~bodyMask);
rc = sub2ind([H W], min(max(round(py),1),H), min(max(round(px),1),W));
dsignal = sd(rc);

% keep adhesions from just inside the edge (base FAs) out to the max reach
keep = (dsignal >= -p.BaseInsideBand) & (dsignal <= p.TipMaxDistFromBody);
if ~any(keep), return; end
px = px(keep); py = py(keep); pa = pa(keep); dk = dsignal(keep);
rk = res(sub2ind([H W], min(max(round(py),1),H), min(max(round(px),1),W)));

n = numel(px);
for i = 1:n
    adh(i).pos = [px(i), py(i)]; adh(i).amp = pa(i);
    adh(i).dist = dk(i); adh(i).res = rk(i); adh(i).frame = t;
end
mi.xCoord = [px, zeros(n,1)];
mi.yCoord = [py, zeros(n,1)];
mi.amp    = [pa, zeros(n,1)];
end

% ===================================================================
function [fi, mi] = detectTips(img, bodyMask, res, theta, sigma, psArgs, p, t)
% Single-frame mode: one tip adhesion per filopodium + shaft trace to body.
fi = struct('tipPos',{}, 'tipAmp',{}, 'basePos',{}, 'baseAmp',{}, ...
    'centerline',{}, 'arc',{}, 'length',{}, 'shaftMeanInt',{}, 'frame',{});
mi = struct('xCoord', [], 'yCoord', [], 'amp', []);
[H, W] = size(bodyMask);

pstruct = pointSourceDetection(img, sigma, psArgs{:});
if isempty(pstruct) || isempty(pstruct.x), return; end
px = pstruct.x(:); py = pstruct.y(:); pa = pstruct.A(:);

rc = sub2ind([H W], min(max(round(py),1),H), min(max(round(px),1),W));
inBody  = bodyMask(rc);
distOut = bwdist(bodyMask);
dPt = distOut(rc);
isTip = (~inBody) & (dPt > p.BaseSearchBand) & (dPt <= p.TipMaxDistFromBody);
if ~any(isTip), return; end
tipXY = [px(isTip), py(isTip)];  tipA = pa(isTip);  tipD = dPt(isTip);

[~, ord] = sort(tipD, 'descend');
tipXY = tipXY(ord,:); tipA = tipA(ord); tipD = tipD(ord);
blockedMask = false(H, W);
absR = max(1, round(p.ShaftAbsorbRadius));

for it = 1:size(tipXY, 1)
    tcol = min(max(round(tipXY(it,1)),1),W);
    trow = min(max(round(tipXY(it,2)),1),H);
    if blockedMask(trow, tcol), continue; end
    Rhint = tipD(it) * 1.8 + 15;
    [cl, arc, len, basePos, ok] = traceFilopodiaToBody(res, theta, ...
        bodyMask, tipXY(it,:), p, blockedMask, Rhint);
    if ~ok || len < p.MinFiloLength, continue; end
    idxLine = sub2ind([H W], min(max(round(cl(:,2)),1),H), min(max(round(cl(:,1)),1),W));
    nPts = size(cl,1);
    keepFrom = max(1, round((1 - p.CarveDistalFrac) * nPts));
    carve = false(H, W); carve(idxLine(keepFrom:end)) = true;
    if absR > 0, carve = imdilate(carve, strel('disk', absR)); end
    blockedMask = blockedMask | carve;
    bpx = min(max(round(basePos(1)),1),W); bpy = min(max(round(basePos(2)),1),H);
    e = numel(fi) + 1;
    fi(e).tipPos = tipXY(it,:); fi(e).tipAmp = tipA(it);
    fi(e).basePos = basePos;    fi(e).baseAmp = img(bpy, bpx);
    fi(e).centerline = cl; fi(e).arc = arc;
    fi(e).length = len;   fi(e).shaftMeanInt = mean(img(idxLine));
    fi(e).frame = t;
end
if ~isempty(fi)
    tips = cat(1, fi.tipPos);
    mi.xCoord = [tips(:,1), zeros(size(tips,1),1)];
    mi.yCoord = [tips(:,2), zeros(size(tips,1),1)];
    mi.amp    = [cat(1,fi.tipAmp), zeros(numel(fi),1)];
end
end
