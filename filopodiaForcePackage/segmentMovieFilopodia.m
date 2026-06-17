function segmentMovieFilopodia(movieData)
%SEGMENTMOVIEFILOPODIA  Process 1 wrapper. Cell body mask + steerable maps.
%
% Body = blurred threshold + opening (despike filopodia roots) + closing
% (round/smooth edge). Shaft enhancement = multiscale steerable ridge
% response with hysteresis. Saves per-frame:
%   bodyMask, res, theta, nms, scaleMap, shaftMask.
% Frames run in parallel (parfor) with a progress bar.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaSegmentationProcess', 1, 0);
assert(~isempty(iProc), ['No FilopodiaSegmentationProcess found. ' ...
    'Create it and call process.run.']);
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);

iChan = p.ChannelIndex;
frames = p.ProcessFrames; if isempty(frames), frames = 1:movieData.nFrames_; end
nProc = numel(frames);

%% I/O
inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = movieData.getChannelPaths{iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outDir;
proc.setOutFilePaths(outFilePaths);

%% preload images (avoid BioFormatsReader contention in parfor)
chan = movieData.channels_(iChan);
fprintf('Segmentation: loading %d frames...\n', nProc);
imgs = cell(1, nProc);
for k = 1:nProc, imgs{k} = double(chan.loadImage(frames(k))); end

%% progress bar
progressText(0, 'Filopodia segmentation', 'Filopodia segmentation');
dq = parallel.pool.DataQueue;
afterEach(dq, @(~) progressTick(nProc));

%% run (parfor over frames)
parfor k = 1:nProc
    t = frames(k);
    img = imgs{k};
    result = segmentOneFrame(img, p);
    fname = fullfile(outDir, sprintf('filoSeg_frame_%04d.mat', t));
    saveSegFrame(fname, result);
    send(dq, 1);
end

progressText(1, 'Filopodia segmentation');
fprintf('Filopodia segmentation done. Frames: %s\n', mat2str(frames));
end

% ===================================================================
function saveSegFrame(fname, r)
% Save one frame's maps. Isolating save() in a function (not the parfor body)
% avoids the parfor transparency restriction on save/load.
bodyMask = r.bodyMask; res = r.res; theta = r.theta;
nms = r.nms; scaleMap = r.scaleMap; shaftMask = r.shaftMask;
save(fname, 'bodyMask', 'res', 'theta', 'nms', 'scaleMap', 'shaftMask');
end

% ===================================================================
function progressTick(N)
persistent c
if isempty(c), c = 0; end
c = c + 1;
progressText(c / N, 'Filopodia segmentation');
if c >= N, c = []; end
end

% ===================================================================
function r = segmentOneFrame(img, p)
% Segment one frame: body mask + steerable ridge maps.

% --- cell body: blurred threshold, then despike + smooth ---
% Filopodia are recovered separately by the steerable filter, so the body
% must exclude them and end smooth/rounded. Opening removes the thin
% filopodia roots; closing rounds the boundary and fills the notches
% between them.
imgB = imgaussfilt(img, p.GaussianBlurSigma);
switch lower(num2str(p.BodyThreshold))
    case 'otsu',  level = thresholdOtsu(imgB);
    case 'rosin', level = thresholdRosin(imgB);
    otherwise,    level = p.BodyThreshold;
end
bodyMask = imgB > level;
bodyMask = imfill(bodyMask, 'holes');
if any(bodyMask(:))
    bodyMask = bwareafilt(bodyMask, 1);              % keep largest object
end
if p.BodyOpenRadius > 0
    bodyMask = imopen(bodyMask, strel('disk', p.BodyOpenRadius));   % despike
end
if p.BodyClosingRadius > 0
    bodyMask = imclose(bodyMask, strel('disk', p.BodyClosingRadius)); % round
end
bodyMask = imfill(bodyMask, 'holes');
bodyMask = bwareaopen(bodyMask, p.BodyMinArea);

% --- steerable ridge maps (on the original unblurred image) ---
[res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, ...
    p.SteerableOrder, p.SigmaArray);

% --- hysteresis shaft mask ---
bg = res(~bodyMask);
hi = p.HysteresisHigh; lo = p.HysteresisLow;
if isempty(hi), hi = median(bg) + 3 * mad(bg,1) * 1.4826; end
if isempty(lo), lo = median(bg) + 1 * mad(bg,1) * 1.4826; end
seed = (res > hi) & ~bodyMask;
cand = (res > lo) & ~bodyMask;
shaftMask = imreconstruct(seed & cand, cand);

r.bodyMask  = bodyMask;
r.res       = res;
r.theta     = theta;
r.nms       = nms;
r.scaleMap  = scaleMap;
r.shaftMask = shaftMask;
end
