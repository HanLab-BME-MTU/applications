function detectMovieFilopodia(movieData)
%DETECTMOVIEFILOPODIA  Process 2 wrapper. Tip puncta + tip->body shaft trace.
%
% pointSourceDetection on the talin channel finds bright puncta; puncta
% outside the body within range are tip candidates. Each tip is traced to the
% body edge with traceFilopodiaToBody, which defines the base geometrically as
% the body-boundary contact point. No punctum-based base detection or LAP
% pairing is needed. Tips whose trace cannot reach the body are dropped.
% Saves a tracker-compatible movieInfo and the rich filoInfo.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaDetectionProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaDetectionProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);

iChan = p.ChannelIndex;

% locate P1 output directory directly from its funParams
iSeg = p.SegProcessIndex;
if isempty(iSeg)
    iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess', 1, 0);
end
assert(~isempty(iSeg), 'FilopodiaSegmentationProcess must be run first.');
segProc   = movieData.processes_{iSeg};
segOutDir = segProc.funParams_.OutputDirectory;
assert(exist(segOutDir, 'dir') == 7, ...
    'P1 output directory not found: %s\nRun FilopodiaSegmentationProcess first.', segOutDir);

% frames = intersection of requested frames and frames P1 produced
segFiles  = dir(fullfile(segOutDir, 'filoSeg_frame_*.mat'));
segFrames = cellfun(@(n) sscanf(n, 'filoSeg_frame_%d.mat'), {segFiles.name});
frames = p.ProcessFrames; if isempty(frames), frames = 1:movieData.nFrames_; end
missing = setdiff(frames, segFrames);
if ~isempty(missing)
    warning(['P1 output missing for %d frame(s) (e.g. frame %d). Skipping. ' ...
        'Re-run P1 with matching ProcessFrames.'], numel(missing), missing(1));
end
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

    function s = loadSegFrame(t)
        fname = fullfile(segOutDir, sprintf('filoSeg_frame_%04d.mat', t));
        assert(exist(fname,'file')==2, 'P1 frame file missing: %s', fname);
        s = load(fname);
    end

%% run
nF = movieData.nFrames_;
filoInfo = cell(1, nF);
movieInfo(1:nF) = struct('xCoord', [], 'yCoord', [], 'amp', []);

for t = frames
    img = double(movieData.channels_(iChan).loadImage(t));
    seg = loadSegFrame(t);
    bodyMask = seg.bodyMask;
    res      = seg.res;
    theta    = seg.theta;
    [H, W] = size(bodyMask);

    % --- bright talin puncta ---
    pstruct = pointSourceDetection(img, sigma, psArgs{:});
    if isempty(pstruct) || isempty(pstruct.x), continue; end
    px = pstruct.x(:); py = pstruct.y(:); pa = pstruct.A(:);

    % --- tip candidates: puncta outside the body, within reach ---
    rc = sub2ind([H W], min(max(round(py),1),H), min(max(round(px),1),W));
    inBody  = bodyMask(rc);
    distOut = bwdist(bodyMask);                       % 0 inside, dist outside
    dPt = distOut(rc);
    isTip = (~inBody) & (dPt > p.BaseSearchBand) & (dPt <= p.TipMaxDistFromBody);
    tipXY = [px(isTip), py(isTip)];  tipA = pa(isTip);
    if isempty(tipXY), continue; end

    % --- trace each tip to the body; base = body-edge contact ---
    fi = struct('tipPos',{}, 'tipAmp',{}, 'basePos',{}, 'baseAmp',{}, ...
        'centerline',{}, 'arc',{}, 'length',{}, 'shaftMeanInt',{}, 'frame',{});
    for it = 1:size(tipXY, 1)
        [cl, arc, len, basePos, ok] = traceFilopodiaToBody(res, theta, ...
            bodyMask, tipXY(it,:), p);
        if ~ok || len < p.MinFiloLength, continue; end
        idxLine = sub2ind([H W], min(max(round(cl(:,2)),1),H), ...
                                 min(max(round(cl(:,1)),1),W));
        bpx = min(max(round(basePos(1)),1),W);
        bpy = min(max(round(basePos(2)),1),H);
        e = numel(fi) + 1;
        fi(e).tipPos      = tipXY(it,:); fi(e).tipAmp  = tipA(it);
        fi(e).basePos     = basePos;     fi(e).baseAmp = img(bpy, bpx);
        fi(e).centerline  = cl;  fi(e).arc = arc;
        fi(e).length      = len; fi(e).shaftMeanInt = mean(img(idxLine));
        fi(e).frame       = t;
    end
    filoInfo{t} = fi;

    % tracker-compatible detections (tips)
    if ~isempty(fi)
        tips = cat(1, fi.tipPos);
        movieInfo(t).xCoord = [tips(:,1), zeros(size(tips,1),1)];
        movieInfo(t).yCoord = [tips(:,2), zeros(size(tips,1),1)];
        movieInfo(t).amp    = [cat(1,fi.tipAmp), zeros(numel(fi),1)];
    end
end

save(outFile, 'movieInfo', 'filoInfo', '-v7.3');
fprintf('Filopodia detection done. Frames: %s\n', mat2str(frames));
end
