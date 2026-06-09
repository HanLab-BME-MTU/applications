function detectMovieFilopodia(movieData)
%DETECTMOVIEFILOPODIA  Process 2 wrapper. Tip/base puncta + shaft tracing.
%
% pointSourceDetection on the talin channel -> bright puncta. Classify by the
% body mask (inside/edge band -> base/FA; outside within range -> tip). Pair
% tip<->base by LAP, then trace each shaft as a min-cost path on the steerable
% response. Saves a tracker-compatible movieInfo and the rich filoInfo.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaDetectionProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaDetectionProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);

iChan = p.ChannelIndex;
frames = p.ProcessFrames; if isempty(frames), frames = 1:movieData.nFrames_; end

% locate P1 output directory directly from its funParams (robust to
% outFilePaths_ being empty after process reconstruction from disk)
iSeg = p.SegProcessIndex;
if isempty(iSeg)
    iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess', 1, 0);
end
assert(~isempty(iSeg), 'FilopodiaSegmentationProcess must be run first.');
segProc  = movieData.processes_{iSeg};
segOutDir = segProc.funParams_.OutputDirectory;
assert(exist(segOutDir, 'dir') == 7, ...
    'P1 output directory not found: %s\nRun FilopodiaSegmentationProcess first.', segOutDir);

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

%% helper: load one frame of P1 output directly from the mat file
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

    % --- bright talin puncta ---
    pstruct = pointSourceDetection(img, sigma, psArgs{:});
    if isempty(pstruct) || isempty(pstruct.x), continue; end
    px = pstruct.x(:); py = pstruct.y(:); pa = pstruct.A(:);

    % --- classify by body mask / distance to body ---
    [H, W] = size(bodyMask);
    rc = sub2ind([H W], min(max(round(py),1),H), min(max(round(px),1),W));
    inBody  = bodyMask(rc);
    distOut = bwdist(bodyMask);
    dPt = distOut(rc);
    isBase = inBody | (dPt <= p.BaseSearchBand);
    isTip  = (~inBody) & (dPt > p.BaseSearchBand) & (dPt <= p.TipMaxDistFromBody);

    tipXY  = [px(isTip),  py(isTip)];   tipA  = pa(isTip);
    baseXY = [px(isBase), py(isBase)];  baseA = pa(isBase);
    if isempty(tipXY) || isempty(baseXY), continue; end

    % --- global tip<->base assignment ---
    baseForTip = pairFilopodiaTipBase(tipXY, baseXY, p.MaxTipBaseDist);

    % --- trace shaft for each matched pair ---
    fi = struct('tipPos',{}, 'tipAmp',{}, 'basePos',{}, 'baseAmp',{}, ...
        'centerline',{}, 'arc',{}, 'length',{}, 'shaftMeanInt',{}, 'frame',{});
    for it = 1:size(tipXY, 1)
        jb = baseForTip(it);
        if isnan(jb), continue; end
        [cl, arc, len, ok] = traceFilopodiaShaft(res, theta, ...
            tipXY(it,:), baseXY(jb,:), p);
        if ~ok || len < p.MinFiloLength, continue; end
        idxLine = sub2ind([H W], min(max(round(cl(:,2)),1),H), ...
                                 min(max(round(cl(:,1)),1),W));
        e = numel(fi) + 1;
        fi(e).tipPos      = tipXY(it,:); fi(e).tipAmp  = tipA(it);
        fi(e).basePos     = baseXY(jb,:); fi(e).baseAmp = baseA(jb);
        fi(e).centerline  = cl;  fi(e).arc    = arc;
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
disp('Filopodia detection done.');
end
