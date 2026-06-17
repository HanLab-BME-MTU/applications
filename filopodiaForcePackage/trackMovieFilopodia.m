function trackMovieFilopodia(movieData)
%TRACKMOVIEFILOPODIA  Process 3 wrapper. Track talin adhesions over time.
%
% All adhesions detected in Process 2 ('all' mode) are linked across frames
% with u-track's Kalman tracker (trackCloseGapsKalmanSparse). No tip/base/shaft
% decision is made here: every adhesion is tracked so that tip adhesions are
% never lost. Tracks shorter than MinTrackLength frames are discarded (this
% removes transient false detections). Tip/base/shaft roles and per-filopodium
% length and velocity are assigned afterwards in Process 4 from these tracks.
%
% (Legacy 'tip' detection mode is single-frame and is not tracked.)
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaTrackingProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaTrackingProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

iDet = p.DetProcessIndex;
if isempty(iDet), iDet = movieData.getProcessIndex('FilopodiaDetectionProcess', 1, 0); end
assert(~isempty(iDet), 'FilopodiaDetectionProcess must be run first.');
detProc = movieData.processes_{iDet};
detFile = detProc.outFilePaths_{1, iChan};
if isempty(detFile) || exist(detFile,'file')~=2
    detFile = fullfile(detProc.funParams_.OutputDirectory, 'filoDetection.mat');
end
assert(exist(detFile,'file')==2, 'Detection output not found. Run Process 2 first.');
S = load(detFile);
movieInfo = S.movieInfo;
detectMode = 'all'; if isfield(S,'detectMode'), detectMode = S.detectMode; end
assert(strcmp(detectMode,'all'), ...
    ['Detection was run in ''%s'' mode. Tracking expects ''all'' mode ' ...
     '(multi-frame). Re-run P2 with DetectMode=''all''.'], detectMode);
adhesionInfo = S.adhesionInfo;

%% I/O
inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = detFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoTracks.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

%% ---- u-track configuration (reuse TrackingProcess defaults) ----
tw = p.MaxGapFrames;
gapCloseParam.timeWindow  = tw;
gapCloseParam.mergeSplit  = 0;
gapCloseParam.minTrackLen = 1;
gapCloseParam.diagnostics = 0;

kalmanFunctions = TrackingProcess.getKalmanFunctions(1);
kf = fieldnames(kalmanFunctions); valid = {'reserveMem','initialize','calcGain','timeReverse'};
kalmanFunctions = rmfield(kalmanFunctions, kf(~ismember(kf, valid)));

costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(movieData, tw, 1);
costMatrices(1).parameters.maxSearchRadius = p.MaxLinkDist;
costMatrices(1).parameters.minSearchRadius = 2;
costMatrices(1).parameters.diagnostics     = [];
costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(movieData, tw, 1);
costMatrices(2).parameters.maxSearchRadius = p.MaxLinkDist;

probDim = 2;
verbose = ~p.BatchMode;

%% ---- run the tracker ----
fprintf('Linking adhesions with u-track (timeWindow=%d, maxSearch=%g)...\n', tw, p.MaxLinkDist);
tracksFinal = trackCloseGapsKalmanSparse(movieInfo, costMatrices, gapCloseParam, ...
    kalmanFunctions, probDim, 0, verbose);

%% ---- convert tracksFinal -> adhesionTracks ----
dt  = movieData.timeInterval_;  if isempty(dt),  dt  = 1;  end
pix = movieData.pixelSize_;     if isempty(pix), pix = 1;  end

adhesionTracks = struct('trackId',{}, 'frames',{}, 'featIdx',{}, 'pos',{}, ...
    'amp',{}, 'dist',{}, 'dist_nm',{}, 'lifetime',{}, 'lifetime_s',{});

kept = 0;
keptOrig = [];
for i = 1:numel(tracksFinal)
    se = tracksFinal(i).seqOfEvents;
    sf = se(1,1); ef = se(2,1);
    featRow = tracksFinal(i).tracksFeatIndxCG(1, :);
    absFrames = sf:ef;
    present = featRow > 0;
    if nnz(present) < p.MinTrackLength, continue; end

    fr = absFrames(present); fIdx = featRow(present);
    nT = numel(fr);
    pos = nan(nT,2); amp = nan(1,nT); dist = nan(1,nT);
    for j = 1:nT
        a = adhesionInfo{fr(j)};
        k = fIdx(j);
        if isempty(a) || k > numel(a), continue; end
        pos(j,:) = a(k).pos; amp(j) = a(k).amp; dist(j) = a(k).dist;
    end
    good = ~isnan(amp);
    if nnz(good) < p.MinTrackLength, continue; end
    fr=fr(good); fIdx=fIdx(good); pos=pos(good,:); amp=amp(good); dist=dist(good);

    kept = kept + 1;
    keptOrig(kept) = i; %#ok<AGROW>
    adhesionTracks(kept).trackId    = kept;
    adhesionTracks(kept).frames     = fr;
    adhesionTracks(kept).featIdx    = fIdx;        % index into adhesionInfo{frame}
    adhesionTracks(kept).pos        = pos;         % [x y] per frame
    adhesionTracks(kept).amp        = amp;
    adhesionTracks(kept).dist       = dist;        % signed dist to body edge (px)
    adhesionTracks(kept).dist_nm    = dist * pix;
    adhesionTracks(kept).lifetime   = numel(fr);
    adhesionTracks(kept).lifetime_s = (fr(end)-fr(1)) * dt;
end

% keep the matching subset of raw u-track output for movieViewer (TracksDisplay)
if isempty(keptOrig)
    tracksFinal = tracksFinal([]);
else
    tracksFinal = tracksFinal(keptOrig);
end

save(outFile, 'adhesionTracks', 'tracksFinal', 'detectMode', '-v7.3');
fprintf('Adhesion tracking done: %d tracks kept (of %d raw) with >= %d frames.\n', ...
    kept, numel(tracksFinal), p.MinTrackLength);
fprintf('Next: Process 4 assigns tip/base/shaft roles and assembles filopodia.\n');
end
