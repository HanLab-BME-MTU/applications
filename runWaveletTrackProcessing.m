%[data] = runTrackProcessing_Wavelets(data, varargin)
%
% INPUTS    data        : array of experiment structures
%           {Buffer}    : length of buffer before/after tracks. Default: 5
%           {Overwrite} : true to overwrite previous results
%           {FileName}  : name of the output
%
% Outputs
%
%           nSeg : number of segments; 1: regular, 2 or more: merging/splitting
%     visibility : 1) complete track, 2) incomplete (cut off), 3) persistent
%      gapStatus : 4) valid gap, 5) invalid gap
%
% Usage example: runTrackProcessing(data, 'Buffer', 3);
%
% Notes: The buffer size influences the number of visible tracks
%
% Francois Aguet, November 2010 (last modified 05/18/2012)

function runWaveletTrackProcessing(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Buffer', [5 5], @(x) numel(x)==2);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackerOutput', 'trackedFeatures_v1.mat', @ischar);
ip.addParamValue('FileName', 'ProcessedTracks_v1.mat', @ischar);
ip.addParamValue('FrameIndexes', arrayfun(@(x) 1:x.movieLength, data, 'UniformOutput', false), @(x) numel(unique(diff(x)))==1); %check that frame rate is constant
ip.addParamValue('Sigma', [], @(x) numel(x)==length(data(1).channels));
ip.addParamValue('Preprocess', true, @islogical);
ip.addParamValue('Postprocess', true, @islogical);
ip.addParamValue('Cutoff_f', 2, @isscalar);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 125 150]); % used in post-proc
ip.parse(data, varargin{:});
filename = ip.Results.FileName;
overwrite = ip.Results.Overwrite;
buffer = ip.Results.Buffer;
trackerOutput = ip.Results.TrackerOutput;
frameIdx = ip.Results.FrameIndexes;
sigma = ip.Results.Sigma;
preprocess = ip.Results.Preprocess;
postprocess = ip.Results.Postprocess;
cohortBounds = ip.Results.CohortBounds_s;

for i = 1:length(data)
    data(i).tracks = [];
    data(i).smTracks = [];
end
parfor i = 1:length(data)
    if ~(exist([data(i).source filesep 'Tracking' filesep filename],'file')==2) || overwrite
        data(i) = main(data(i), buffer, trackerOutput, filename, frameIdx{i}, sigma, preprocess, postprocess, cohortBounds);
    else
        fprintf('Tracks from %s have already been processed.\n', getShortPath(data(i)));
    end
end



function [data] = main(data, buffer, trackerOutput, filename, frameIdx, sigmaV, preprocess, postprocess, cohortBounds)
cutoff_f = 2;
minLft = cutoff_f*data.framerate;
cohortBounds(cohortBounds<=minLft) = [];
cohortBounds = [minLft cohortBounds data.movieLength*data.framerate];

detection = load([data.source 'Detection' filesep 'detection_v1.mat']);
frameInfo = detection.frameInfo;

ny = data.imagesize(1);
nx = data.imagesize(2);
nFrames = length(frameIdx);

%=================================
% Identify master/slave channels
%=================================
nCh = length(data.channels);
mCh = strcmp(data.source, data.channels);

if isempty(sigmaV)
    sigmaV = zeros(nCh, 1);
    for k = 1:nCh
        sigmaV(k) = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k});
        data.framePaths{k} = data.framePaths{k}(frameIdx);
    end
end

if isempty(data.maskPaths)
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    data.maskPaths = arrayfun(@(x) [maskPath x.name], dir([maskPath '*.tif']), 'UniformOutput', false);
end
data.maskPaths = data.maskPaths(frameIdx);

data.framerate = data.framerate*(frameIdx(2)-frameIdx(1));
data.movieLength = length(data.framePaths{1});

sigma = sigmaV(mCh);
w4 = ceil(4*sigma);

%======================================================================
% Read and convert tracker output
%======================================================================
tPath = [data.source 'Tracking' filesep trackerOutput];
if exist(tPath, 'file')==2
    trackinfo = load(tPath);
    trackinfo = trackinfo.tracksFinal;
    nTracks = length(trackinfo);
else
    error('No valid tracker output found.');
end


%======================================================================
% Preprocessing
%======================================================================
if preprocess
% Remove single-frame tracks
bounds = arrayfun(@(i) i.seqOfEvents([1 end],1), trackinfo, 'UniformOutput', false);
rmIdx = diff(horzcat(bounds{:}), [], 1)==0;
trackinfo(rmIdx) = [];
nTracks = size(trackinfo, 1);

%----------------------------------------------------------------------
% Merge compound tracks with overlapping ends/starts
%----------------------------------------------------------------------
for i = 1:nTracks
    nSeg = size(trackinfo(i).tracksFeatIndxCG,1);
    if nSeg > 1
        seqOfEvents = trackinfo(i).seqOfEvents;
        tracksCoordAmpCG = trackinfo(i).tracksCoordAmpCG;
        tracksFeatIndxCG = trackinfo(i).tracksFeatIndxCG;
        
        rmIdx = [];
        for s = 1:nSeg
            iEvent = seqOfEvents(seqOfEvents(:,3)==s,:);
            parentSeg = iEvent(2,4);
            parentStartIdx = seqOfEvents(seqOfEvents(:,2)==1 & seqOfEvents(:,3)==parentSeg,1);
            
            % conditions for merging:
            % -current segment merges at end
            % -overlap between current segment and 'parent' it merges into: 1 frame
            if ~isnan(iEvent(2,4)) &&...
                    iEvent(2,1)-1==parentStartIdx
   
                % replace start index of parent with start index of current segment
                seqOfEvents(seqOfEvents(:,3)==parentSeg & seqOfEvents(:,2)==1,1) = iEvent(1,1);
                % remove current segment
                seqOfEvents(seqOfEvents(:,3)==s,:) = [];
                % assign segments that merge/split from current to parent
                seqOfEvents(seqOfEvents(:,4)==s,4) = parentSeg;
                
                % use distance of points at overlap to assign
                xMat = tracksCoordAmpCG(:,1:8:end);
                yMat = tracksCoordAmpCG(:,2:8:end);
                                
                % indexes in the 8-step matrices
                iMat = repmat(1:size(xMat,2), [nSeg 1]).*~isnan(xMat);
                
                overlapIdx = setdiff(intersect(iMat(parentSeg,:), iMat(s,:)), 0);
                if overlapIdx(1)>1 && overlapIdx(end)<seqOfEvents(end,1) && overlapIdx(1)~=(iEvent(1,1)-seqOfEvents(1,1)+1)
                    idx = [overlapIdx(1)-1 overlapIdx(end)+1];
                    if isnan(xMat(s,idx(1)))
                        idx(1) = overlapIdx(1);
                    end
                    if isnan(xMat(parentSeg,idx(2)))
                        idx(2) = overlapIdx(end);
                    end
                elseif overlapIdx(1)==1 || overlapIdx(1)==(iEvent(1,1)-seqOfEvents(1,1)+1)
                    idx = [overlapIdx(1) overlapIdx(end)+1];
                    if isnan(xMat(parentSeg,idx(2)))
                        idx(2) = overlapIdx(end);
                    end
                else
                    idx = [overlapIdx(1)-1 overlapIdx(end)];
                    if isnan(xMat(s,idx(1)))
                        idx(1) = overlapIdx(1);
                    end
                end
                xRef = interp1(idx, [xMat(s,idx(1)) xMat(parentSeg,idx(2))], overlapIdx);
                yRef = interp1(idx, [yMat(s,idx(1)) yMat(parentSeg,idx(2))], overlapIdx);
                
                d = sqrt((xMat([s parentSeg],overlapIdx)-xRef).^2 + (yMat([s parentSeg],overlapIdx)-yRef).^2);
                % remove overlap
                rm = [s parentSeg];
                rm = rm(d~=min(d));
                iMat(rm,overlapIdx) = 0;
                tracksCoordAmpCG(rm,(overlapIdx-1)*8+(1:8)) = NaN;
                tracksFeatIndxCG(rm,overlapIdx) = 0;
                tracksFeatIndxCG(parentSeg,iMat(s,:)~=0) = tracksFeatIndxCG(s,iMat(s,:)~=0);
                
                % concatenate segments
                range8 = iMat(s,:);
                range8(range8==0) = [];
                range8 = (range8(1)-1)*8+1:range8(end)*8;
                tracksCoordAmpCG(parentSeg, range8) = tracksCoordAmpCG(s, range8);
                
                rmIdx = [rmIdx s]; %#ok<AGROW>
            end % segment loop
        end
        rmIdx = unique(rmIdx);
        tracksFeatIndxCG(rmIdx,:) = [];
        tracksCoordAmpCG(rmIdx,:) = [];
        
        % re-order seqOfEvents
        [~,ridx] = sort(seqOfEvents(:,1));
        %[~,lidx] = sort(ridx);
        seqOfEvents = seqOfEvents(ridx,:);
        
        % indexes in seqOfEvents must be in order of segment appearance
        % replace with unique(seqOfEvents(:,3), 'stable') in future versions (>= 2012a)
        oldIdx = seqOfEvents(:,3);
        [~, m] = unique(oldIdx, 'first');
        % mapping: oldIdx(sort(m)) -> 1:nSeg
        idxMap = oldIdx(sort(m));
        [~,newIdx] = ismember(oldIdx, idxMap);
        seqOfEvents(:,3) = newIdx;
        % replace parent indexes
        [~,newIdx] = ismember(seqOfEvents(:,4), idxMap);
        seqOfEvents(:,4) = newIdx;
        seqOfEvents(seqOfEvents(:,4)==0,4) = NaN;
        
        % re-assign to trackinfo, re-arrange with new index
        [~,ridx] = sort(idxMap);
        [~,lidx] = sort(ridx);
        trackinfo(i).seqOfEvents = seqOfEvents;
        trackinfo(i).tracksCoordAmpCG = tracksCoordAmpCG(lidx,:);
        trackinfo(i).tracksFeatIndxCG = tracksFeatIndxCG(lidx,:);
    end
end
end
%======================================================================



tracks(1:nTracks) = struct('t', [], 'f', [],...
    'x', [], 'y', [], 'A', [],...
    'tracksFeatIndxCG', [], 'gapVect', [], 'gapStatus', [], 'gapIdx', [], 'seqOfEvents', [],...
    'nSeg', [], 'visibility', [], 'lifetime_s', [], 'start', [], 'end', []);


%==============================
% Loop through tracks
%==============================
fprintf('Processing tracks (%s) - converting tracker output:     ', getShortPath(data));
for k = 1:nTracks
    
    % convert/assign structure fields
    seqOfEvents = trackinfo(k).seqOfEvents;
    tracksFeatIndxCG = trackinfo(k).tracksFeatIndxCG; % index of the feature in each frame
    nSeg = size(tracksFeatIndxCG,1);
    
    segLengths = NaN(1,nSeg);
    
    % discarding rules: single frame segments w/ merge/split. ADD: longer segments with merge & split
    msIdx = NaN(1,nSeg);
    for s = 1:nSeg
        idx = seqOfEvents(:,3)==s;
        ievents = seqOfEvents(idx, :);
        bounds = ievents(:,1); % beginning & end of this segment
        if ~isnan(ievents(2,4))
            bounds(2) = bounds(2)-1; % correction if end is a merge
        end
        segLengths(s) = bounds(2)-bounds(1)+1;
        
        % remove short (<4 frames) merging/splitting branches if:
        % -the segment length is a single frame
        % -the segment is splitting and merging from/to the same parent
        % -short segment merges, segment starts after track start
        % -short segment splits, segment ends before track end
        msIdx(s) = segLengths(s)==1 || (segLengths(s)<4 && ( diff(ievents(:,4))==0 ||...
            (isnan(ievents(1,4)) && ~isnan(ievents(2,4)) && ievents(1,1)>seqOfEvents(1,1)) ||...
            (isnan(ievents(2,4)) && ~isnan(ievents(1,4)) && ievents(2,1)<seqOfEvents(end,1)) ));
    end
    if preprocess
        if nSeg>1
            segIdx = find(msIdx==0); % segments to retain
            nSeg = numel(segIdx); % update segment #
            msIdx = find(msIdx);
            if ~isempty(msIdx)
                tracksFeatIndxCG(msIdx,:) = [];
                seqOfEvents(ismember(seqOfEvents(:,3), msIdx),:) = [];
            end
            segLengths = segLengths(segIdx);
        else
            segIdx = 1;
        end
    else
        segIdx = 1:nSeg;
    end
            
    tracks(k).nSeg = nSeg;
    firstIdx = trackinfo(k).seqOfEvents(1,1);
    lastIdx = trackinfo(k).seqOfEvents(end,1);
    
    tracks(k).lifetime_s = (lastIdx-firstIdx+1)*data.framerate;
    tracks(k).start = firstIdx;
    tracks(k).end = lastIdx;
    
    tracks(k).seqOfEvents = seqOfEvents;
    tracks(k).tracksFeatIndxCG = tracksFeatIndxCG; % index of the feature in each frame
    
    if (buffer(1)<tracks(k).start) && (tracks(k).end<=nFrames-buffer(2)) % complete tracks
        tracks(k).visibility = 1;
    elseif tracks(k).start==1 && tracks(k).end==nFrames % persistent tracks
        tracks(k).visibility = 3;
    else
        tracks(k).visibility = 2; % incomplete tracks
    end
    
    %==============================================================================
    % Initialize arrays
    %==============================================================================
    
    % Segments are concatenated into single arrays, separated by NaNs.
    fieldLength = sum(segLengths)+nSeg-1;
    tracks(k).x = NaN(1, fieldLength);
    tracks(k).y = NaN(1, fieldLength);
    tracks(k).A = NaN(1, fieldLength);
    tracks(k).t = NaN(1, fieldLength);
    tracks(k).f = NaN(1, fieldLength);
    
    %==============================================================================
    % Read amplitude & background from detectionResults.mat (localization results)
    %==============================================================================
    delta = [0 cumsum(segLengths(1:end-1))+(1:nSeg-1)];
    
    for s = 1:nSeg
        ievents = seqOfEvents(seqOfEvents(:,3)==segIdx(s), :);
        bounds = ievents(:,1);
        if ~isnan(ievents(2,4))
            bounds(2) = bounds(2)-1;
        end
        
        nf = bounds(2)-bounds(1)+1;
        frameRange = frameIdx(bounds(1):bounds(2)); % relative to movie (also when movie is subsampled)
        
        for i = 1:length(frameRange)
            idx = tracksFeatIndxCG(s, frameRange(i) - tracks(k).start + 1); % -> relative to IndxCG
            if idx ~= 0 % if not a gap, get detection values
                tracks(k).x(:,i+delta(s)) = frameInfo(frameRange(i)).xCoord(idx,1);
                tracks(k).y(:,i+delta(s)) = frameInfo(frameRange(i)).yCoord(idx,1);
                tracks(k).A(:,i+delta(s)) = frameInfo(frameRange(i)).amp(idx,1);
            end
        end
        tracks(k).t(delta(s)+(1:nf)) = (bounds(1)-1:bounds(2)-1)*data.framerate;
        tracks(k).f(delta(s)+(1:nf)) = frameRange;
    end
    
    fprintf('\b\b\b\b%3d%%', round(100*k/nTracks));
end
fprintf('\n');

% remove tracks that fall into image boundary
minx = round(arrayfun(@(t) min(t.x(:)), tracks));
maxx = round(arrayfun(@(t) max(t.x(:)), tracks));
miny = round(arrayfun(@(t) min(t.y(:)), tracks));
maxy = round(arrayfun(@(t) max(t.y(:)), tracks));

idx = minx<=w4 | miny<=w4 | maxx>nx-w4 | maxy>ny-w4;
tracks(idx) = [];
nTracks = numel(tracks);

%=======================================
% Interpolate gaps and clean up tracks
%=======================================
fprintf('Processing tracks (%s) - classification:     ', getShortPath(data));
for k = 1:nTracks
    
    % gap locations in 'x' for all segments
    gapVect = isnan(tracks(k).x(mCh,:)) & ~isnan(tracks(k).t);
    tracks(k).gapVect = gapVect;
    
    %=================================
    % Determine track and gap status
    %=================================
    sepIdx = isnan(tracks(k).t);
    
    gapCombIdx = diff(gapVect | sepIdx);
    gapStarts = find(gapCombIdx==1)+1;
    gapEnds = find(gapCombIdx==-1);
    gapLengths = gapEnds-gapStarts+1;
    
    segmentIdx = diff([0 ~(gapVect | sepIdx) 0]); % these variables refer to segments between gaps
    segmentStarts = find(segmentIdx==1);
    segmentEnds = find(segmentIdx==-1)-1;
    segmentLengths = segmentEnds-segmentStarts+1;
    
    % loop through track segments with gaps
    %if nanIdx{s}(1)==1 || nanIdx{s}(end)==1 % temporary fix for segments that begin or end with gap
    %    tracks(k).gapStatus{s} = 5;
    %else
    
    % loop over gaps
    nGaps = numel(gapLengths);
    if nGaps>0
        gv = 1:nGaps;
        gapStatus = 5*ones(1,nGaps);
        % gap valid if segments that precede/follow are > 1 frame or if gap is a single frame
        gapStatus(segmentLengths(gv)>1 & segmentLengths(gv+1)>1 | gapLengths(gv)==1) = 4;
        
        sepIdx = sepIdx(gapStarts)==1;
        gapStatus(sepIdx) = [];
        gapStarts(sepIdx) = [];
        gapEnds(sepIdx) = [];
        nGaps = numel(gapStatus);
        
        % fill position information for valid gaps using linear interpolation
        for g = 1:nGaps
            borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
            gacombIdx = gapStarts(g):gapEnds(g);
            for c = 1:nCh
                tracks(k).x(c, gacombIdx) = interp1(borderIdx, tracks(k).x(c, borderIdx), gacombIdx);
                tracks(k).y(c, gacombIdx) = interp1(borderIdx, tracks(k).y(c, borderIdx), gacombIdx);
                tracks(k).A(c, gacombIdx) = interp1(borderIdx, tracks(k).A(c, borderIdx), gacombIdx);
            end
        end
        tracks(k).gapStatus = gapStatus;
        tracks(k).gapIdx = arrayfun(@(i) gapStarts(i):gapEnds(i), 1:nGaps, 'UniformOutput', false);
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/nTracks));
end
fprintf('\n');

%====================================================================================
% Generate buffers before and after track, estimate gap values
%====================================================================================
% Gap map for fast indexing
gapMap = zeros(nTracks, data.movieLength);
for k = 1:nTracks
    gapMap(k, tracks(k).f(tracks(k).gapVect==1)) = 1;
end


%============================================================================
% Run post-processing
%============================================================================
if postprocess
%----------------------------------------------------------------------------
% I. Assign category to each track
%----------------------------------------------------------------------------
% Categories:
% Ia)  Single tracks with valid gaps
% Ib)  Single tracks with invalid gaps
% Ic)  Single tracks cut at beginning or end
% Id)  Single tracks, persistent
% IIa) Compound tracks with valid gaps
% IIb) Compound tracks with invalid gaps
% IIc) Compound tracks cut at beginning or end
% IId) Compound tracks, persistent

% The categories correspond to index 1-8, in the above order

validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
singleIdx = [tracks.nSeg]==1;
vis = [tracks.visibility];

mask_Ia = singleIdx & validGaps & vis==1;
mask_Ib = singleIdx & ~validGaps & vis==1;
idx_Ia = find(mask_Ia);
idx_Ib = find(mask_Ib);
trackLengths = [tracks.end]-[tracks.start]+1;

C = [mask_Ia;
     2*mask_Ib;
     3*(singleIdx & vis==2);
     4*(singleIdx & vis==3);
     5*(~singleIdx & validGaps & vis==1);
     6*(~singleIdx & ~validGaps & vis==1);
     7*(~singleIdx & vis==2);
     8*(~singleIdx & vis==3)];

C = num2cell(sum(C,1));
% assign category
[tracks.catIdx] = deal(C{:});    


%----------------------------------------------------------------------------
% III. Map gaps
%----------------------------------------------------------------------------
% Reference distribution: class Ia tracks
% Determine critical max. intensity values from class Ia tracks, per lifetime cohort

% # cohorts
nc = numel(cohortBounds)-1;

% max intensities of all 'Ia' tracks
maxInt = arrayfun(@(i) max(i.A(mCh,:)), tracks(idx_Ia));
maxIntDistr = cell(1,nc);
mappingThresholdMaxInt = zeros(1,nc);
lft_Ia = [tracks(idx_Ia).lifetime_s];
for i = 1:nc
    maxIntDistr{i} = maxInt(cohortBounds(i)<=lft_Ia & lft_Ia<cohortBounds(i+1));
    % critical values for test
    mappingThresholdMaxInt(i) = prctile(maxIntDistr{i}, 2.5);
end

% get lifetime histograms before change
processingInfo.lftHists.before = getLifetimeHistogram(data, tracks);

% Criteria for mapping:
% - max intensity must be within 2.5th percentile of max. intensity distribution for 'Ia' tracks
% - lifetime >= 5 frames (at 4 frames: track = [x o o x])

% assign category I to tracks that match criteria
for k = 1:numel(idx_Ib);
    i = idx_Ib(k);

    % get cohort idx for this track (logical)
    cIdx = cohortBounds(1:nc)<=tracks(i).lifetime_s & tracks(i).lifetime_s<cohortBounds(2:nc+1);
    
    if max(tracks(i).A(mCh,:)) >= mappingThresholdMaxInt(cIdx) && trackLengths(i)>4
        tracks(i).catIdx = 1;
    end
end

processingInfo.lftHists.after = getLifetimeHistogram(data, tracks);


% remove tracks with more gaps than frames
nGaps = arrayfun(@(i) sum(i.gapVect), tracks);
trackLengths = [tracks.end]-[tracks.start]+1;
[tracks(nGaps./trackLengths>=0.5).catIdx] = deal(2);

fprintf('Processing for %s complete - valid/total tracks: %d/%d (%.1f%%).\n',...
    getShortPath(data), sum([tracks.catIdx]==1), numel(tracks), sum([tracks.catIdx]==1)/numel(tracks)*100);

end % postprocessing


%==========================================
% Save results
%==========================================
if ~(exist([data.source 'Tracking'], 'dir')==7)
    mkdir([data.source 'Tracking']);
end
if isunix
    cmd = ['svn info ' mfilename('fullpath') '.m | grep "Last Changed Rev"'];
    [status,rev] = system(cmd);
    if status==0
        rev = regexp(rev, '\d+', 'match');
        processingInfo.revision = rev{1};
    end
end
processingInfo.procFlag = [preprocess postprocess];

save([data.source 'Tracking' filesep filename], 'tracks', 'processingInfo');
