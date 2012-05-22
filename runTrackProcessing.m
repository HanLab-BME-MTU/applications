%[data] = runTrackProcessing(data, varargin)
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

function runTrackProcessing(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Buffer', 5, @isscalar);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackerOutput', 'trackedFeatures.mat', @ischar);
ip.addParamValue('FileName', 'ProcessedTracks.mat', @ischar);
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
for i = 1:length(data)
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

detection = load([data.source 'Detection' filesep 'detection_v2.mat']);
frameInfo = detection.frameInfo;

ny = data.imagesize(1);
nx = data.imagesize(2);
nFrames = length(frameIdx);

alpha = 0.05;
kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background

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
w2 = ceil(2*sigma);
w3 = ceil(3*sigma);
w4 = ceil(4*sigma);

[x,y] = meshgrid(-w4:w4);
r = sqrt(x.^2+y.^2);
annularMask = zeros(size(r));
annularMask(r<=w4 & r>=w3) = 1;

%======================================================================
% Read and convert tracker output
%======================================================================
tPath = [data.source 'Tracking' filesep trackerOutput];
if exist(tPath, 'file')==2
    trackinfo = load(tPath);
    trackinfo = trackinfo.tracksFinal;
    nTracks = length(trackinfo);
elseif exist([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat'], 'file')==2
    % (for old tracker. oldest version: trackInfo.mat)
    trackinfo = load([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat']);
    trackinfo = trackinfo.trackedFeatureInfo;
    nTracks = size(trackinfo, 1);
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
    'x', [], 'y', [], 'A', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'sigma_r', [], 'SE_sigma_r', [],...
    'pval_Ar', [], 'isPSF', [],...
    'tracksFeatIndxCG', [], 'gapVect', [], 'gapStatus', [], 'gapIdx', [], 'seqOfEvents', [],...
    'nSeg', [], 'visibility', [], 'lifetime_s', [], 'start', [], 'end', [],...
    'startBuffer', [], 'endBuffer', [], 'MotionParameters', []);
%    'alphaMSD', [], 'MSD', [], 'MSDstd', [], 'totalDisplacement', [], 'D', [], ...

% track field names
idx = structfun(@(i) size(i,2)==size(frameInfo(1).x,2), frameInfo(1));
mcFieldNames = fieldnames(frameInfo);
mcFieldNames = mcFieldNames(idx);
[~,loc] = ismember({'x_init', 'y_init', 'RSS'}, mcFieldNames);
mcFieldNames(loc) = [];

bufferFieldNames = {'t', 'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar'};

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
    
    if (buffer<tracks(k).start) && (tracks(k).end<=nFrames-buffer) % complete tracks
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
    for f = 1:length(mcFieldNames)
        tracks(k).(mcFieldNames{f}) = NaN(nCh, fieldLength);
    end
    tracks(k).t = NaN(1, fieldLength);
    tracks(k).f = NaN(1, fieldLength);
    
    if fieldLength>1
        
        % start buffer size for this track
        sb = firstIdx - max(1, firstIdx-buffer);
        eb = min(lastIdx+buffer, data.movieLength)-lastIdx;
        
        if sb>0 && tracks(k).visibility==1
            for f = 1:length(bufferFieldNames)
                tracks(k).startBuffer.(bufferFieldNames{f}) = NaN(nCh, sb);
            end
        end
        if eb>0 && tracks(k).visibility==1
            for f = 1:length(bufferFieldNames)
                tracks(k).endBuffer.(bufferFieldNames{f}) = NaN(nCh, eb);
            end
        end
    end
    
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
                for f = 1:length(mcFieldNames)
                    tracks(k).(mcFieldNames{f})(:,i+delta(s)) = frameInfo(frameRange(i)).(mcFieldNames{f})(:,idx);
                end
            end
        end
        tracks(k).t(delta(s)+(1:nf)) = (bounds(1)-1:bounds(2)-1)*data.framerate;
        tracks(k).f(delta(s)+(1:nf)) = frameRange;
    end
    %tracks(k).pval_Ar = 1-tracks(k).pval_Ar; % bug fix for detection
    
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

% for buffers:
trackStarts = [tracks.start];
trackEnds = [tracks.end];
fullTracks = [tracks.visibility]==1;

fprintf('Processing tracks (%s) - gap interpolation, buffer readout:     ', getShortPath(data));
for f = 1:data.movieLength
    
    mask = double(imread(data.maskPaths{f}));
    % binarize
    mask(mask~=0) = 1;
    labels = bwlabel(mask);
    
    for ch = 1:nCh
        frame = double(imread(data.framePaths{ch}{f}));
        
        %------------------------
        % Gaps
        %------------------------
        % tracks with valid gaps visible in current frame
        currentGapsIdx = find(gapMap(:,f));
        for ki = 1:numel(currentGapsIdx)
            k = currentGapsIdx(ki);
            
            % index in the track structure (.x etc)
            idxList = find(tracks(k).f==f & tracks(k).gapVect==1);
            
            for l = 1:numel(idxList)
                idx = idxList(l);
                xi = round(tracks(k).x(ch,idx));
                yi = round(tracks(k).y(ch,idx));
                
                % window/masks (see psfLocalization.m for details)
                maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
                
                ci = mean(window(cmask==1));
                window(maskWindow~=0) = NaN;
                
                x0 = tracks(k).x(ch,idx)-xi;
                y0 = tracks(k).y(ch,idx)-yi;
                npx = sum(isfinite(window(:)));
                [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
                dx = prm(1);
                dy = prm(2);
                if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2)
                    tracks(k).x(ch,idx) = xi+dx;
                    tracks(k).y(ch,idx) = yi+dy;
                    tracks(k).A_pstd(ch,idx) = prmStd(3);
                    tracks(k).c_pstd(ch,idx) = prmStd(4);
                else
                    [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                    tracks(k).A_pstd(ch,idx) = prmStd(1);
                    tracks(k).c_pstd(ch,idx) = prmStd(2);
                end
                tracks(k).A(ch,idx) = prm(3);
                tracks(k).c(ch,idx) = prm(5);
                
                tracks(k).sigma_r(ch,idx) = res.std;
                tracks(k).SE_sigma_r(ch,idx) = res.std/sqrt(2*(npx-1));
                
                SE_r = tracks(k).SE_sigma_r(ch,idx) * kLevel;
                
                tracks(k).hval_AD(idx) = res.hAD;
                
                df2 = (npx-1) * (tracks(k).A_pstd(ch,idx).^2 + SE_r.^2).^2 ./...
                    (tracks(k).A_pstd(ch,idx).^4 + SE_r.^4);
                
                scomb = sqrt((tracks(k).A_pstd(ch,idx).^2 + SE_r.^2)/npx);
                T = (tracks(k).A(ch,idx) - res.std*kLevel) ./ scomb;
                tracks(k).pval_Ar(ch,idx) = tcdf(-T, df2);
            end
        end
        
        %------------------------
        % start buffer
        %------------------------
        % tracks with start buffers in this frame
        cand = max(1, trackStarts-buffer)<=f & f<trackStarts;
        % corresponding tracks, only if status = 1
        currentBufferIdx = find(cand & fullTracks);
        
        for ki = 1:length(currentBufferIdx)
            k = currentBufferIdx(ki);
            
            xi = round(tracks(k).x(ch,1));
            yi = round(tracks(k).y(ch,1));
            
            % window/masks (see psfLocalization.m for details)
            maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
            maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
            ci = mean(window(cmask==1));
            window(maskWindow~=0) = NaN;
            
            x0 = tracks(k).x(1)-xi;
            y0 = tracks(k).y(1)-yi;
            
            [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
            bi = f - max(1, tracks(k).start-buffer) + 1;
            dx = prm(1);
            dy = prm(2);
            if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                tracks(k).startBuffer.x(ch,bi) = xi+dx;
                tracks(k).startBuffer.y(ch,bi) = yi+dy;
                tracks(k).startBuffer.A_pstd(ch,bi) = prmStd(3);
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).startBuffer.x(ch,bi) = tracks(k).x(ch,1);
                tracks(k).startBuffer.y(ch,bi) = tracks(k).y(ch,1);
                tracks(k).startBuffer.A_pstd(ch,bi) = prmStd(1);
            end
            tracks(k).startBuffer.A(ch,bi) = prm(3);
            tracks(k).startBuffer.c(ch,bi) = prm(5);
            tracks(k).startBuffer.sigma_r(ch,bi) = res.std;
            
            npx = sum(isfinite(window(:)));
            
            tracks(k).startBuffer.SE_sigma_r(ch,bi) = res.std/sqrt(2*(npx-1));
            SE_r = tracks(k).startBuffer.SE_sigma_r(ch,bi) * kLevel;
            
            df2 = (npx-1) * (tracks(k).startBuffer.A_pstd(ch,bi).^2 + SE_r.^2).^2 ./...
                (tracks(k).startBuffer.A_pstd(ch,bi).^4 + SE_r.^4);
            
            scomb = sqrt((tracks(k).startBuffer.A_pstd(ch,bi).^2 + SE_r.^2)/npx);
            T = (tracks(k).startBuffer.A(ch,bi) - res.std*kLevel) ./ scomb;
            tracks(k).startBuffer.pval_Ar(ch,bi) = tcdf(-T, df2);
        end
        
        %------------------------
        % end buffer
        %------------------------
        % segments with end buffers in this frame
        cand = trackEnds<f & f<=min(data.movieLength, trackEnds+buffer);
        % corresponding tracks
        currentBufferIdx = find(cand & fullTracks);
        
        for ki = 1:length(currentBufferIdx)
            k = currentBufferIdx(ki);
            xi = round(tracks(k).x(ch,end));
            yi = round(tracks(k).y(ch,end));
            
            % window/masks (see psfLocalization.m for details)
            maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
            maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
            ci = mean(window(cmask==1));
            window(maskWindow~=0) = NaN;
            
            x0 = tracks(k).x(end)-xi;
            y0 = tracks(k).y(end)-yi;
            
            [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
            bi = f - tracks(k).end;
            dx = prm(1);
            dy = prm(2);
            if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                tracks(k).endBuffer.x(ch,bi) = xi+dx;
                tracks(k).endBuffer.y(ch,bi) = yi+dy;
                tracks(k).endBuffer.A_pstd(ch,bi) = prmStd(3);
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).endBuffer.x(ch,bi) = tracks(k).x(ch,end);
                tracks(k).endBuffer.y(ch,bi) = tracks(k).y(ch,end);
                tracks(k).endBuffer.A_pstd(ch,bi) = prmStd(1);
            end
            tracks(k).endBuffer.A(ch,bi) = prm(3);
            tracks(k).endBuffer.c(ch,bi) = prm(5);
            tracks(k).endBuffer.sigma_r(ch,bi) = res.std;
            
            npx = sum(isfinite(window(:)));
            
            tracks(k).endBuffer.SE_sigma_r(ch,bi) = res.std/sqrt(2*(npx-1));
            SE_r = tracks(k).endBuffer.SE_sigma_r(ch,bi) * kLevel;
            
            df2 = (npx-1) * (tracks(k).endBuffer.A_pstd(ch,bi).^2 + SE_r.^2).^2 ./...
                (tracks(k).endBuffer.A_pstd(ch,bi).^4 + SE_r.^4);
            
            scomb = sqrt((tracks(k).endBuffer.A_pstd(ch,bi).^2 + SE_r.^2)/npx);
            T = (tracks(k).endBuffer.A(ch,bi) - res.std*kLevel) ./ scomb;
            tracks(k).endBuffer.pval_Ar(ch,bi) = tcdf(-T, df2);
        end
        fprintf('\b\b\b\b%3d%%', round(100*(ch + (f-1)*nCh)/(nCh*data.movieLength)));
    end
end
fprintf('\n');


% %==========================================
% % Compute displacements
% %==========================================
% % Only on tracks with no/valid gaps
% trackIdx = find(arrayfun(@(t) isempty([t.gapStatus{:}]) || max([t.gapStatus{:} 4])==4, tracks));
% fprintf('TrackProcessing - Properties:     ');
% for ki = 1:length(trackIdx)
%     k = trackIdx(ki);
%     ns = numel(tracks(k).x);
%     msdVect = cell(1,ns);
%     msdStdVect = cell(1,ns);
%     for s = 1:ns
%
%         x = tracks(k).x{s}(mCh,:);
%         y = tracks(k).y{s}(mCh,:);
%         tracks(k).totalDisplacement{s} = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
%         % MSD
%         L = 10;
%         msdVect{s} = NaN(1,L);
%         msdStdVect{s} = NaN(1,L);
%         for l = 1:min(L, numel(x)-1)
%             tmp = (x(1+l:end)-x(1:end-l)).^2 + (y(1+l:end)-y(1:end-l)).^2;
%             msdVect{s}(l) = mean(tmp);
%             msdStdVect{s}(l) = std(tmp);
%         end
%         tracks(k).MSD = msdVect;
%         tracks(k).MSDstd = msdStdVect;
%
%         %if L > 1 % min 2 points to fit
%         %    [D c alpha] = fitMSD(MSDvect(1:L), [MSDvect(L)/(4*L) 0 1], 'Dc');
%         %    tracks(k).D = D;
%         %    tracks(k).c = c;
%         %    tracks(k).alpha = alpha;
%         %end
%
%         % add buffer time vectors
%         %b = size(tracks(k).startBuffer.x,2);
%         %tracks(k).startBuffer.t = ((-b:-1) + tracks(k).segmentStarts(s)-1) * data.framerate;
%         %b = size(tracks(k).endBuffer.x,2);
%         %tracks(k).endBuffer.t = (tracks(k).segmentEnds(s) + (1:b)-1) * data.framerate;
%     end
%     fprintf('\b\b\b\b%3d%%', round(100*k/nTracks));
% end
% fprintf('\n');

for k = 1:nTracks
    % add buffer time vectors
    if ~isempty(tracks(k).startBuffer)
        b = size(tracks(k).startBuffer.x,2);
        tracks(k).startBuffer.t = ((-b:-1) + tracks(k).start-1) * data.framerate;
    end
    if ~isempty(tracks(k).endBuffer)
        b = size(tracks(k).endBuffer.x,2);
        tracks(k).endBuffer.t = (tracks(k).end + (1:b)-1) * data.framerate;
    end
end

% save([data.source 'Tracking' filesep 'tracksTMP.mat'], 'tracks');
% load([data.source 'Tracking' filesep 'tracksTMP.mat']);


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
% II. Identify diffraction-limited tracks (CCPs)
%----------------------------------------------------------------------------
% Criterion: if all detected points pass AD-test, then track is a CCP.
% (gaps in the track are not considered in this test)

% # diffraction-limited points per track (can be different from track length for compound tracks!)
nPl = arrayfun(@(i) nansum(i.hval_AD(mCh,:) .* ~i.gapVect), tracks);

isCCP = num2cell(nPl==0);
[tracks.isCCP] = deal(isCCP{:});
isCCP = [isCCP{:}];

% average mask area per track
% meanMaskAreaCCP = arrayfun(@(i) nanmean(i.maskN), tracks(isCCP));
% meanMaskAreaNotCCP = arrayfun(@(i) nanmean(i.maskN), tracks(~isCCP));

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

% v = hist([tracks.catIdx], 1:8);
% v = v/numel(tracks);
% plotTrackClasses(v');
%----------------------------------------------------------------------------
% IV. Apply threshold on buffer intensities
%----------------------------------------------------------------------------
% Condition: at least 2 frames must be below background noise level in both start and end buffer
Tbuffer = 2;


% loop through cat. Ia tracks
idx_Ia = find([tracks.catIdx]==1);
for k = 1:numel(idx_Ia)
    i = idx_Ia(k);
    
    % H0: A = background (p-value >= 0.05) 
    if sum(tracks(i).startBuffer.pval_Ar(mCh,:)>=0.05) < Tbuffer ||...
            sum(tracks(i).endBuffer.pval_Ar(mCh,:)>=0.05) < Tbuffer
        tracks(i).catIdx = 2;
    end
end
% v = hist([tracks.catIdx], 1:8);
% v = v/numel(tracks);
% plotTrackClasses(v');

%----------------------------------------------------------------------------
% V. Assign Cat. Ib to tracks that are not diffraction-limited CCPs
%----------------------------------------------------------------------------
[tracks([tracks.catIdx]==1 & ~isCCP).catIdx] = deal(2);

%----------------------------------------------------------------------------
% VI. Cut tracks with sequential events (hotspots) into individual tracks
%----------------------------------------------------------------------------
splitCand = find([tracks.catIdx]==1 & arrayfun(@(i) ~isempty(i.gapIdx), tracks) & trackLengths>4);

% Loop through tracks and test whether gaps are at background intensity
rmIdx = []; % tracks to remove from list after splitting
newTracks = [];
for i = 1:numel(splitCand);
    k = splitCand(i); 
    
    % all gaps
    gapIdx = [tracks(k).gapIdx{:}];
    
    % # residual points
    npx = round((tracks(k).sigma_r(mCh,:) ./ tracks(k).SE_sigma_r(mCh,:)).^2/2+1);
    npx = npx(gapIdx);
    
    % perform t-test on all gaps
    A = tracks(k).A(mCh, gapIdx);
    sigma_A = tracks(k).A_pstd(mCh, gapIdx);
    %sigma_c = tracks(k).A_pstd(mCh, gapIdx);
    
    %SE_sigma_r = tracks(k).SE_sigma_r(gapIdx) * kLevel;
    %sigma_r = tracks(k).sigma_r(gapIdx) * kLevel;
    
    % df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
    % scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
    % % H1: sigma_r > A
    % T = (sigma_r - A) ./ scomb;
    % pval = tcdf(-T, df2);
    
    %T = A./(sigma_A./sqrt(npx));
    %pval = tcdf(T, npx-1);
    
    T = (A-sigma_A)./(sigma_A./sqrt(npx)); % instead of comparing against 0: 0+sigma_A
    pval = tcdf(T, npx-1);
    
    %df2 = (npx-1) .* (sigma_A.^2 + sigma_c.^2).^2 ./ (sigma_A.^4 + sigma_c.^4);
    %scomb = sqrt((sigma_A.^2 + sigma_c.^2)./npx);
    % % H1: sigma_r > A
    %T = ( - A) ./ scomb;
    %pval = tcdf(-T, df2);
    
    % gaps with signal below background level: candidates for splitting
    splitIdx = pval<0.05;
    gapIdx = gapIdx(splitIdx==1);
    
    % new segments must be at least 5 frames
    delta = diff([1 gapIdx trackLengths(k)]);
    gapIdx(delta(1:end-1)<5 | delta(2:end)<5) = [];
    
    ng = numel(gapIdx);
    splitIdx = zeros(1,ng);
    
    for g = 1:ng

        % split track at gap position
        x1 = tracks(k).x(mCh, 1:gapIdx(g)-1);
        y1 = tracks(k).y(mCh, 1:gapIdx(g)-1);
        x2 = tracks(k).x(mCh, gapIdx(g)+1:end);
        y2 = tracks(k).y(mCh, gapIdx(g)+1:end);
        mux1 = median(x1);
        muy1 = median(y1);
        mux2 = median(x2);
        muy2 = median(y2);
        
        % projections
        v = [mux2-mux1; muy2-muy1];
        v = v/norm(v);
        
        % x1 in mux1 reference
        X1 = [x1-mux1; y1-muy1];
        sp1 = sum(repmat(v, [1 numel(x1)]).*X1,1);
        
        % x2 in mux1 reference
        X2 = [x2-mux1; y2-muy1];
        sp2 = sum(repmat(v, [1 numel(x2)]).*X2,1);
        
        % test whether projections are distinct distributions of points
        % may need to be replaced by outlier-robust version
        %splitIdx(g) = ttest2(sp1, sp2, 0.01);
        if mean(sp1)<mean(sp2) && prctile(sp1,95)<prctile(sp2,5)
            splitIdx(g) = 1;
        elseif mean(sp1)>mean(sp2) && prctile(sp1,5)>prctile(sp2,95)
            splitIdx(g) = 1;
        else
            splitIdx(g) = 0;
        end
    end    
    gapIdx = gapIdx(splitIdx==1);
    
    if ~isempty(gapIdx)
        rmIdx = [rmIdx k]; % store index of parent track, to be removed at end

        % new tracks
        splitTracks = cutTrack(tracks(k), gapIdx);
        newTracks = [newTracks splitTracks];
    end
end
% final assignment
% fprintf('# tracks cut: %d\n', numel(rmIdx));
tracks(rmIdx) = [];
tracks = [tracks newTracks];

% remove tracks with more gaps than frames
nGaps = arrayfun(@(i) sum(i.gapVect), tracks);
trackLengths = [tracks.end]-[tracks.start]+1;

% fprintf('# tracks with >50%% gaps: %d\n', sum(nGaps./trackLengths>=0.5));
[tracks(nGaps./trackLengths>=0.5).catIdx] = deal(2);

% v = hist([tracks.catIdx], 1:8);
% v = v/numel(tracks);
% plotTrackClasses(v');

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



% function [D c alpha] = fitMSD(MSD, prmVect, prmSel)
%
% opts = optimset('Jacobian', 'off', ...
%     'MaxFunEvals', 1e4, ...
%     'MaxIter', 1e4, ...
%     'Display', 'off', ...
%     'TolX', 1e-6, ...
%     'Tolfun', 1e-6);
%
%
% estIdx = false(1,3); % [D c alpha]
% estIdx(regexp('Dca', ['[' prmSel ']'])) = true;
%
% x = lsqnonlin(@costMSD, prmVect(estIdx), [], [], opts, MSD, prmVect, estIdx);
% prmVect(estIdx) = deal(abs(x));
%
% D = prmVect(1);
% c = prmVect(2); % constant offset
% alpha = prmVect(3);


% function [v] = costMSD(p, MSD, prmVect, estIdx)
% prmVect(estIdx) = deal(abs(p));
%
% D = prmVect(1);
% c = prmVect(2); % constant offset
% alpha = prmVect(3);
% %d % dimensionality
%
% t = 1:length(MSD);
%
% v = MSD - 4*D*t.^alpha - c;