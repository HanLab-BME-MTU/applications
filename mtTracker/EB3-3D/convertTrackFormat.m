function tracks=convertTrackFormat(trackInfo,movieInfo)
% A more readable structures for tracks 
% Philippe Roudot 2014. After F. Aguet 

preprocess=0;
nTracks=length(trackInfo);

% Set up track structure
tracks(1:nTracks) = struct('t', [], 'f', [],...
                           'x', [], 'y', [],'z',[], 'A', [],...
                           'lifetime',[]...
                           );

frameIdx=1:length(movieInfo);           % No option for redudction now.

%==============================
% Loop through tracks
%==============================
fprintf('Processing tracks (%s) - converting tracker output:     ');
for k = 1:nTracks
    % convert/assign structure fields
    seqOfEvents = trackInfo(k).seqOfEvents;
    tracksFeatIndxCG = trackInfo(k).tracksFeatIndxCG; % index of the feature in each frame
    nSeg = size(tracksFeatIndxCG,1);
    
    segLengths = NaN(1,nSeg);
    
    % Remove short merging/splitting branches
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
    if preprocess && nSeg>1
        segIdx = find(msIdx==0); % index segments to retain (avoids re-indexing segments)
        nSeg = numel(segIdx); % update segment #
        msIdx = find(msIdx);
        if ~isempty(msIdx)
            tracksFeatIndxCG(msIdx,:) = [];
            seqOfEvents(ismember(seqOfEvents(:,3), msIdx),:) = [];
        end
        segLengths = segLengths(segIdx);
    else
        segIdx = 1:nSeg;
    end
    
    tracks(k).nSeg = nSeg;
    firstIdx = trackInfo(k).seqOfEvents(1,1);
    lastIdx = trackInfo(k).seqOfEvents(end,1);
    tracks(k).lifetime = (lastIdx-firstIdx+1);
    tracks(k).start = firstIdx;
    tracks(k).end = lastIdx; 
    
    tracks(k).seqOfEvents = seqOfEvents;
    tracks(k).tracksFeatIndxCG = tracksFeatIndxCG; 

    % index of the feature in each frame

    %==============================================================================
    % Initialize arrays
    %==============================================================================

    % Segments are concatenated into single arrays, separated by NaNs.
    % !!!
    mFieldNames={'x','y','z','A'};
    mFieldSizes = structfun(@(i) size(i,1), movieInfo(1));

    fieldLength = sum(segLengths)+nSeg-1;
    for f = 1:length(mFieldNames)
        tracks(k).(mFieldNames{f}) = nan(1, fieldLength);
    end

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
                tracks(k).x(:,i+delta(s)) = movieInfo(frameRange(i)).xCoord(idx,1);
                tracks(k).y(:,i+delta(s)) = movieInfo(frameRange(i)).yCoord(idx,1);
                tracks(k).z(:,i+delta(s)) = movieInfo(frameRange(i)).zCoord(idx,1);
                tracks(k).A(:,i+delta(s)) = movieInfo(frameRange(i)).amp(idx,1);
            end
        end
        tracks(k).t(delta(s)+(1:nf)) = (bounds(1)-1:bounds(2)-1);%*data.framerate;
        tracks(k).f(delta(s)+(1:nf)) = frameRange;
    end

    fprintf('\b\b\b\b%3d%%', round(100*k/nTracks));
end
fprintf('\n');