function out = mergeShortTracks(intersection,minTrackLength)
%MERGESHORTRACKS Merge together partition inside/outside segments that are 
%shorter than a specified length. E.g., if a particle is inside for only a
%few frames before moving back outside, get rid of that short inside
%segment
if minTrackLength == 1
    % minimum length of 1 doesn't do anything
    return
end
nFrames = numel(intersection);
out = intersection;
if nFrames <= minTrackLength
    % If the whole track is too short, just output a vector equal to the
    % mode of its values
    out = ones(1,nFrames)*mode(intersection);
else
    % Start with merging smaller lengths and work up, so that
    % [out out out in in out in in out]
    % doesn't turn into 
    % [out out out out out out out out out]
    % By starting with smaller lengths, this example would be correctly
    % merged into
    % [out out out in in in in in out]
    for l = 2:minTrackLength
        remaining = true;
        % The first segment shouldn't be merged with anything
        startFrame = find(out,1,'first');
        startValue = out(startFrame);
        if startValue == 1
            startFrame = find(out(startFrame:end) == 100,1,'first')+startFrame-1;
        else
            startFrame = find(out(startFrame:end) == 1,1,'first')+startFrame-1;
        end
        if isempty(startFrame)
            remaining = false;
        end
        while remaining
            startValue = out(startFrame);
            % Find where this inside/outside segment stops
            if startValue == 1
                endFrame = find(out(startFrame:end) == 100,1,'first')+startFrame-1;
            else
                endFrame = find(out(startFrame:end) == 1,1,'first')+startFrame-1;
            end
            if isempty(endFrame)
                % End of track has been reached
                remaining = false;
            elseif (endFrame-startFrame < l) && (out(endFrame) == out(startFrame-1))
                % If the segment is too short AND it is surrounded by
                % segments of the opposite type, change this segment to
                % their type
                out(startFrame:endFrame-1) = out(startFrame-1);
            end
            startFrame = endFrame;
        end
    end
end
end