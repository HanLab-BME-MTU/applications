function out = mergeShortTracks(intersection,minTrackLength)
sz = numel(intersection);
gapSize = minTrackLength;
out = intersection;
if sz <= 3*gapSize
    % If the original track is too short anyway
    out = ones(1,sz)*mode(intersection);
else
    parfor i = 1:sz
        beforeInd = max(gapSize+1,i-ceil(gapSize/2));
        afterInd = min(sz-gapSize,i+ceil(gapSize/2));
        before = intersection((beforeInd-gapSize):beforeInd);
        after = intersection(afterInd:(afterInd+gapSize));
        if intersection(i) ~= mode(after)
        end
        % If there are segments of the same value before and after this location
        if mode(before) == mode(after) 
            out(i) = mode([before,after]);
        else
            out(i) = intersection(i);
        end

    end
end



end