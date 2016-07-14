function out = mergeShortTracks(intersection,minTrackLength)
%MERGESHORTRACKS Merge together partition tracks that are shorter than a
%specified length

nFrames = numel(intersection);
gapSize = minTrackLength;
out = intersection;
if nFrames <= 3*gapSize
    % If the whole track is too short, just output a vector equal to the
    % mode of its values
    out = ones(1,nFrames)*mode(intersection);
else
    % Note: this loop probably won't parallelize since trackPartition uses
    % it within another parallelized operation
    for i = 1:nFrames
        beforeInd = max(gapSize+1,i-floor(gapSize/2));
        afterInd = min(nFrames-gapSize,i+floor(gapSize/2));
        before = intersection((beforeInd-gapSize):beforeInd);
        after = intersection(afterInd:(afterInd+gapSize));
        
        if mode(before) == mode(after) 
            % If there are sections with the same mode before and after this
            % frame, then this frame should equal that value
            % e.g. 1 1 1 1...[100]...1 1 1 1 -> 1 1 1 1...[1]...1 1 1 1
            % Use the mode to account for cases like this:
            % 1 1 1 1 [100] 1 100 1 1 -> 1 1 1 1 [1] 1 100 1 1
            out(i) = mode(before);
        else
            % If the sections before and after this frame have different
            % modes, then this frame may be on an edge
            % e.g. 1 1 1 1 [1] 0 0 0 0 
            % But it could also be something like this:
            % 1 1 1 1...[1]...0 0 0 1 ([0 0 0] constitutes another short  
            % track that should be removed later on)
            
            % Check for the second case by examining sections twice as long
            % before and after this frame
            beforeInd2 = max(2*gapSize+1,i-floor(gapSize/2));
            afterInd2 = min(nFrames-2*gapSize,i+floor(gapSize/2));
            before2 = intersection((beforeInd2-2*gapSize):beforeInd2);
            after2 = intersection(afterInd2:(afterInd2+2*gapSize));
            if mode(before2) == mode(after2)
                % Not actually an edge frame
                out(i) = mode(before2);
            else
                % Is actually an edge frame
                out(i) = intersection(i);
            end 
        end
    end
end
end