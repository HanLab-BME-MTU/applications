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
        if out(i) ~= 0 % don't fill frames where track doesn't exist
%             beforeInd = max(gapSize+1,i-floor(gapSize/2));
%             afterInd = min(nFrames-gapSize,i+floor(gapSize/2));
%             before = intersection((beforeInd-gapSize):beforeInd);
%             after = intersection(afterInd:(afterInd+gapSize));

            beforeInd = i-floor(gapSize/2);
            afterInd = i+floor(gapSize/2);
            before = intersection(max(beforeInd-gapSize,1):max(beforeInd,1));
            after = intersection(min(afterInd,nFrames):min(afterInd+gapSize,nFrames));

            if ((mode(before)==1)&&(mode(after)==1)) || ((mode(before)==100)&&(mode(after)==100))
                % If there are sections with the same mode before and after this
                % frame, then this frame should equal that value
                % e.g. 1 1 1 1...[100]...1 1 1 1 -> 1 1 1 1...[1]...1 1 1 1
                % Use the mode to account for cases like this:
                % 1 1 1 1 [100] 1 100 1 1 -> 1 1 1 1 [1] 1 100 1 1
                out(i) = mode(before);
            else
                % If the sections before and after this frame have different
                % modes, then this frame may be on an edge
                % e.g. 1 1 1 1 [1] 100 100 100 100
                % But it could also be something like this:
                % 1 1 1 1...[1]...100 100 100 1 ([100 100 100] constitutes another short  
                % track that should be removed later on)

                % Check for the second case by examining sections twice as long
                % before and after this frame
%                 beforeInd2 = max(2*gapSize+1,i-floor(gapSize/2));
%                 afterInd2 = min(nFrames-2*gapSize,i+floor(gapSize/2));
%                 before2 = intersection((beforeInd2-2*gapSize):beforeInd2);
%                 after2 = intersection(afterInd2:(afterInd2+2*gapSize));

%                 beforeInd2 = i-floor(gapSize/2);
%                 afterInd2 = i+floor(gapSize/2);
                before2 = intersection(max(beforeInd-2*gapSize,1):max(beforeInd,1));
                after2 = intersection(min(afterInd,nFrames):min(afterInd+2*gapSize,nFrames));
                
                if ((mode(before2)==1)&&(mode(after2)==1)) || ((mode(before2)==100)&&(mode(after2)==100))
                    % Not actually an edge frame
                    out(i) = mode(before2);
                else
                    % Check the immediate neighbors of this frame
                    before3 = intersection(max(beforeInd,1):max(i-1,1));
                    after3 = intersection(min(afterInd,nFrames):min(i+1,nFrames));
                    if ((mode(before3)==1)&&(mode(after3)==1)) || ((mode(before3)==100)&&(mode(after3)==100))
                        out(i) = mode(before3);
                    else
                        % Is actually an edge frame
                        out(i) = intersection(i);
                    end
                end 
            end
        end
    end
end
end