function check = trackPartitionResultCheck(tracksPart)
% Make sure all seqOfEvent's frame numbers are within the length of their
% tracks. Returns 1 if everything is good
badTrack = zeros(numel(tracksPart,1));
for i = 1:numel(tracksPart,1)
    trackFirstEvent = min(tracksPart(i).seqOfEvents(:,1));
    trackLastEvent = max(tracksPart(i).seqOfEvents(:,1));
    if (trackLastEvent-trackFirstEvent+1) ~= size(tracksPart(i).tracksFeatIndxCG,2)
        if trackLastEvent ~= size(tracksPart(i).trackFeatIndxCG,2)
            badTrack(i) = 1;
        end
    end
end
check = sum(badTrack) == 0;
end