function [ out ] = findCandidateIntersections( idxList, angleR )
%findCandidateIntersections Finds the local maxima along the the linear
%indices such as in PixelIdxList from bwconncomp
%
% indxList list of linear indices for 2D matrix
% angleR angle response Y x X x Theta, from steerable filter
%
% assumes the idxList is ordered such as by orderPixelList

angleR = shiftdim(angleR,2);
if(iscell(idxList))
    % vectorize if given a cell array
    out = cellfun(@findCandidateIntersections_scalar,idxList,'UniformOutput',false);
    % combine vertices by angle
    out = cellfun(@vertcat,out{:},'UniformOutput',false);
else
    % a simple cell array if an idxList vector given
    out = findCandidateIntersections_scalar(idxList);
end


    function outScalar = findCandidateIntersections_scalar(idxList)
        % get orientation space for idxList
        M = angleR(:,idxList);
        % zero the first value
        M = bsxfun(@minus,M,M(:,1));
        % smooth it out
%         K = normpdf(-3:3);
%         M = conv2(1,K,M,'same');

        maxima = M(:,2:end-1) > M(:,1:end-2) & M(:,2:end-1) > M(:,3:end);

        outScalar = cellfun(@(m) idxList([false m false]),num2cell(maxima,2),'Unif',false);
    end


end

