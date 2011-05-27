function [normWeights,listNormWeights]=getNormWeights(forceMesh)
weights    =vertcat(forceMesh.basis(:).unitVolume); % volume of the basis function
repWeights =repmat(weights(:),2,1); % the basis function for x/y comp. have the same weight
normWeights=repWeights/max(repWeights); % normalize it with max value.

if nargout>1
    listNormWeights=unique(normWeights);
end

