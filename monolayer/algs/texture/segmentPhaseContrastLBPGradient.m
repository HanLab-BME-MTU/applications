%% Input:
% I - input image
% patchSize - size of patches
% centers for histogram
% 
% Output: binary image at original resolution

function [BIN] = segmentPhaseContrastLBPGradient(I,lpbResult)

% get histogram for each patch in the image
mask = [[1 2 1];[0 0 0];[-1 -2 -1]];
Deriv = highDimensionDerivative(lpbResult.IHist,mask);

BIN = imresize(Deriv,size(I));
end

function [DERIV] = highDimensionDerivative(I,mask)
[sizeY,sizeX,nDim] = size(I);
allDeriv  = nan(sizeY,sizeX,nDim);
for d = 1 : nDim
    Id = I(:,:,d);
    energy1 = abs(imfilter(Id,mask,'same'));
    energy2 = abs(imfilter(Id,mask','same'));
    allDeriv(:,:,d) = sqrt(energy1.^2+energy2.^2);
end
DERIV = sum(allDeriv,3);
end