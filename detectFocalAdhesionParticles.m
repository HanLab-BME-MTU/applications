function [x y amp] = detectFocalAdhesionParticles(ima, mask, sigmaPSF, kSigma)

ima = double(ima);
[nrows ncols] = size(ima);

% Filter image with laplacian
bandPassIso = filterLoG(ima,sigmaPSF);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(bandPassIso, [5 5]);

bw = blobSegmentThreshold(bandPassIso,0,0,mask);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(ima), indMax);
amp = ima(indMax);

% Subresolution detection
radius = kSigma * sigmaPSF;
hside = ceil(radius);
xmin = x - hside;
xmax = x + hside;
ymin = y - hside;
ymax = y + hside;

validFeaturesIdx = find(xmin >= 1 & xmax <= ncols & ymin >= 1 & ymax <= nrows);

[X,Y] = meshgrid(-hside:hside);
disk = X.^2 + Y.^2 - (kSigma * sigmaPSF)^2 <= 0;



for iiFeature = 1:numel(validFeaturesIdx)
    
    iFeature = validFeaturesIdx(iiFeature);
    
    crop = ima(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    crop(~disk) = NaN;
    
    [params stdParams] = fitGaussian2D(crop, [0, 0, amp(iFeature), sigmaPSF, min(crop(:))], 'xyAC');
    
    if stdParams(1) >= radius || stdParams(2) 
    
    x(iFeature) = x(iFeature) + params(1);
    y(iFeature) = y(iFeature) + params(2);
    amp(iFeature) = params(3);
end
