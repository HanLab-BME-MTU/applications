function [ orientationAtScale ] = getResponseAtScale( I, scale, angularOrder , npts)
%getResponseAtScale Get orientation space at a different scale per pixel at
%a particular angularOrder
%
% INPUT
% I - image
% scale - scale map the size of I
% angularOrder - scalar value indicating angular order to use
% npts - number of points to use with Chebyshev-Lobatto grid
%
% OUTPUT
% orientationAtScale - orientationSpace matrix at the scale specified


    if(nargin < 4)
        npts = 5;
    end
    mask = scale > 0;
    maskedScale = scale(mask);
    minScale = min(maskedScale);
    maxScale = max(maskedScale);
    scales = minScale:min(maxScale-minScale,1):maxScale;
    nAngles = 2*angularOrder+1;
    orientationAtScale = zeros([nAngles size(I)]);

    for s = 1:length(scales)-1
        if(s == length(scales) - 1)
            inRange = scale >= scales(s) & scale <= scales(s+1);
        else
            inRange = scale >= scales(s) & scale < scales(s+1);
        end
        if(sum(inRange(:)) == 0)
            continue;
        end
        % Scan scales along the Chebyshev-Lobatto grid between pairs of
        % scale values
        F = OrientationSpaceRidgeFilter(1./(2*pi)./chebpts(npts,scales([0 1]+s)),[],angularOrder);
        R = F*I;
        orientationScaleSpace = permute(R.getArraySpace,[3 4 1 2]);
        orientationScaleSpace = orientationScaleSpace(:,:,inRange);

        clear F R A
        interval = scales(s+1) - scales(s);
        % map to be within range [-interval/2 interval/2]
        mappedScales = scale(inRange) - scales(s) - interval/2;
        % map to be within range [-1 1]
        mappedScales = mappedScales/interval*2;
        % map to be [0 pi]
        mappedScales = acos(mappedScales);
        orientationAtScale(:,inRange) = orientationSpace.interpolateScale(orientationScaleSpace,mappedScales');
    end
    
    orientationAtScale = permute(orientationAtScale,[2 3 1]);
end