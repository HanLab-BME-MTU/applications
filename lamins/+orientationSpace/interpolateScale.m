function orientationSpaceAtScale = interpolateScale(orientationScaleSpace,scale)
% Interpolate scale
%
% INPUT
% orientationScaleSpace -  orientatation x scale x spatial dimensions
% scale - scale at which to find each pixel
%
% OUTPUT
% orientationSpace - at scale

% Rearrange orientationScaleSpace to be scale x orientation x space
orientationScaleSpace = permute(orientationScaleSpace,[2 1 3]);
% Mirror scale into Chebyshev coordinates
orientationScaleSpace = orientationScaleSpace([end:-1:1 2:end-1],:,:);
% Replicate scale across the the orientation dimension
scale = shiftdim(repmat(scale,[size(orientationScaleSpace,2) 1]),-1);
% Use Horner's method to evaluate scale at a particular point
orientationSpaceAtScale = interpft1([0 2*pi],orientationScaleSpace,scale,'horner',NaN);
%Squeeze so that output is orientation x space
orientationSpaceAtScale = squeeze(orientationSpaceAtScale);
    
end
