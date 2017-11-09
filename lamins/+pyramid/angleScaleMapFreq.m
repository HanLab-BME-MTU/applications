function [ map ] = angleScaleMapFreq( I_hat, angularFilters, scaleFilters , show)
%ANGLESCALEMAP Produces a map that represents the overall tiling of
%frequency space by angular and scale filters
%
% INPUT
% I_hat - image data transformed by fft2, X by Y
% angularFilters - frequency filter over angular space, X by Y by A
% scaleFilters - frequency filter over scale space (radial frequency), X by Y
%                by S
%
% OUTPUT
% map - tile map of frequency space in complex form, A by S

if(nargin < 2 || isempty(angularFilters))
    angularFilters = orientationSpace.angularKernel(5,[],size(I_hat));
end
if(nargin < 3 || isempty(scaleFilters))
    scaleFilters = pyramid.multiplierFilterKernel(size(I_hat));
end
    

P = numel(I_hat);
angularFilters = reshape(angularFilters,P,size(angularFilters,3));
scaleFilters = reshape(scaleFilters,P,size(scaleFilters,3));
map = bsxfun(@times,angularFilters',I_hat(1:end))*scaleFilters;

if(nargin > 3 && show)
    imagesc(abs(map));
    xlabel('Scale Index');
    ylabel('Angle Index');
end

end

