function [ map ] = angleScaleMap( I, varargin )
%ANGLESCALEMAP Produces a map that represents the overall tiling of
%frequency space by angular and scale filters
%
% INPUT
% I - the original image data, untranformed by fft2, X by Y
% angularFilters - frequency filter over angular space, X by Y by A
% scaleFilters - frequency filter over scale space (radial frequency), X by Y
%                by S
%
% OUTPUT
% map - tile map of frequency space in complex form, A by S

I_hat = fft2(I);
map = pyramid.angleScaleMapFreq(I_hat, varargin{:});

end

