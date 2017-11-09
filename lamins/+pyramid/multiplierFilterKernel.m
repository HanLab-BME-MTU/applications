function [ multF , multA , maxima] = multiplierFilterKernel( imSize, coords,indices )
%multiplierFilter Obtain a filter set of multiplier filters
%
% size of the image as a scalar or two-element vector
% coords as from orientationSpace.getFrequencySpaceCoordinates
%    .f - frequency magnitude, from 0 to sqrt(2)/2
%    .theta - frequency angle (-pi,pi]
% indices index of multiplier

if(nargin < 2)
    coords = orientationSpace.getFrequencySpaceCoordinates(imSize);
end
if(nargin < 3)
%     indices = -9:27;
    indices = 0:15;
end

omega = 2;
N = 9;

Mi = pyramid.nopi.multiplierFxnIsolate(omega,N);

multF = Mi(coords.f,indices);

if(nargout > 1)
    maxima = pyramid.nopi.multiplierFxnMaxima(indices,N);
    multA = squeeze(Mi(maxima,indices));
end

end

