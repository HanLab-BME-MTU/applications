function [ response ] = multiplierFilter( I, coords,indices )
%multiplierFilter Obtain a filter set of multiplier filters

imSize = size(I);


if(nargin < 2)
    coords = orientationSpace.getFrequencySpaceCoordinates(imSize);
end
if(nargin < 3)
    indices = -9:27;
end


[mF, mA] = pyramid.multiplierFilterKernel(imSize, coords, indices);

response = ifft2(bsxfun(@times,fft2(I),mF));

end

