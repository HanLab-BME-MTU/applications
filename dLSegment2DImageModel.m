function Im = dLSegment2DImageModel(params, sigmaPSF, imSize)
% Compute an image model which is the sum of 2D diffraction-limited segment
% models defined by params (see dLSegment2D.m for more details).
%
% Im = dLSegment2DImageModel(params, sigmaPSF, imSize)
%
% parameters:
% params          nx5 matrix where n is the number of segments and their
%                 parameters, i.e. xC, yC, A, l, t are stored column-wise.
%
% sigmaPSF        half width of the gaussian PSF model
%
% imSize          image size

[n,p] = size(params);

if p ~= 5
    error('Invalid number of segment parameters.');
end

% Generate the image segments
Im = zeros(imSize);

for i = 1:n
    xC = params(i,1);
    yC = params(i,2);
    A = params(i,3);
    l = params(i,4);
    t = params(i,5);
    
    [xRange yRange] = dLSegment2DSupport(xC, yC, sigmaPSF, l, t);

    xRange = max(xRange(1),1):min(xRange(end),imSize(2));
    yRange = max(yRange(1),1):min(yRange(end),imSize(1));
    
    S = dLSegment2D(xRange, yRange, xC, yC, A, sigmaPSF, l, t);

    Im(yRange,xRange) = Im(yRange,xRange) + S;
end