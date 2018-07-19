function [ I ] = drawPartialGaussianLine( angle, sigma, startRadius, stopRadius)
%drawPartialGaussianLine Summary of this function goes here
%   Detailed explanation goes here

radius = (stopRadius-startRadius)/2;
offset = (startRadius+stopRadius)/2;
sz = [101 101];
angle = shiftdim(angle(:),-2);
center = bsxfun(@plus,[cos(angle) sin(angle)]*offset,sz/2);

I = zeros([sz length(angle)]);

for i=1:size(center,3)
    I(:,:,i) = intersections.drawGaussianLine(angle(i),sigma,center(:,:,i),radius);
end

I = sum(I,3);


end

