function [ I ] = drawGaussianLine( angle, sigma, center, radius )
%drawGaussianLine Draw Gaussian Line

sz = [101 101];

[X,Y] = meshgrid(1:sz(1),1:sz(2));

if(nargin < 3)
    center = ceil(sz/2);
end

X = X - center(1);
Y = Y - center(2);

% ax + by + c = 0
% y = mx + d    => -mx + y - d = 0
% m = tan(theta)
% d = y - mx

m = shiftdim(tan(angle),-1);
dist = bsxfun(@times,-m,X);
dist = bsxfun(@plus,dist,Y);
% hypot(m,1) == sqrt(m.^2+1);
dist = bsxfun(@rdivide,dist,hypot(m,1));
dist = abs(dist);

if(nargin > 3)
%     keyboard;
    distToStart = hypot(X,Y);
    nearestPointToEndpoint = sqrt(max(distToStart.^2 - dist.^2,0)) - radius;
    s = nearestPointToEndpoint > 0;
    dist(s) = hypot(dist(s),nearestPointToEndpoint(s));
end

I = exp(-dist.^2/2/sigma.^2);
% I = erf((dist+0.5)/sigma/sqrt(2)) - erf((dist-0.5)/sigma/sqrt(2));

I = sum(I,3);
I = mat2gray(I);

% keyboard;


end

