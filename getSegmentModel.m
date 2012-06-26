function [prm res] = getSegmentModel(x,y)
% TODO: handle basically vertical case

% Ordinary Linear Regression
% V = [x, ones(size(x))];
% [Q,R] = qr(V, 0);
% p = R \ (Q' * y);

% Deming Regression y = b0 + b1 * x, where both x and y have an error
n = numel(x);
x0 = mean(x);
y0 = mean(y);
sxx = (1/(n-1)) * sum((x - x0).^2);
sxy = (1/(n-1)) * sum((x - x0) .* (y - y0));
syy = (1/(n-1)) * sum((y - y0).^2);
delta = 1;
b1 = (syy - delta * sxx + sqrt((syy - delta * sxx)^2 + 4 * delta * sxy^2)) ...
    / (2 * sxy);

% angle
theta = atan2(b1, 1);
ct = cos(theta);
st = sin(theta);

% position
distAlong = st * (y - y0) + ct * (x - x0);
[dMin,iMin] = min(distAlong);
[dMax,iMax] = max(distAlong);
x0 = .5 * (x(iMin) + x(iMax));
y0 = .5 * (y(iMin) + y(iMax));

% length
length = dMax - dMin;

prm = [x0, y0, length, theta];

% residuals
res = st * (x - x0) - ct * (y - y0);
