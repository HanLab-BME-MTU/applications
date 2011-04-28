function [prm res] = getSegmentModel(x,y)
% TODO: handle basically vertical case

V = [x, ones(size(x))];
[Q,R] = qr(V, 0);
p = R \ (Q' * y);

% angle
theta = atan2(p(1), 1);
ct = cos(theta);
st = sin(theta);

% position
x0 = mean(x);
y0 = mean(y);
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
