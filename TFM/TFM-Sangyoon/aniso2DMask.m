function F = aniso2DMask(x0, y0, rx, ry, theta, xRange, yRange)
% anisotropic Gaussian 2D model defined by 6 parameters:
%    xy      : position of the segment's center
%    amp     : mean amplitude along the segment
%    rx      : radius along the main axis
%    ry      : radius aside the main axis
%    theta   : orientation [-pi/2, pi/2)
%
% F = aniso2DMask(x, y, rx, ry, theta, xRange, yRange)
%
% parameters:
% (x0,y0)            center of the model (in the image domain)
%
% sigmaX             dispersion along the main axis
%
% sigmaY             dispersion aside the main axis
%
% theta              orientation of the model 
%
% imSize             image size
%
% output:
% F                  the model defined on nzIdx pixels.
%
% Sylvain Berlemont, 2011

X = xRange - x0;
Y = yRange - y0;

ct2 = cos(theta)^2;
st2 = sin(theta)^2;

sx2 = rx^2;
sy2 = ry^2;
s2t = sin(2 * theta);

a = ct2 / (2 * sx2) + st2 / (2 * sy2);
b = s2t / (4 * sx2) - s2t / (4 * sy2);
c = st2 / (2 * sx2) + ct2 / (2 * sy2);

F = heaviside(1-(a * X.^2 + 2 * b * X .* Y + c * Y.^2));
