function F = anisoGaussian2Dkernel(x0, y0, amp, sigmaX, sigmaY, theta, xRange, yRange, nzIdx)
% anisotropic Gaussian 2D model defined by 6 parameters:
%    xy      : position of the segment's center
%    amp     : mean amplitude along the segment
%    sigmaX  : dispersion along the main axis
%    sigmaY  : dispersion aside the main axis
%    theta   : orientation [-pi/2, pi/2)
%
% F = anisoGaussian2D(x, y, amp, sigmaX, sigmaY, theta, xRange, yRange, nzIdx)
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

xRange = xRange - x0;
yRange = yRange - y0;

N = numel(yRange);
M = numel(xRange);

if nargin < 9 || isempty(nzIdx)
    nzIdx = 1:N*M;
end

[X Y] = meshgrid(xRange, yRange);
% X = X(nzIdx);
% Y = Y(nzIdx);

ct2 = cos(theta)^2;
st2 = sin(theta)^2;

sx2 = sigmaX^2;
sy2 = sigmaY^2;
s2t = sin(2 * theta);

a = ct2 / (2 * sx2) + st2 / (2 * sy2);
b = s2t / (4 * sx2) - s2t / (4 * sy2);
c = st2 / (2 * sx2) + ct2 / (2 * sy2);

F = amp * exp(-(a * X.^2 + 2 * b * X .* Y + c * Y.^2));
