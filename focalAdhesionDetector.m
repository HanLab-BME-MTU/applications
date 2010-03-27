function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% Make sure I is type double
I = double(I);

% Get a first coarse segmentation
BW = logical(blobSegmentThreshold(I,1,0));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Remove background from image
Ib = filterGauss2D(I,10);
I = I - Ib;

% Get initial rod parameter from the the connected component properties
CCstats = regionprops(BW, I, 'Centroid','Orientation','MajorAxisLength', ...
    'PixelValues');

% params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
params = zeros(numel(CCstats), 5);
% Segment centre is initialized by the connected component (CC) centre
params(:,1:2) = vertcat(CCstats(:).Centroid);
% Segment amplitude is initialized by the mean of the CC
params(:,3) = cellfun(@(px) mean(px), {CCstats(:).PixelValues});
% Segment length is initialized by the length of the CC major axis
params(:,4) = vertcat(CCstats(:).MajorAxisLength);
% Segment orientation is initialized by the CC's main orientation
params(:,5) = -vertcat(CCstats(:).Orientation) * pi/180;

% Fit
initParams = params;

% Define bounds
lb = zeros(size(params));
% min value of center coordinates
lb(:,1:2) = 1;
% min value of segment amplitude
lb(:,3) = cellfun(@(px) min(px), {CCstats(:).PixelValues}) - 3 * sigmaPSF;
% min length
lb(:,4) = .25*params(:,4);
% min orientation value
lb(:,5) = -inf;

ub = zeros(size(params));
% max x-coordinate
ub(:,1) = size(I,2);
% max y-coordinate
ub(:,2) = size(I,1);
% max value of signal intensity
ub(:,3) = cellfun(@(px) max(px), {CCstats(:).PixelValues}) + 3 * sigmaPSF;
% max length
ub(:,4) = 2*params(:,4);
% max orientation
ub(:,5) = +inf;

options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
fun = @(x) dLSegment2DFit(x, I, sigmaPSF);
[params, ~,residual,exitflag] = lsqnonlin(fun, initParams, lb, ub, options);

Im = I - reshape(residual, size(I));

% DEBUG
imshow(I, []); hold on;
xC = initParams(:,1);
yC = initParams(:,2);
l = initParams(:,4);
t = initParams(:,5);

quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'g');
quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'g');
plot(xC, yC, 'g.');

xC = params(:,1);
yC = params(:,2);
l = params(:,4);
t = params(:,5);

quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'r');
quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'r');
plot(xC, yC, 'r.');
