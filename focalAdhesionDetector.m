function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% make sure I is type double
I = double(I);

% Get a first crude segmentation
BW = logical(blobSegmentThreshold(I,1,0));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Get initial rod parameter from the the connected component properties
CCstats = regionprops(BW, I, 'Centroid','Orientation','MajorAxisLength','PixelValues');

% params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
params = zeros(numel(CCstats), 6);
% Segment centre is initialized by the connected component (CC) centre
params(:,1:2) = vertcat(CCstats(:).Centroid);
% Segment amplitude is initialized by the mean of the CC
params(:,3) = cellfun(@(px) mean(px), {CCstats(:).PixelValues});
% Segment background intensity is initialized by the min value of the CC
params(:,4) = cellfun(@(px) min(px), {CCstats(:).PixelValues});
% Segment length is initialized by the length of the CC major axis
params(:,5) = vertcat(CCstats(:).MajorAxisLength);
% Segment orientation is initialized by the CC's main orientation
params(:,6) = -vertcat(CCstats(:).Orientation) * pi/180;

% Fit
initParams = params;
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
fun = @(x) dLSegment2DFit(x, I, sigmaPSF);
[params, ~, residual] = lsqnonlin(fun, initParams, [], [], options);

Im = I - reshape(residual, size(I));

% DEBUG
imshow(I, []); hold on;
xC = initParams(:,1);
yC = initParams(:,2);
l = initParams(:,5);
t = initParams(:,6);

quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'g');
quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'g');
plot(xC, yC, 'g.');

xC = params(:,1);
yC = params(:,2);
l = params(:,5);
t = params(:,6);

quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'r');
quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'r');
plot(xC, yC, 'r.');
