function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% Make sure I is type double
I = double(I);

% Get a first coarse segmentation
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

% Define bounds
lb = zeros(size(params));
lb(:,1:2) = 1;  % min value of center coordinates
lb(:,3:4) = 0;  % min value of signal and background intensity
lb(:,5) = eps;  % min length
lb(:,6) = -inf; % min orientation value

ub = zeros(size(params));
ub(:,1) = size(I,2); % max x-coordinate
ub(:,2) = size(I,1); % max y-coordinate
ub(:,3:4) = max(I(:)); % max value of signal and background intensity
ub(:,5) = max(size(I)); % max length
ub(:,6) = +inf; % max orientation

options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
fun = @(x) dLSegment2DFit(x, I, sigmaPSF);
[params, ~,residual,exitflag,output,lambda,Jon] = lsqnonlin(fun, initParams, lb, ub, options);

% options = optimset('Jacobian', 'off', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
%     'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
% [params, ~,residual,exitflag,output,lambda,Joff] = lsqnonlin(fun, initParams, [], [], options);

% for j = 1:size(params,1)
%     figure,
%     for i = 1:6
%         subplot(6, 2, 2*i-1); imshow(reshape(full(Jon(:, j + (i-1) * size(params,1))), size(I)),[]);
%         subplot(6, 2, 2*i); imshow(reshape(full(Joff(:, j + (i-1) * size(params,1))), size(I)),[]);
%     end
% end

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
