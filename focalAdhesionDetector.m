function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% make sure I is type double
I = double(I);

% Get a first crude segmentation
BW = logical(blobSegmentThreshold(I,1,0));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Remove noise by filtering image with a Gaussian whose sigma = sigmaPSF
% NOTE: this is not optimal since this step has already been done in the
% blobSegmentThreshold function.
If = Gauss2D(I,sigmaPSF,1);

% Estimate background by filtering image with a Gaussian whose sigma = 10
% NOTE: this is not optimal since this step has already been done in the
% blobSegmentThreshold function.
BG = Gauss2D(I,10,1);

% Update the image
Iu = If - BG;

% Get initial rod parameter from the the connected component properties
CCstats = regionprops(BW, Iu, 'Centroid','Orientation','MajorAxisLength','PixelValues');

% params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
params = zeros(numel(CCstats), 5);
params(:,1:2) = vertcat(CCstats(:).Centroid);
params(:,3) = cellfun(@(px) mean(px), {CCstats(:).PixelValues});
params(:,4) = vertcat(CCstats(:).MajorAxisLength);
params(:,5) = -vertcat(CCstats(:).Orientation) * pi/180;

% Fit
initParams = params;
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
fun = @(x) dlSegment2DFit(x, Iu, sigmaPSF);
[params, ~, residual] = lsqnonlin(fun, initParams, [], [], options);

Im = Iu - residual;

% DEBUG
imshow(Iu, []); hold on;
quiver(initParams(:,1), initParams(:,2), (initParams(:,4) / 2) .* cos(initParams(:,5)), ...
    (initParams(:,4) / 2) .* sin(initParams(:,5)), 0, 'g');
quiver(initParams(:,1), initParams(:,2), (initParams(:,4) / 2) .* cos(initParams(:,5) + pi), ...
    (initParams(:,4) / 2) .* sin(initParams(:,5) + pi), 0, 'g');
plot(initParams(:,1), initParams(:,2), 'r.');

quiver(params(:,1), params(:,2), (params(:,4) / 2) .* cos(params(:,5)), ...
    (params(:,4) / 2) .* sin(params(:,5)), 0, 'g');
quiver(params(:,1), params(:,2), (params(:,4) / 2) .* cos(params(:,5) + pi), ...
    (params(:,4) / 2) .* sin(params(:,5) + pi), 0, 'g');
plot(params(:,1), params(:,2), 'r.'); hold off;

