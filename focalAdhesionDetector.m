function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)
% [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% Options used by lsqnonlin.
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);

% Make sure I is type double
I = double(I);

% Estimate background intensity (TODO: 10 is arbitrary)
ImBG = I - filterGauss2D(I,10);

% Get a first coarse segmentation
BW = logical(blobSegmentThreshold(I,1,0));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Get the local orientations
[R,T] = steerableFiltering(I,2,sigmaPSF); % TODO: Use M=4
R(R < 0) = 0; % Should not happen
R(BW == 0) = 0;

% Get the connected component properties
CCstats = regionprops(BW, 'Orientation','MajorAxisLength', 'PixelIdxList', ...
    'BoundingBox');
nCC = numel(CCstats);
indCC = (1:nCC)';

% params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
params = zeros(nCC, 5);

% Segment amplitude is initialized by the mean intensity of I - BG
params(:,3) = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});

% Segment length is initialized by the length of the CC major axis
params(:,4) = vertcat(CCstats(:).MajorAxisLength);

% Segment orientation is initialized by the CC's main orientation
params(:,5) = -vertcat(CCstats(:).Orientation) * pi/180;
% Optimize initial segment orientation using local orientation map
params(:,5) = arrayfun(@(i) lsqnonlin(@getOrientationCost, params(i,5),...
    [],[],options, T(CCstats(i).PixelIdxList), R(CCstats(i).PixelIdxList)), ...
    indCC);

% Segment centre is initialized using Radon transform
params(:,1:2) = cell2mat(arrayfun(@(i) getSegmentInitialPosition(R,CCstats(i),...
    params(i,5)),indCC, 'UniformOutput', false));

% Optimize all segment parameters together
%initParams = params;

% Define bounds
lb = zeros(size(params));
% min value of center coordinates
lb(:,1:2) = 1;
% min value of segment amplitude
lb(:,3) = cellfun(@(ind) min(ImBG(ind)), {CCstats(:).PixelIdxList}) - 3 * sigmaPSF;
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
ub(:,3) = cellfun(@(ind) max(ImBG(ind)), {CCstats(:).PixelIdxList}) + 3 * sigmaPSF;
% max length
ub(:,4) = 2*params(:,4);
% max orientation
ub(:,5) = +inf;

fun = @(x) dLSegment2DFit(x, ImBG, sigmaPSF);
[params, ~,residual] = lsqnonlin(fun, params, lb, ub, options);

Im = I - reshape(residual, size(I));

% DEBUG
% imshow(I, []); hold on;
% 
% xC = initParams(:,1);
% yC = initParams(:,2);
% l = initParams(:,4);
% t = initParams(:,5);
% 
% quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'g');
% quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'g');
% plot(xC, yC, 'g.');
% 
% xC = params(:,1);
% yC = params(:,2);
% l = params(:,4);
% t = params(:,5);
% 
% quiver(xC, yC, (l / 2) .* cos(t), (l / 2) .* sin(t), 0, 'r');
% quiver(xC, yC, (l / 2) .* cos(t + pi), (l / 2) .* sin(t + pi), 0, 'r');
% plot(xC, yC, 'r.');

function [F J] = getOrientationCost(x,varargin)

phi = varargin{1};
w = varargin{2};
w = w / sum(w);

F = w .* sin(phi-x);

if nargout > 1
    J = -w .* cos(x-phi);
end

function pos = getSegmentInitialPosition(R,propCC,theta)

% Floor bounding box
bb = ceil(propCC.BoundingBox);

% Crop R using bounding box and restrinct to pixels in the cc
Rcrop = zeros(bb(4),bb(3));
[y x] = ind2sub(size(R),propCC.PixelIdxList);
x = x - bb(1) + 1;
y = y - bb(2) + 1;
indLocal = sub2ind(size(Rcrop), y, x);
Rcrop(indLocal) = R(propCC.PixelIdxList);

% Compute Radon transform
[ccR ccXp] = radon(Rcrop,(-theta + pi/2)*180/pi);

% Find the distance from the origin that maximize Radon coefficients. The
% origin is floor((size(I)+1)/2).
[~, maxInd] = max(ccR);
d = ccXp(maxInd);

pos = bb(:,1:2) + floor((bb(:,3:4)+1)/2) - 1 + ...
    bsxfun(@times,[cos(theta-pi/2), sin(theta-pi/2)], d);

