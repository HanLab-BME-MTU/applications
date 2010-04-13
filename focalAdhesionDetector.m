function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)
% [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% Make sure I is type double
I = double(I);

% Estimate background intensity (TODO: 10 is arbitrary)
ImBG = I - filterGauss2D(I,10);

% Get a first coarse segmentation
minSize = 2;
BW = logical(blobSegmentThreshold(I,minSize,1));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Get the local orientations
[R,T] = steerableFiltering(I,2,sigmaPSF); % TODO: Use M=4
R(R < 0) = 0; % Should not happen
R(BW == 0) = 0;

% Get the connected component properties
CCstats = regionprops(BW, 'Orientation', 'PixelIdxList', 'BoundingBox');
nCC = numel(CCstats);
indCC = (1:nCC)';

% Get initial segment parameters
params = getInitialSegmentParams(ImBG, R, T, sigmaPSF, CCstats);

% % params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
% params = zeros(nCC, 5);
% 
% % Segment amplitude is initialized by the mean intensity of I - BG
% params(:,3) = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});
% 
% % Segment length is initialized by the length of the CC major axis
% params(:,4) = vertcat(CCstats(:).MajorAxisLength);
% 
% % Segment orientation is initialized by the CC's main orientation
% params(:,5) = -vertcat(CCstats(:).Orientation) * pi/180;
% % Optimize initial segment orientation using local orientation map
% params(:,5) = arrayfun(@(i) lsqnonlin(@getOrientationCost, params(i,5),...
%     [],[],options, T(CCstats(i).PixelIdxList), R(CCstats(i).PixelIdxList)), ...
%     indCC);
% 
% % Segment centre is initialized using Radon transform. Each initial CC's
% % can yield multiple centers.
% 
% winSize = 2*ceil(2*sigmaPSF)+1;
% centers = arrayfun(@(i) getSegmentInitialPosition(R,CCstats(i),...
%     params(i,5),winSize),indCC, 'UniformOutput', false);
% 
% params = cell2mat(arrayfun(@(i) [centers{i}, repmat(params(i,3:5), ...
%     size(centers{i},1),1)], indCC, 'UniformOutput', false));

%Optimize all segment parameters together
%initParams = params;

% Define bounds
lb = zeros(size(params));
% min value of center coordinates
lb(:,1:2) = 1;
% min value of segment amplitude
lb(:,3) = cell2mat(arrayfun(@(i) repmat(min(ImBG(CCstats(i).PixelIdxList)), ...
    size(centers{i},1),1), indCC, 'UniformOutput', false));
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
% +1 is to prevent that lower bound equals upper bound.
ub(:,3) = cell2mat(arrayfun(@(i) repmat(max(ImBG(CCstats(i).PixelIdxList)), ...
    size(centers{i},1),1), indCC, 'UniformOutput', false)) + 1;
% max length
ub(:,4) = 2*params(:,4);
% max orientation
ub(:,5) = +inf;

% Set options for lsqnonlin
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);

fun = @(x) dLSegment2DFit(x, ImBG, sigmaPSF);
%[params, ~,residual] = lsqnonlin(fun, params, lb, ub, options);

%Im = I - reshape(residual, size(I));

function [F J] = getOrientationCost(x,varargin)

phi = varargin{1};
w = varargin{2};
w = w / sum(w);

F = w .* sin(phi-x);

if nargout > 1
    J = -w .* cos(x-phi);
end

function params = getInitialSegmentParams(ImBG,R,T,sigmaPSF,CCstats)

% Set options for lsqnonlin
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);

winSize = 2*ceil(2*sigmaPSF)+1;

nCC = numel(CCstats);
indCC = (1:nCC)';

%
% Amplitude
%

A = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});

%
% Theta
%

% First initialize the orientation using the CC's one
theta0 = -vertcat(CCstats(:).Orientation) * pi/180;
% Then optimize it using the local response nonlsq fitting
theta0 = arrayfun(@(i) lsqnonlin(@getOrientationCost, theta0(i),[],[],...
    options, T(CCstats(i).PixelIdxList), R(CCstats(i).PixelIdxList)), ...
    indCC);

%
% xC, yC and L
%

centers  = cell(nCC,1);
L = cell(nCC,1);

for i = 1:nCC
    % Floor bounding box
    bb = ceil(CCstats(i).BoundingBox);

    % Get the footprint of the CC
    BW = zeros(bb(4),bb(3));
    [y x] = ind2sub(size(R),CCstats(i).PixelIdxList);
    x = x - bb(1) + 1;
    y = y - bb(2) + 1;
    indLocal = sub2ind(size(BW), y, x);
    BW(indLocal) = 1;
    
    % Crop R using bounding box and restrinct to pixels in the CC
    Rcrop = zeros(bb(4),bb(3));
    Rcrop(indLocal) = R(CCstats(i).PixelIdxList);

    % Compute Radon transform on Rcrop
    thetaRadon = mod(pi-theta0(i)+pi/2,pi)*180/pi;
    
    % The Radon origin is floor((size(I)+1)/2).
    cRadon = floor((bb(:,3:4)+1)/2);
    
    [ccR ccXp] = radon(Rcrop,thetaRadon);

    % Find the distances ASIDE from the origin that maximize locally Radon
    % coefficients. 
    
    indMax = locmax1d(ccR,winSize);
    indMax = indMax(ccR(indMax) > 0);
    distAside = ccXp(indMax);

    p = repmat(bb(:,1:2) + cRadon - 1,numel(distAside),1) + ...
        distAside(:) * [cos(theta-pi/2), sin(theta-pi/2)];
    
    % Find the distance ALONG each line from which p will need to be
    % shifted so that it represents the center of the footprint of the CC
    % along that line.
    
    tt = tan(theta0(i));
    D = zeros(bb(4),bb(3));
    D(indLocal) = (tt * (cRadon(2) - x) + y - cRadon(1)) ./ sqrt(1 + tt^2);
    
    ccShift = radon(D,thetaRadon);
    ccShift = ccShift(indMax);
    
    centers{i} = p + ccShift;

    % Compute Radon transform on the CC's footprint to get the length along
    % each line.
    ccL = radon(BW, thetaRadon);
    L{i} = ccL(indMax);
end


