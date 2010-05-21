function [params, ImBG] = getInitialSegmentParams(ima,mask,sigmaPSF,minSize)

% Make sure image's class is double
if ~isa(ima,'double')
    ima = double(ima);
end

% Get a first coarse segmentation
BW = logical(blobSegmentThreshold(ima,minSize,0,mask));

% Get the local orientations
[R,T] = steerableFiltering(ima,2,sigmaPSF); % TODO: Use M=4
R(R < 0) = 0; % Should not happen
R(BW == 0) = 0;

% Get the connected component properties
CCstats = regionprops(BW, 'Orientation', 'PixelIdxList', 'BoundingBox');

winSize = 2*ceil(2*sigmaPSF)+1;

nCC = numel(CCstats);
indCC = (1:nCC)';

%
% Amplitude
%

% Estimate background intensity (TODO: 10 is arbitrary)
ImBG = ima - filterGauss2D(ima,10);

A = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});

%
% Theta
%

% Set options for lsqnonlin
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4,...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);

% First initialize the orientation using the CC's one
theta0 = -vertcat(CCstats(:).Orientation) * pi/180;
% Then optimize it using the local response nonlsq fitting
theta0 = arrayfun(@(i) lsqnonlin(@getOrientationCost, theta0(i),[],[],...
    options, T(CCstats(i).PixelIdxList), R(CCstats(i).PixelIdxList)), ...
    indCC);

% lsqnonlin could yield theta value outside [-pi/2, pi/2] due to the fact
% we don't bound the optimization. Make sure theta0 is in [-pi/2,pi/2]
theta0 = rem(theta0,pi);
ind = theta0 < -pi/2 & theta0 >= -pi;
theta0(ind) = theta0(ind) + pi;
ind = theta0 > pi/2 & theta0 <= pi;
theta0(ind) = theta0(ind) - pi;

% t90 is theta +/- pi/2 (still belongs to -pi/2, pi/2)
ind = theta0 >= -pi/2 & theta0 < 0;
t90 = zeros(size(theta0));
t90(ind) = theta0(ind) + pi/2;
t90(~ind) = theta0(~ind) - pi/2;

ct_90 = cos(t90);
st_90 = sin(t90);
tt_90 = tan(t90);

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
    
    % Crop R using bounding box and restrinct to CC's footprint
    Rcrop = zeros(bb(4),bb(3));
    Rcrop(indLocal) = R(CCstats(i).PixelIdxList);

    % strategy: so far, we have the angle (theta0) along which every FA is
    % supposed to be aligned in the CC. We have to find the number of FA in
    % the CC and their position along a line perpendicular to theta0 and
    % passing through the center of the CC. We use the Radon transform to
    % achieve that goal. We use the center defined by Radon as the center
    % of the CC (i.e center = floor((bb(:,3:4)+1)/2). See Radon's help).
    % Along the line perpendicular to theta0 (i.e. along t90), we want to
    % find the position (called distAside) that are local maxima of the
    % mean integral of R along theta0. This correspond to the find every
    % local maximum of the 1-dimentional signal radon(Rcrop,theta0) ./
    % radon(BW,theta0).
    
    % The Radon origin is floor((size(BW)+1)/2).
    cRadon = floor((bb(:,3:4)+1)/2);
    
    % Radon angle is defined in degree in a counterclockwise axis
    thetaRadon = -t90(i)*180/pi;
    
    % Compute Radon transform on Rcrop
    [ccR ccXp] = radon(Rcrop,thetaRadon);
    
    % Compute Radon transform on the CC's footprint (ccXp will be the
    % same as above).
    ccL = radon(BW, thetaRadon);

    % Compute the mean integral of Rcrop along the line oriented along t90.
    ccMeanR = ccR ./ (ccL + 1e-10);
    
    % Find local maxima, under the constraint of a minimal length.
    indMax = locmax1d(ccMeanR, winSize);
    indMax = indMax(ccMeanR(indMax) > 0 & ccL(indMax) >= minSize);
    
    distAside = ccXp(indMax);
    
    % Store the length of each segment
    L{i} = ccL(indMax);
    
    % Now that we have the distance aside the main CC's orientation, we are
    % able to defines lines that are oriented along theta0 and passing
    % through the center of the CC shifted by distAside. We hope this set
    % of lines overlap real FA's in bundle. Still, each point traversed by
    % a line may not represent the real center of a FA. we need to find the
    % center of each FA along each line. Again, we use Radon to compute
    % these positions. distAlong represents the distance along each line
    % (those along FAs) away from the perpendicular line (perpendicular to
    % FAs) passing through the center of the CC.
    
    % D is the signed distance transform from the perpendicular axis.
    D = zeros(bb(4),bb(3));
    D(indLocal) = (tt_90(i) * (cRadon(1) - x) + y - cRadon(2)) ./ ...
        sqrt(1 + tt_90(i)^2);
    
    % The integration of that distance transform along lines and restricted
    % to the CC's footprint will gives the center of each FA.
    
    distAlong = radon(D,thetaRadon);
    distAlong = distAlong(indMax) ./ L{i};
    
    % distAside and distAlong can be empty. We want to generate 0x2 matrix
    % so that pts remains compliant with the next set of operations.
    if isempty(distAside)
        pts = zeros(0,2);
    else
        pts = [distAside distAlong];
    end
    
    % Rotate pts
    rotT90 = [ct_90(i) st_90(i); -st_90(i) ct_90(i)];
    
    % Centers of FAs = center of the CC + Rot90(pts)
    centers{i} = repmat(bb(:,1:2) + cRadon - 1,numel(distAside),1) + ...
        pts * rotT90;
end

% Concatenate all parameters
nEmpty = cellfun(@(x) ~isempty(x), centers);

params = cell2mat(arrayfun(@(i) [...
    centers{i},...                               % Xc, Yc
    repmat(A(i), size(centers{i},1),1),...       % A
    L{i},...                                     % L
    repmat(theta0(i), size(centers{i},1),1)],... % theta
    indCC(nEmpty), 'UniformOutput', false));

function [F J] = getOrientationCost(x,varargin)

phi = varargin{1};
w = varargin{2};
w = w / sum(w);

F = w .* sin(phi-x);

if nargout > 1
    J = -w .* cos(x-phi);
end
