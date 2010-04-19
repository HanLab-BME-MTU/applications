function [params, ImBG] = getInitialSegmentParamsNEW(I,sigmaPSF,minSize)

% Make sure image's class is double
if ~isa(I,'double')
    I = double(I);
end

% Get a first coarse segmentation
BW = logical(blobSegmentThreshold(I,minSize,1));

% Get the connected component properties
CCstats = regionprops(BW, 'Orientation', 'PixelIdxList', 'BoundingBox');

% Estimate amplitude

% Estimate background intensity (TODO: 10 is arbitrary)
ImBG = I - filterGauss2D(I,10);

amp0 = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});

nCC = numel(CCstats);
indCC = (1:nCC)';

% Define filter banks for 2nd order ridge detector
hside = ceil(6 * sigmaPSF);
x = -hside:hside;
gKernel = exp(-x.^2 / (2 * sigmaPSF^2));
g0 = (1 / (sqrt(2 * pi) * sigmaPSF)) * gKernel;
g1 = (- x ./ (sqrt(2 * pi) * sigmaPSF^3)) .* gKernel;
g2 = (x.^2 - sigmaPSF^2) ./ (sqrt(2 * pi) * sigmaPSF^5) .* gKernel;
% M=2 specific constant
a20 = sigmaPSF / (2 * sqrt(3 * pi));
a22 = - sqrt(3 / (4 * pi)) * sigmaPSF;

% Filter image with basis filters
Ixx = conv2(g0,g2,I,'same');
Ixy = conv2(g1,g1,I,'same');
Iyy = conv2(g2,g0,I,'same');

Ixx = Ixx .* BW;
Ixy = Ixy .* BW;
Iyy = Iyy .* BW;

% Radon angle range: this angle correspond to the orientation of the
% projector, i.e. the line where the image will be projected on. It
% corresponds then to the perpendicular of the projection orientation.
thetaRadon = 0:179;

% Equivalent in radian in [0...-pi/2] U [pi/2...0] (clockwise)
theta = [0:-1:-90, 89:-1:1] * pi / 180;

% t90 is theta +/- pi/2 (still belongs to -pi/2, pi/2)
ind = theta >= -pi/2 & theta < 0;
theta90 = zeros(size(theta));
theta90(ind) = theta(ind) + pi/2;
theta90(~ind) = theta(~ind) - pi/2;

% Pre-compute trigonometric functions
ct = cos(theta);
st = sin(theta);
tt = tan(theta);

ct90 = cos(theta90);
st90 = sin(theta90);
ct90_2 = ct90.^2;
st90_2 = st90.^2;

% Initialize output arrays
theta0 = cell(nCC,1);
length0 = cell(nCC,1);
centers0 = cell(nCC,1);

% Minimal distance between local maxima
winSize = 2*ceil(2*sigmaPSF)+1;

for iCC = 1:nCC    
    bb = ceil(CCstats(iCC).BoundingBox);
    
    % The Radon origin is floor((size(I)+1)/2).
    cRadon = floor((bb(:,3:4)+1)/2);

    % Get the footprint of the CC
    BWcrop = false(bb(4),bb(3));
    indGlobal = CCstats(iCC).PixelIdxList;
    [y x] = ind2sub(size(I),indGlobal);
    x = x - bb(1) + 1;
    y = y - bb(2) + 1;
    indLocal = sub2ind(size(BWcrop), y, x);
    BWcrop(indLocal) = true;
    
    % Filter Icrop using filter banks of M=2
    IxxCrop = zeros(bb(4),bb(3));
    IxxCrop(indLocal) = Ixx(indGlobal);
    IxyCrop = zeros(bb(4),bb(3));
    IxyCrop(indLocal) = Ixy(indGlobal);
    IyyCrop = zeros(bb(4),bb(3));
    IyyCrop(indLocal) = Iyy(indGlobal);
    
    % Compute Radon transform on filtered images
    [Rxx,xp] = radon(IxxCrop,thetaRadon);
    Rxy = radon(IxyCrop,thetaRadon);
    Ryy = radon(IyyCrop,thetaRadon);
    
    % Combine RTs
    R = repmat(ct90_2,numel(xp),1) .* (Rxx * a20 + Ryy * a22) + ...
        repmat(st90_2,numel(xp),1) .* (Ryy * a20 + Rxx * a22) + ...
        2 * repmat(ct90,numel(xp),1) .* repmat(st90,numel(xp),1) .* ...
        Rxy * (a20 - a22);
    
    % Compute Radon transform on BWcrop
    L = radon(BWcrop,thetaRadon);
    
    % Compute the mean line integral over the filtered image
    R = R ./ sqrt(L);

    % Set a threshold based on the variance of the data
    R(R <= eps | L < minSize) = NaN;
    th = 3 * nanstd(R(:));

    % Get the local maxima along theta
    RlocMax = zeros(size(R));
    for iTheta = 1:numel(theta)
        indMax = locmax1d(R(:,iTheta), winSize);
        indMax = indMax(R(indMax,iTheta) > eps & L(indMax,iTheta) >= minSize);
        RlocMax(indMax,iTheta) = R(indMax,iTheta);
    end
    
    % Get the peak in the Radon transform
    [~,iThetaMax] = max(sum(RlocMax));
    Rmax = max(RlocMax(:, iThetaMax));
    
    % If the peak is significant, get every local maxima along iThetaMax
    if Rmax >= th     
        indMax = find(RlocMax(:,iThetaMax));
        indMax = indMax(L(indMax,iThetaMax) >= minSize);
        
        % Save the initial orientation
        theta0{iCC} = theta90(iThetaMax);
        
        % Save the initial set of lengths
        length0{iCC} = L(indMax, iThetaMax);
        
        distAside = xp(indMax);
        
        D = zeros(bb(4),bb(3));
        D(indLocal) = (tt(iThetaMax) * (cRadon(1) - x) + y - ...
            cRadon(2)) ./ sqrt(1 + tt(iThetaMax)^2);
        
        distAlong = radon(D,thetaRadon(iThetaMax));
        distAlong = distAlong(indMax) ./ length0{iCC};
        
        pts = [distAside distAlong];
        
        % Rotate pts
        rotT = [ct(iThetaMax) st(iThetaMax); ...
            -st(iThetaMax) ct(iThetaMax)];
    
        % Centers of FAs = center of the CC + Rot(pts)
        centers0{iCC} = repmat(bb(:,1:2) + cRadon - 1,numel(distAside),1) + ...
            pts * rotT;
    end
end

% Concatenate all parameters
nEmpty = cellfun(@(x) ~isempty(x), centers0);

params = cell2mat(arrayfun(@(i) [...
    centers0{i},...                               % Xc, Yc
    repmat(amp0(i), size(centers0{i},1),1),...    % amp
    length0{i},...                                % length
    repmat(theta0{i}, size(centers0{i},1),1)],... % theta
    indCC(nEmpty), 'UniformOutput', false));
