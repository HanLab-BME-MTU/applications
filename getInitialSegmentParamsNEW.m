function [params, ImBG] = getInitialSegmentParamsNEW(I,sigmaPSF,minSize)

% Make sure image's class is double
if ~isa(I,'double')
    I = double(I);
end

% Get a first coarse segmentation
BW = logical(blobSegmentThreshold(I,minSize,0));

% Get the connected component properties
CCstats = regionprops(BW, 'Orientation', 'PixelIdxList', 'BoundingBox');

% Estimate amplitude

% Estimate background intensity (TODO: 10 is arbitrary)
ImBG = I - filterGauss2D(I,10);

amp0 = cellfun(@(idx) mean(ImBG(idx)), {CCstats(:).PixelIdxList});

nCC = numel(CCstats);

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

% Radon angle range
thetaRadon = 0:179;

% Equivalent in radian in [0...-pi/2] U [pi/2...0] (clockwise)
theta = [0:-1:-90, 89:-1:1] * pi / 180;

% Pre-compute trigonometric functions
ct = cos(theta);
ct2 = ct.^2;
st = sin(theta);
st2 = st.^2;

% t90 is theta +/- pi/2 (still belongs to -pi/2, pi/2)
ind = theta >= -pi/2 & theta < 0;
t90 = zeros(size(theta));
t90(ind) = theta(ind) + pi/2;
t90(~ind) = theta(~ind) - pi/2;

% Pre-compute trigonometric functions
ct_90 = cos(t90);
st_90 = sin(t90);
tt_90 = tan(t90);

% Initialize output arrays
theta0 = cell(nCC,1);
length0 = cell(nCC,1);
centers0 = cell(nCC,1);

% Temporary arrays
accIndMax = cell(1, numel(theta));

% Minimal distance between local maxima
winSize = 2*ceil(2*sigmaPSF)+1;

for iCC = 1:nCC
    % Floor bounding box
    bb = ceil(CCstats(iCC).BoundingBox);
    
    % The Radon origin is floor((size(I)+1)/2).
    cRadon = floor((bb(:,3:4)+1)/2);

    % Get the footprint of the CC
    BWcrop = zeros(bb(4),bb(3));
    [y x] = ind2sub(size(I),CCstats(iCC).PixelIdxList);
    x = x - bb(1) + 1;
    y = y - bb(2) + 1;
    indLocal = sub2ind(size(BWcrop), y, x);
    BWcrop(indLocal) = 1;
    
    % Get the image crop
    Icrop = I(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1);
    
    % Filter Icrop using filter banks of M=2
    Ixx = conv2(g0,g2,Icrop,'same');
    Ixy = conv2(g1,g1,Icrop,'same');
    Iyy = conv2(g2,g0,Icrop,'same');
    
    Ixx = Ixx .* BWcrop;
    Ixy = Ixy .* BWcrop;
    Iyy = Iyy .* BWcrop;
    
    % Compute Radon transform on filtered images
    [Rxx,xp] = radon(Ixx,thetaRadon);
    Rxy = radon(Ixy,thetaRadon);
    Ryy = radon(Iyy,thetaRadon);
    
    % Combine RTs
    R = arrayfun(@(iTheta) ...
        ct2(iTheta) * (Rxx(:,iTheta) * a20 + Ryy(:,iTheta) * a22) + ...
        st2(iTheta) * (Ryy(:,iTheta) * a20 + Rxx(:,iTheta) * a22) + ...
        2 * ct(iTheta) * st(iTheta) * Rxy(:,iTheta) * (a20 - a22), ...
        1:numel(theta), 'UniformOutput', false);
    
    R = cell2mat(R);
    
    % Compute Radon transform on BWcrop
    L = radon(BWcrop,thetaRadon);
    
    % Compute the mean line integral over the filtered image
    R = R ./ L;

    % Set a threshold using the variance of the data
    R(R <= 0) = NaN;
    th = 3 * nanstd(R(:));
    
    % For each angle (column), compute the indices (rows) whose value in R
    % is greater than th and the length is greater that minSize. Store
    % these indices in a cell array accIndMax.
    
    for iTheta = 1:numel(theta)        
        indMax = locmax1d(R(:,iTheta), winSize);
        indMax = indMax(R(indMax,iTheta) >= th & L(indMax,iTheta) >= minSize);
        
        accIndMax{iTheta} = indMax;
    end
    
    if any([accIndMax{:}])
  
        % TODO !!!!!
        
        %maxSignif = arrayfun(@(iTheta) max(R(accIndMax{iTheta}, iTheta)),...
        %1:numel(theta), 'UniformOutput', false);
    

        
        [~, iThetaMax] = max(nSignifIndMax);

        % Save the initial orientation
        theta0{iCC} = theta(iThetaMax);
        
        % Save the initial set of lengths
        length0{iCC} = L(accIndMax{iThetaMax}, iThetaMax);
        
        distAside = xp(accIndMax{iThetaMax}, iThetaMax);
        
        D = zeros(bb(4),bb(3));
        D(indLocal) = (tt_90(iThetaMax) * (cRadon(1) - x) + y - ...
            cRadon(2)) ./ sqrt(1 + tt_90(iThetaMax)^2);
        
        distAlong = radon(D,thetaRadon(iThetaMax));
        distAlong = distAlong(accIndMax{iThetaMax}) ./ length0{iCC};
        
        pts = [distAside distAlong];
        
        % Rotate pts
        rotT90 = [ct_90(iThetaMax) st_90(iThetaMax); ...
            -st_90(iThetaMax) ct_90(iThetaMax)];
    
        % Centers of FAs = center of the CC + Rot90(pts)
        centers0{iCC} = repmat(bb(:,1:2) + cRadon - 1,numel(distAside),1) + ...
            pts * rotT90;
    end
end

% Concatenate all parameters
nEmpty = cellfun(@(x) ~isempty(x), centers);

params = cell2mat(arrayfun(@(i) [...
    centers0{i},...                               % Xc, Yc
    repmat(amp0(i), size(centers0{i},1),1),...    % amp
    length0{i},...                                % length
    repmat(theta0(i), size(centers0{i},1),1)],... % theta
    indCC(nEmpty), 'UniformOutput', false));
