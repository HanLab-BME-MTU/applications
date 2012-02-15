function computeOrientation(obj,filterLength,filterAngularSampling,nBins,minBinResponse)
% function computeOrientation(obj,filterLength,filterAngularSampling)
% SYNOPSIS:
% Computes the orientation of each point by filtering the data with a
% rotated template centered at the data points.
%
% REQUIRED INPUTS:
% - filterLength
% The length of the orientation filter. The size of the filter support
%
% - filterAngularSampling
% The angular sampling step for phi and theta in degrees.
%
% - nBins
% Number of bins along the filter. The responses of these bins are
% multiplied with each other.
%
% - minBinResponse
% Minimum response of a bin in case it should be empty.
%
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES:
% - obj.data.points
% - obj.data.nPoints
% - obj.data.error
%
% MODIFIED PROPERTIES:
% - obj.data.orientation
% - obj.data.magnitude
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

% The mean deviation from the filament center
meanLocPrec = mean(obj.data.error,1);

% Init data
obj_data_orientation = zeros(obj.data.nPoints,3);
obj_data_magnitude = zeros(obj.data.nPoints,1);
obj_data_nPoints = obj.data.nPoints;
obj_data_points = obj.data.points;

% Find the neighbors of the points
[idxs,~] = KDTreeBallQuery(obj.data.points,obj.data.points,filterLength/2);

% Sampling angles
phiAngles = 0:filterAngularSampling:180-filterAngularSampling;
thetaAngles = 0:filterAngularSampling:180-filterAngularSampling;
[phi,theta] = meshgrid(phiAngles,thetaAngles);

% Sampling orientations
directionVector = zeros(1,3,size(theta,1),size(theta,2));
directionVector(1,1,:,:) = sind(theta).*cosd(phi)/meanLocPrec(1);
directionVector(1,2,:,:) = sind(theta).*sind(phi)/meanLocPrec(2);
directionVector(1,3,:,:) = cosd(theta)/meanLocPrec(3);

parfor i=1:obj_data_nPoints
    % The point that is currently processed
    centerPoint = obj_data_points(i,:)./meanLocPrec;
    centerPoints = repmat(centerPoint,[size(idxs{i},1),1,size(phiAngles,2),size(thetaAngles,2)]);
    
    % All orientations that are evaluated for this point
    directionVectorGrid = repmat(directionVector,[size(idxs{i},1),1,1,1]);
    
    % The neighbor points including the point itself
    points = obj_data_points(idxs{i},:);
    points = points./repmat(meanLocPrec,size(points,1),1);
    points = repmat(points,[1,1,size(phiAngles,2),size(thetaAngles,2)]);
    
    % Project the points onto the orientation vector
    t = sum(directionVectorGrid.*(points-centerPoints),2)./sqrt(sum(directionVectorGrid.^2,2));
    
    % Normalize the position values of the points to [0,1]
    len = sqrt(sum((directionVectorGrid*filterLength/2).^2,2));
    tNorm = t./len;
    tNormRep = repmat(tNorm,[1,3,1,1]);

    % Compute the exponential value of this distance
    orthogonalVectors = points-(centerPoints+tNormRep.*directionVectorGrid*filterLength/2);
    distanceSq = sum(orthogonalVectors.^2,2); % Squared distances
    
    % CHOICE: Compute the contribution to the response. The distances have a
    % variance of 2. Two times the standard deviation (95%) (four times the 
    % variance) is used here. => 2*(2*sqrt(2))^2 = 16
    contribution = exp(-distanceSq/16);
    
    % If the localization precision is anisotropic tNorm can be bigger than
    % 1. Set the contribution of these points to 0.
    contribution(abs(tNorm) > 1) = 0;
    
    % Find the first and last projection that have a non zero contribution
    tNormMin = min(min(tNorm.*double(contribution>0),[],1),0);
    tNormMax = max(max(tNorm.*double(contribution>0),[],1),0);
    
    % CHOICE: The effective filter length should be at least half as long as the 
    % defined length
    m = tNormMin+0.5 > tNormMax;
    delta = (tNormMin(m)+0.5-tNormMax(m))/2;
    tNormMin(m) = tNormMin(m)-delta;
    tNormMax(m) = tNormMax(m)+delta;

    % Distribute the points to bins
    classification = min(round(bsxfun(@rdivide,bsxfun(@minus,tNorm.*double(contribution>0),tNormMin)*nBins,(tNormMax-tNormMin))+0.5),nBins);
    
    % Compute the response of the different orientations
    sz = size(contribution);
    [~,i2,i3,i4] = ind2sub(sz,1:numel(contribution));
    subs = {classification(:),i2',i3',i4'};
    response = accumarray(subs,contribution(:),[nBins,sz(2:4)]);
    response = squeeze(prod(response + minBinResponse,1));
    
    % Find the max response of all the different orientations
    [rowVal rowIdx] = max(response,[],1);
    [colVal colIdx] = max(rowVal);
    rowIdx = rowIdx(colIdx);
    phiMax = phiAngles(colIdx);
    thetaMax = thetaAngles(rowIdx);
    
    % Save the result
    obj_data_orientation(i,:) = [sind(thetaMax)*cosd(phiMax) sind(thetaMax)*sind(phiMax) cosd(thetaMax)];
    obj_data_magnitude(i) = colVal;
end

obj.data.orientation = obj_data_orientation;
obj.data.magnitude = obj_data_magnitude;

disp('Process: Orientation computed!');

end

