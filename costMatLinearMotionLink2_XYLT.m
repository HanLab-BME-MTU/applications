function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
    errFlag] = costMatLinearMotionLink2_XYLT(movieInfo,kalmanFilterInfoFrame1,...
    costMatParam,nnDistFeatures,probDim,prevCost,featLifetime)
%COSTMATLINEARMOTIONLINK2_XYLT provides a cost matrix for linking features
%based on competing linear motion models, asuming alongated features
%defined by x, y position, length and orientation.
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = costMatLinearMotionLink2_XYLT(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,nnDistFeatures,probDim,prevCost,featLifetime)
%
%INPUT  movieInfo             : A 2x1 array (corresponding to the 2 frames of
%                               interest) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features.
%                                   1st column for values and 2nd column
%                                   for standard deviations.
%             .length             : Length of features
%             .angle              : Orientation of features
%
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%      kalmanFilterInfoFrame1 : Structure with at least the following fields:
%             .stateVec           : Kalman filter state vector for each
%                                   feature in 1st frame.
%             .stateCov           : Kalman filter state covariance matrix
%                                   for each feature in 1st frame.
%             .noiseVar           : Variance of state noise for each
%                                   feature in 1st frame.
%      costMatParam           : Structure with fields:
%             .linearMotion       : 1 to propagate using a linear motion
%                                   model, 0 otherwise.
%             .minSearchRadius    : Minimum allowed search radius.
%             .maxSearchRadius    : Maximum allowed search radius.
%             .brownStdMult       : Factor multiplying Brownian
%                                   displacement std to get search radius.
%             .maxLengthRatio     : Maximum ratio between the length of two
%                                   features in two consecutive time points
%                                   that  allows linking them.
%             .maxAnglePenalty    : Maximum penalty between the angle of two
%                                   features in two consecutive time points
%                                   that  allows linking them.
%             .useLocalDensity    : Logical variable indicating whether to use
%                                   local density in search radius estimation.
%             .nnWindow           : Number of past frames for calculating
%                                   nearest neighbor distance.
%             .diagnostics        : Row vector indicating frames at which
%                                   histogram of linking distances (from
%                                   the beginning till that frame) are to
%                                   be plotted. Does not work for 1st or
%                                   last frame of a movie.
%                                   Optional. Default: None.
%      nnDistFeatures         : Matrix of nearest neighbor distances of
%                               features in first frame as well as of
%                               features in previous frames that they are
%                               connected to.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%      prevCost               : Structure with fields:
%             .all                : Matrix of previous linking costs.
%             .max                : Maximum previous linking cost.
%      featLifetime           : Lengths of tracks that features in
%                               first frame belong to.
%
%OUTPUT costMat               : Cost matrix.
%       propagationScheme     : Propagation scheme corresponding to each
%                               cost in the cost matrix.
%       kalmanFilterInfoFrame2: Structure with at least the following fields:
%             .stateVec           : Kalman filter prediction of state
%                                   vector in 2nd frame based on all 3
%                                   motion models.
%             .stateCov           : Kalman filter prediction of state
%                                   covariance matrix in 2nd frame based on
%                                   all 3 motion models.
%             .obsVec             : Kalman filter prediction of the
%                                   observed variables in 2nd frame based
%                                   on all 3 motion models.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%REMARKS Three competing linear motion models: 1, 2 and 3.
%        1: forward drift, 2: backward drift, 3: zero drift (Brownian).
%
% Khuloud Jaqaman, March 2007
% Sylvain Berlemont, April 2010

%% Output

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinearMotionLink2_XYLT')
    disp('--costMatLinearMotionLink2: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost function parameters
linearMotion = costMatParam.linearMotion;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult    = costMatParam.brownStdMult;
maxLengthRatio = costMatParam.maxLengthRatio;
maxAnglePenalty = costMatParam.maxAnglePenalty;
useLocalDensity = costMatParam.useLocalDensity;
nnWindow = costMatParam.nnWindow;

if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
end

if isfield(costMatParam,'diagnostics')
    diagnostics = costMatParam.diagnostics;
else
    diagnostics = 0;
end

%calculate nearest neighbor distance given feature history
frameNum = size(nnDistFeatures,2);
tmpNN = max(1,frameNum-nnWindow);
nnDistTracks = min(nnDistFeatures(:,tmpNN:end),[],2);

%% Motion propagation

%specify number of propagation schemes used
numSchemes = 3;

%calculate vector sizes
vecSize = 2 * probDim;

%construct transition matrices
if linearMotion
    transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    transMat(:,:,3) = eye(vecSize); %zero drift transition matrix
else
    transMat(:,:,3) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    transMat(:,:,1) = eye(vecSize); %zero drift transition matrix
end

%construct observation matrix
observationMat = [eye(probDim) zeros(probDim)]; %observation matrix

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%reserve memory for "kalmanFilterInfoframe2"
kalmanFilterInfoFrame2 = struct('stateVec',zeros(numFeaturesFrame1,vecSize,numSchemes),...
    'stateCov',zeros(vecSize,vecSize,numFeaturesFrame1,numSchemes),...
    'obsVec',zeros(numFeaturesFrame1,probDim,numSchemes));

%apply Kalman filters to each feature in 1st frame
for iFeature = 1 : numFeaturesFrame1
    
    %get state vector and its covariance matrix of feature in 1st frame
    stateOld = kalmanFilterInfoFrame1.stateVec(iFeature,:)';
    stateCovOld = kalmanFilterInfoFrame1.stateCov(:,:,iFeature);
    noiseVar = kalmanFilterInfoFrame1.noiseVar(:,:,iFeature);
    
    %go over all possible propagation schemes
    for iScheme = 1 : numSchemes
        
        %predict state vector of feature in 2nd frame
        stateVec = transMat(:,:,iScheme)*stateOld;
        
        %predict state covariance matrix of feature in 2nd frame
        stateCov = transMat(:,:,iScheme)*stateCovOld*transMat(:,:,iScheme)' ...
            + noiseVar;
        
        %determine observation vector of feature in 2nd frame (i.e. the
        %propagated position of the feature)
        obsVec = observationMat*stateVec;
        
        %save information in kalmanFilterInfoFrame2
        kalmanFilterInfoFrame2.stateVec(iFeature,:,iScheme) = stateVec';
        kalmanFilterInfoFrame2.stateCov(:,:,iFeature,iScheme) = stateCov;
        kalmanFilterInfoFrame2.obsVec(iFeature,:,iScheme) = obsVec';
    end
end

%get the propagated positions of features in 1st frame based on the three propagation schemes
propagatedPos = kalmanFilterInfoFrame2.obsVec;

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

%calculate the cost matrices for all three propagation schemes
for iScheme = 1 : numSchemes
    
    %put the propagated x and y coordinates of features from 1st frame in
    %one matrix
    coord1 = propagatedPos(:,:,iScheme);
    
    %calculate the distances between features
    costMatTmp(:,:,iScheme) = createDistanceMatrix(coord1,coord2);
    
end

%find the minimum cost for the link between every pair, which also
%determines the best propagation scheme to perform that link
[costMat,propagationScheme] = min(costMatTmp,[],3);

%% Search radius

%get the Kalman standard deviation of all features in frame 1
kalmanStd = sqrt(probDim * squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)));

%copy brownStdMult into vector
stdMultInd = repmat(brownStdMult,numFeaturesFrame1,1);

%if local density information is used to expand search radius ...
if useLocalDensity
    
    %divide each feature's nearest neighbor distance/closestDistScale by kalmanStd
    ratioDist2Std = nnDistTracks./kalmanStd/closestDistScale;
    
    %make ratios larger than maxStdMult equal to maxStdMult
    ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;
    
    %expand search radius multiplication factor if possible
    stdMultInd = max([stdMultInd ratioDist2Std],[],2);
    
end

%get the search radius of each feature in frame 1 and make sure it falls
%within reasonable limits
searchRadius = stdMultInd .* kalmanStd;
searchRadius(searchRadius>maxSearchRadius) = maxSearchRadius;
searchRadius(searchRadius<minSearchRadius) = minSearchRadius;

%replicate the search radius to compare to cost matrix
searchRadius = repmat(searchRadius,1,numFeaturesFrame2);

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%square the cost matrix to make the cost = distance squared
costMat = costMat.^2;

%% Length penalty
l1 = repmat(movieInfo(1).length(:,1),1,movieInfo(2).num);
l2 = repmat(movieInfo(2).length(:,1)',movieInfo(1).num,1);
lengthRatio = l1./l2;
ind = lengthRatio < 1;
if ~isempty(ind)
    lengthRatio(ind) = 1 ./ lengthRatio(ind);
end

%assign NaN to all pairs whose length ratio is larger than the
% maximum allowed
lengthRatio(lengthRatio > maxLengthRatio) = NaN;

%multiply the distance between pairs with the ratio between their
%length
costMat = costMat.*lengthRatio;

%% Angle penalty
theta1 = repmat(movieInfo(1).angle(:,1),1,movieInfo(2).num);
theta2 = repmat(movieInfo(2).angle(:,1)',movieInfo(1).num,1);
theta12 = min(abs(theta1 - theta2), abs(theta2 - theta1)); % theta12 in [0 +pi/2]
tt12 = tan(theta12) + 1; % tt12 in [1, +Inf]
tt12(tt12 > maxAnglePenalty) = NaN;

%multiply the distance between pairs with the penalty on their angle
costMat = costMat.*tt12;

%% Birth and death

%append matrix to allow birth and death
% jonas, 10/09: fix for non-sparse tracker
if isstruct(prevCost)
    prevCostMax = prevCost.max;
else
    prevCostMax = max(prevCost(:));
end

if ~isnan(prevCostMax) && prevCostMax ~= 0
    maxCost = 1.05*prevCostMax;
else
    maxCost = 1.05*max(prctile(costMat(:),100),eps);
    
end

deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
% costLR = min(min(min(costMat))-1,-1);
costLR = maxCost;
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;

%% Histogram of linking distances

%get current frame
% jonas, 10/09: fix for non-sparse tracker
if isstruct(prevCost)
    currentFrame = size(prevCost.all,2);
else
    currentFrame = size(prevCost,2);
end

%check whether current frame matches any of the diagnostics frames
if currentFrame ~= 1 && any(diagnostics == currentFrame)
    
    %get linking distances
    % jonas, 10/09: fix for non-sparse tracker
    if isstruct(prevCost)
        prevCostNoCol1 = prevCost.all(:,2:end);
    else
        prevCostNoCol1 = prevCost(:,2:end);
    end
    linkingDistances = sqrt(prevCostNoCol1(~isnan(prevCostNoCol1)));
    
    %plot histogram
    figure('Name',['frame # ' num2str(currentFrame)],'NumberTitle','off');
    try
        histogram(linkingDistances,[],0);
        xlabel('Linking distance');
        ylabel('Counts');
    catch
        disp('histogram plot failed');
    end
    
end


%% ~~~ the end ~~~
