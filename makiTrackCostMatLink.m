function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
     errFlag] = makiTrackCostMatLink(movieInfo,kalmanFilterInfoFrame1,...
     costMatParam,nnDistFeatures,probDim,prevCost,featLifetime)
%MAKITRACKCOSTMATLINK provides a cost matrix for frame-to-frame linking of HeLa kinetochores
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = makiTrackCostMatLink(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,nnDistFeatures,probDim,prevCost,featLifetime)
%
%INPUT  movieInfo             : A 2x1 array (corresponding to the 2 frames of 
%                               interest) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features. 
%                                   1st column for values and 2nd column 
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%             .kinType            : Kinetochore type: 0 - inlier, 1 -
%                                   unaligned, 2 - lagging.
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
%             .minSearchRadius    : Minimum allowed search radius. 3
%                                   entries for inlier, unaligned and
%                                   lagging kinetochores.
%             .maxSearchRadius    : Maximum allowed search radius. 3
%                                   entries for inlier, unaligned and
%                                   lagging kinetochores.
%             .brownStdMult       : Factor multiplying Brownian
%                                   displacement std to get search radius.
%             .lftCdf             : Lifetime cumulative density function.
%                                   Column vector, specifying cdf for
%                                   lifetime = 0 to movie length.
%                                   Enter [] if cdf is not to be used. 
%                                   Optional. Default: [].
%             .useLocalDensity    : Logical variable indicating whether to use
%                                   local density in search radius estimation.
%             .nnWindow           : Number of past frames for calculating
%                                   nearest neighbor distance.
%      nnDistFeatures         : Matrix of nearest neighbor distances of 
%                               features in first frame as well as of 
%                               features in previous frames that they are
%                               connected to.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%      prevCost               : Matrix of previous linking costs.
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
%Khuloud Jaqaman, July 2008

%% Output

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('makiTrackCostMatLink')
    disp('--makiTrackCostMatLink: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost function parameters
linearMotion = costMatParam.linearMotion;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult    = costMatParam.brownStdMult;
useLocalDensity = costMatParam.useLocalDensity;
nnWindow = costMatParam.nnWindow;
if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
end
if isfield('costMatParam','lftCdf')
    lftCdf = costMatParam.lftCdf;
else
    lftCdf = [];
end

%calculate nearest neighbor distance given feature history
frameNum = size(nnDistFeatures,2);
tmpNN = max(1,frameNum-nnWindow);
nnDistTracks = min(nnDistFeatures(:,tmpNN:end),[],2);

%get kinetochore types in first frame
kinType = movieInfo(1).kinType;
type0Indx = find(kinType==0); %inlier kinetochores
type1Indx = find(kinType==1); %unaligned kinetochores
type2Indx = find(kinType==2); %lagging kinetochores

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

%get the search radius of each feature in frame 1
searchRadius = stdMultInd .* kalmanStd;

%put the search radii into groups based on kinetochore type
searchRadiusType0 = searchRadius(type0Indx);
searchRadiusType1 = searchRadius(type1Indx);
searchRadiusType2 = searchRadius(type2Indx);

%make sure that search radii are within the appropriate limites
searchRadiusType0(searchRadiusType0>maxSearchRadius(1)) = maxSearchRadius(1);
searchRadiusType0(searchRadiusType0<minSearchRadius(1)) = minSearchRadius(1);
searchRadiusType1(searchRadiusType1>maxSearchRadius(2)) = maxSearchRadius(2);
searchRadiusType1(searchRadiusType1<minSearchRadius(2)) = minSearchRadius(2);
searchRadiusType2(searchRadiusType2>maxSearchRadius(3)) = maxSearchRadius(3);
searchRadiusType2(searchRadiusType2<minSearchRadius(3)) = minSearchRadius(3);

%write search radii back into one vector
searchRadius(type0Indx) = searchRadiusType0;
searchRadius(type1Indx) = searchRadiusType1;
searchRadius(type2Indx) = searchRadiusType2;

%replicate the search radius to compare to cost matrix
searchRadius = repmat(searchRadius,1,numFeaturesFrame2);

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%square the cost matrix to make the cost = distance squared
costMat = costMat.^2;

%% Lifetime penalty

if ~isempty(lftCdf)

    %specify 1 - lifetime cumulative probability
    oneMinusLftCdf = 1 - lftCdf;

    %calculate 1 / (lifetime penalty), which is 1 / (1-cumulative probability
    %of lifetime of feature in first frame)
    oneOverLftPen = oneMinusLftCdf(featLifetime+1);

    %multiple each cost by the lifetime penalty
    costMat = costMat ./ repmat(oneOverLftPen,1,numFeaturesFrame2);

    %replace infinite costs by NaN
    costMat(isinf(costMat)) = NaN;
    
end

%% Birth and death

%append matrix to allow birth and death
if any(~isnan(prevCost(:)))
    maxCost = 1.05*max(prevCost(:));
else
    maxCost = prctile(costMat(:),80);
end
deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;


%% ~~~ the end ~~~


