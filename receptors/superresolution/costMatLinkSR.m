function [costMat,nonlinkMarker,errFlag] = costMatLinkSR(movieInfo,searchRadius)
%costMatLinkSR provides a cost matrix for linking static features for super-resolution applications
%
%SYNOPSIS [costMat,nonlinkMarker,errFlag] = costMatLinkSR(movieInfo,searchRadius)
%
%INPUT  movieInfo             : A 2x1 array containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features.
%                                   1st column for values and 2nd column
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%      searchRadius           : Search radius for linking features.
%
%OUTPUT costMat               : Cost matrix.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, August 2011

%% Output

costMat = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinkSR')
    disp('--costMatLinkSR: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%% Inter-particle distance calculation

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%put the coordinates of features in the 1st frame in one matrix
coord1 = movieInfo(1).allCoord(:,1:2:end);

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

%calculate the distances between features
costMat = createDistanceMatrix(coord1,coord2);

%% Search radius

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%square the cost matrix to make the cost = distance squared
costMat = costMat.^2;

%% Birth and death

maxCost = 2*max(max(costMat(:)),eps);

deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(costMat(:))-1,-1);
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
