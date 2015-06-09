function [distMaha,pValue] = compareClusterRatesAndDensityPair(ratesDensity1,ratesDensity2)
%COMPARECLUSTERRATESANDDENSITYPAIR compares cluster on and off rates and densities between two conditions
%
%   SYNOPSIS:
%       [testStat,pValue] = compareClusterRatesAndDensityPair(ratesDensity1,ratesDensity2)
%
%   INPUT:   
%       ratesDensity1: Structure of first condition, with fields
%                      rateOnPerClust, rateOffPerClust, densityPerClust and
%                      paramVarCovMat, as output by
%                      combineClusterRatesAndDensity.
%       ratesDensity2: Same but for second condition.
%
%   OUTPUT:
%       distMaha     : Mahalanobis distance between the two conditions.
%       pValue       : P-value from statistical test comparing the two
%                      conditions.
%
%   Khuloud Jaqaman, June 2015

%% Calculation

%define minimum number of datapoints per measurement (i.e. number of movies
%contributing to measurement)
MIN_NUM = 5;

%get indices of good measurements
indxGoodOn1 = find(ratesDensity1.rateOnPerClust(:,3)>=MIN_NUM);
indxGoodOff1 = find(ratesDensity1.rateOffPerClust(:,3)>=MIN_NUM);
indxGoodDen1 = find(ratesDensity1.densityPerClust(:,3)>=MIN_NUM);
indxGoodOn2 = find(ratesDensity2.rateOnPerClust(:,3)>=MIN_NUM);
indxGoodOff2 = find(ratesDensity2.rateOffPerClust(:,3)>=MIN_NUM);
indxGoodDen2 = find(ratesDensity2.densityPerClust(:,3)>=MIN_NUM);

% %get indices of bad measurements
% indxBadOn1 = find(ratesDensity1.rateOnPerClust(:,3)<MIN_NUM);
% indxBadOff1 = find(ratesDensity1.rateOffPerClust(:,3)<MIN_NUM);
% indxBadDen1 = find(ratesDensity1.densityPerClust(:,3)<MIN_NUM);
% indxBadOn2 = find(ratesDensity2.rateOnPerClust(:,3)<MIN_NUM);
% indxBadOff2 = find(ratesDensity2.rateOffPerClust(:,3)<MIN_NUM);
% indxBadDen2 = find(ratesDensity2.densityPerClust(:,3)<MIN_NUM);

%define maximum cluster size for each measurement
numOn1 = size(ratesDensity1.rateOnPerClust,1);
numOn2 = size(ratesDensity2.rateOnPerClust,1);
numOn = max([numOn1 numOn2]);
numOff1 = size(ratesDensity1.rateOffPerClust,1);
numOff2 = size(ratesDensity2.rateOffPerClust,1);
numOff = max([numOff1 numOff2]);
numDen1 = size(ratesDensity1.densityPerClust,1);
numDen2 = size(ratesDensity2.densityPerClust,1);
numDen = max([numDen1 numDen2]);
numTot1 = numOn1 + numOff1 + numDen1;
numTot2 = numOn2 + numOff2 + numDen2;
numTot = numOn + numOff + numDen;

%initialize parameter vectors and variance-covariance matrices
[paramVec1,paramVec2] = deal(NaN(numOn+numOff+numDen,1));

%fill the vectors
paramVec1(indxGoodOn1) = ratesDensity1.rateOnPerClust(indxGoodOn1,1);
paramVec1(numOn+indxGoodOff1) = ratesDensity1.rateOffPerClust(indxGoodOff1,1);
paramVec1(numOn+numOff+indxGoodDen1) = ratesDensity1.densityPerClust(indxGoodDen1,1);
paramVec2(indxGoodOn2) = ratesDensity2.rateOnPerClust(indxGoodOn2,1);
paramVec2(numOn+indxGoodOff2) = ratesDensity2.rateOffPerClust(indxGoodOff2,1);
paramVec2(numOn+numOff+indxGoodDen2) = ratesDensity2.densityPerClust(indxGoodDen2,1);

%find NaN entries to make the same in variance-covariance matrix
indxNaN1 = find(isnan(paramVec1));
indxNaN2 = find(isnan(paramVec2));

%fill the matrices
paramCovMat1 = ratesDensity1.paramVarCovMat;
paramCovMat1 = [paramCovMat1(:,1:numOn1) NaN(numTot1,numOn-numOn1) ...
    paramCovMat1(:,numOn1+1:numOn1+numOff1) NaN(numTot1,numOff-numOff1) ...
    paramCovMat1(:,numOn1+numOff1+1:numTot1) NaN(numTot1,numDen-numDen1)];
paramCovMat1 = [paramCovMat1(1:numOn1,:); NaN(numOn-numOn1,numTot); ...
    paramCovMat1(numOn1+1:numOn1+numOff1,:); NaN(numOff-numOff1,numTot); ...
    paramCovMat1(numOn1+numOff1+1:numTot1,:); NaN(numDen-numDen1,numTot)];
paramCovMat1(indxNaN1,:) = NaN;
paramCovMat1(:,indxNaN1) = NaN;
paramCovMat2 = ratesDensity2.paramVarCovMat;
paramCovMat2 = [paramCovMat2(:,1:numOn2) NaN(numTot2,numOn-numOn2) ...
    paramCovMat2(:,numOn2+1:numOn2+numOff2) NaN(numTot2,numOff-numOff2) ...
    paramCovMat2(:,numOn2+numOff2+1:numTot2) NaN(numTot2,numDen-numDen2)];
paramCovMat2 = [paramCovMat2(1:numOn2,:); NaN(numOn-numOn2,numTot); ...
    paramCovMat2(numOn2+1:numOn2+numOff2,:); NaN(numOff-numOff2,numTot); ...
    paramCovMat2(numOn2+numOff2+1:numTot2,:); NaN(numDen-numDen2,numTot)];
paramCovMat2(indxNaN2,:) = NaN;
paramCovMat2(:,indxNaN2) = NaN;

%identify entries that are not NaN in both vectors
indxNoNaN = find(~isnan(paramVec1) & ~isnan(paramVec2));

%remove these entries from parameter vectors and matrices
paramVec1 = paramVec1(indxNoNaN);
paramVec2 = paramVec2(indxNoNaN);
paramCovMat1 = paramCovMat1(:,indxNoNaN);
paramCovMat1 = paramCovMat1(indxNoNaN,:);
paramCovMat2 = paramCovMat2(:,indxNoNaN);
paramCovMat2 = paramCovMat2(indxNoNaN,:);

%calculate Mahalanobis distance between parameter vectors
distMaha = ((paramVec2-paramVec1)'/(paramCovMat1+paramCovMat2))*(paramVec2-paramVec1);


%% ~~~ the end ~~~