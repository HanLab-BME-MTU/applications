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
indxOn1 = find(ratesDensity1.rateOnPerClust(:,3)>=MIN_NUM);
indxOff1 = find(ratesDensity1.rateOffPerClust(:,3)>=MIN_NUM);
indxDen1 = find(ratesDensity1.densityPerClust(:,3)>=MIN_NUM);
indxOn2 = find(ratesDensity2.rateOnPerClust(:,3)>=MIN_NUM);
indxOff2 = find(ratesDensity2.rateOffPerClust(:,3)>=MIN_NUM);
indxDen2 = find(ratesDensity2.densityPerClust(:,3)>=MIN_NUM);

%thus define maximum cluster size for each measurement
numOn = max([indxOn1; indxOn2]);
numOff = max([indxOff1; indxOff2]);
numDen = max([indxDen1; indxDen2]);

%initialize parameter vectors and variance-covariance matrices
[paramVec1,paramVec2] = deal(NaN(numOn+numOff+numDen,1));
[paramCovMat1,paramCovMat2] = deal(NaN(numOn+numOff+numDen));

%fill the vectors
paramVec1(indxOn1) = ratesDensity1.rateOnPerClust(indxOn1,1);
paramVec1(numOn+indxOff1) = ratesDensity1.rateOffPerClust(indxOff1,1);
paramVec1(numOn+numOff+indxDen1) = ratesDensity1.densityPerClust(indxDen1,1);
paramVec2(indxOn2) = ratesDensity2.rateOnPerClust(indxOn2,1);
paramVec2(numOn+indxOff2) = ratesDensity2.rateOffPerClust(indxOff2,1);
paramVec2(numOn+numOff+indxDen2) = ratesDensity2.densityPerClust(indxDen2,1);

%fill the matrices
paramCovMat1(indxOn1,indxOn1)

%identify entries that are not NaN in at least one of the two vectors
indxNoNaN = find(~isnan(paramVec1) | ~isnan(paramVec2));

%remove these entries from parameter vectors
paramVec1 = paramVec1(indxNoNaN);
paramVec2 = paramVec2(indxNoNaN);

%expand parameter variance-covariance matrices to accommodate additional
%parameter rows if necessary



%% ~~~ the end ~~~