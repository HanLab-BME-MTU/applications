function [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat] = ...
    combineClusterRatesAndDensity(ratesDensityPerMovie)
%COMBINECLUSTERRATESANDDENSITY combines cluster on and off rates and densities for a group of equivalent movies
%
%   SYNOPSIS:
%       [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat] = ...
%           combineClusterRatesAndDensity(ratesDensityPerMovie)
%
%   INPUT:   
%       ratesDensityPerMovie: Cell array of output of
%                            clusterOnOffRatesAndDensity for each movie in
%                            group.
%
%   OUTPUT:
%       rateOnPerClust     : Mean, std and number of data points of on rate
%                            for clusters of size 1, 2, 3, ...
%       rateOffPerClust    : Mean, std and number of data points of off
%                            rate for clusters of size 1, 2, 3, ...
%       densityPerClust    : Mean, std and number of data points of density
%                            for clusters of size 1, 2, 3, ...
%       paramVarCovMat     : Full parameter variance-covariance matrix
%                            where parameter order is on rate, off rate,
%                            then density.
%
%   Khuloud Jaqaman, June 2015

%% Input

%get number of movies in group
numMovies = length(ratesDensityPerMovie);

%get maximum cluster size in group of movies
maxSizeRates = zeros(numMovies,1);
maxSizeDensity = zeros(numMovies,1);
for iMovie = 1 : numMovies
    maxSizeRates(iMovie) = length(ratesDensityPerMovie{iMovie}.rateOnPerClust);
    maxSizeDensity(iMovie) = length(ratesDensityPerMovie{iMovie}.densityPerClust);
end
maxSizeRatesAll = max(maxSizeRates);
maxSizeDensityAll = max(maxSizeDensity);

%% Calculation

%reserve memory
paramMatrix = [zeros(maxSizeRatesAll,numMovies); ... %on rate (0 if cluster not observed)
    NaN(maxSizeRatesAll,numMovies); ... %off rate (NaN if cluster not observed)
    zeros(maxSizeDensityAll,numMovies)]; %density (0 if cluster not observed)

%collect values for each movie
for iMovie = 1 : numMovies
    
    %on rates
    paramMatrix(1:maxSizeRates(iMovie),iMovie) = ratesDensityPerMovie{iMovie}.rateOnPerClust;
    
    %off rates
    paramMatrix(maxSizeRatesAll+1:maxSizeRatesAll+maxSizeRates(iMovie),iMovie) = ratesDensityPerMovie{iMovie}.rateOffPerClust;
    
    %densities
    paramMatrix(2*maxSizeRatesAll+1:2*maxSizeRatesAll+maxSizeDensity(iMovie),iMovie) = ratesDensityPerMovie{iMovie}.densityPerClust;
    
end

%calculate mean, number of data points, variance-covariance matrix and
%standard deviation
paramMean = nanmean(paramMatrix,2);
paramNumMov = sum( (~isnan(paramMatrix) & paramMatrix~=0), 2);
iRowNaN = find(isnan(max(paramMatrix,[],2)));
matTmp = paramMatrix;
matTmp(iRowNaN,:) = 0;
paramVarCovMat = nancov(matTmp');
paramVarCovMat(iRowNaN,:) = NaN;
paramVarCovMat(:,iRowNaN) = NaN;
paramStd = sqrt(diag(paramVarCovMat));

%% Output
paramCombined = [paramMean paramStd paramNumMov];
rateOnPerClust = paramCombined(1:maxSizeRatesAll,:);
rateOffPerClust = paramCombined(maxSizeRatesAll+1:2*maxSizeRatesAll,:);
densityPerClust = paramCombined(2*maxSizeRatesAll+1:end,:);

%% ~~~ the end ~~~