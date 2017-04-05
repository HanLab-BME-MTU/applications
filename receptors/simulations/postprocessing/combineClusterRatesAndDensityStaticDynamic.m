function [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat,paramMatrix] = ...
    combineClusterRatesAndDensityStaticDynamic(ratesDensityPerMovie,systemState)
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
%       paramMatrix        : intermediate statistics array, where the rows 
%                            follow the sequence: on rate, off rate and density
%                            for clusters of size 1, 2, 3, .... and the columns
%                            are the # of movie (simulation).
%                                 
%   Khuloud Jaqaman, June 2015
%   Luciana de Oliveira, August 2016, February 2017.

%% Input

%get number of movies in group
numMovies = length(ratesDensityPerMovie);
maxSizeDensity = zeros(numMovies,1);

%test if it is static or dynamic data

if systemState==1
%% Dynamic Calculation

%get maximum cluster size in group of movies
maxSizeRates = zeros(numMovies,1);
for iMovie = 1 : numMovies
    maxSizeRates(iMovie) = length(ratesDensityPerMovie{iMovie}.rateOnPerClust);
    maxSizeDensity(iMovie) = length(ratesDensityPerMovie{iMovie}.densityPerClust);
end
maxSizeRatesAll = max(maxSizeRates);
maxSizeDensity = max(maxSizeDensity);

%Only use values coming from at least 5 datapoints
MIN_CLUST = 5;

%% Calculation

%Modification Luciana de Oliveira, August 2016
% To have paramMatrix with size multiple of 3, now the calculations follow
% the density size instead on/off sizes. When the sizes are different it fills 
% NaNs in the empty rows.


%reserve memory
paramMatrix = [NaN(maxSizeDensity,numMovies); ... %on rate (NaN if cluster not observed)
    NaN(maxSizeDensity ,numMovies); ... %off rate (NaN if cluster not observed)
    zeros(maxSizeDensity,numMovies)]; %density (0 if cluster not observed)

%collect values for each movie
for iMovie = 1 : numMovies
    
    %get number of clusters used for off rate calculation
    numClust = ratesDensityPerMovie{iMovie}.numClustForRateCalc(:,1);
    
    %identify rates calculated from less than MIN_CLUST clusters
    %these rates will be ignored
    indxBad = find(numClust < MIN_CLUST);
    
    %on rates
    tmp = ratesDensityPerMovie{iMovie}.rateOnPerClust;
    tmp(indxBad) = NaN;
    paramMatrix(1:maxSizeRates(iMovie),iMovie) = tmp;
    
    
    %off rates
    tmp = ratesDensityPerMovie{iMovie}.rateOffPerClust;
    tmp(indxBad) = NaN;
    paramMatrix(maxSizeDensity+1:maxSizeDensity+maxSizeRates(iMovie),iMovie) = tmp;
    
    %densities
    paramMatrix(2*maxSizeDensity+1:2*maxSizeDensity+maxSizeDensity(iMovie),iMovie) = ratesDensityPerMovie{iMovie}.densityPerClust;
    
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
rateOnPerClust = paramCombined(1:maxSizeDensity,:);
rateOffPerClust = paramCombined(maxSizeDensity+1:2*maxSizeRatesAll,:);
densityPerClust = paramCombined(2*maxSizeDensity+1:end,:);

elseif systemState==0
%% static calculations

for iMovie = 1 : numMovies
    maxSizeDensity(iMovie) = length(ratesDensityPerMovie{iMovie}.densityPerClust);
end
maxSizeDensity = max(maxSizeDensity);
       
% Allocate space for the matrix with the densities for all 
       
       matrixDensity=zeros(maxSizeDensity,numMovies);
       
       % replace values for each movie
       
        for iMovie = 1 : numMovies     
  matrixDensity(1:length(ratesDensityPerMovie{iMovie}.densityPerClust),iMovie) = ratesDensityPerMovie{iMovie}.densityPerClust;
        end
paramMatrix=matrixDensity;   
rateOnPerClust=[];
rateOffPerClust=[];
densityPerClust=[];
paramVarCovMat=[];

end
%% ~~~ the end ~~~