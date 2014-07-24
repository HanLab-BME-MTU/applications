function [hotSpotRadius] = determineHotSpotRadius(experiment,rest);

% determineHotSpotRadius calculates radius of hot spot based on pair
% correlation (takes mean value of pair correlation from pixels 10 to 20
% and finds the radius at which the pair correlation is equal to this value
% plus one)
%
% INPUT:   experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the 
%                       function loadConditionData, and it needs to contain
%                       at least the field 
%                       .source, which is the (path) location of the 
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%                       
% OUTPUT
%           hotSpotRadius: mean value of pair correlation from pixels 10 to   
%           20 plus one)
%
% Uses: RipleysKfunction
%       plotPairCorrelation
%
% Daniel Nunez, January 20, 2009
% Daniel Nunez, March 24, 2011

dist = 1:20;
experiment = plotPairCorrelation(experiment,dist);

%make average pairCorrelation for exp (NOTE THAT THIS MEANS THAT THE RADIUS
%WILL CHANGE IF SET OF MOVIES IS CHANGED)
pitDen = mean([experiment.pairCorrelation],2);

%fit spline to productive data
pp = csaps(dist,pitDen);
denFitValues = fnplt(pp);
%find max value of fit
maxNormDen = max(denFitValues(2,:));
%find density value where radius is 10 pixels
findMidRad = find(denFitValues(1,:) >= 10,1,'first');
%take the mean density from 10 pixels to the end as the far distance value
farDen = mean(denFitValues(2,findMidRad:end));
%get point at 10 percent of difference in between maximum and long distance
%densities
findRad = find(denFitValues(2,:) < 1+farDen,1,'first');
hotSpotRadius = denFitValues(1,findRad);

end %of function