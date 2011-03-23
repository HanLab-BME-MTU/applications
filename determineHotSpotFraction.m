function [experiment] = determineHotSpotFraction(experiment);

% plotPairCorrelation calculates the density of pits defined by rest
% and that fall within an inputMask
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
%                       .status (optional), which contains a number for
%                       each pit which identifies it as belonging to a
%                       certain population, so that all pits that belong to
%                       population1 have a status value 1 (this is useful
%                       in calculating pair correlations for the different
%                       populations)
%
%
% OUTPUT
%           .hotSpotFraction =  pits found in hotspots divided by all pits 
%               analyzed by QTcluster function
% Uses:
%       addFieldClusterResults
%
% Daniel Nunez, updated Maech 23, 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add clustering results
experiment = addFieldClusterResults(experiment);


for iexp = 1:length(experiment)
    %store pair in each experiment structure
    experiment(iexp).hotSpotFraction = size(clusterResults.clusterResults(),1)...
        /size(clusterResults.clusterResults,1);    
end

end %of function