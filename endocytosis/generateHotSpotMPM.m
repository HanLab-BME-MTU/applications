function [mpm] = generateHotSpotMPM(clusterResults);

% generateHotSpotMPM converts results from QTcluster into an MPM structure
% that contains mpm's for the cluster centroids, clustered initiations, and
% unclustered initiations.
%
% INPUT clusterResults: output of QTcluster. Must contain .clusterCentroids 
%                       and.cluster Results
%
% OUTPUT mpm.mpmCentroids: [xpos ypos] for all centroids (first colum is x
%           position, second colum is y position; each entry in mpm matrix 
%           is a different movie)
%        mpm.mpmClustered: [xpos ypos] for all clustered initiations
%        mpm.mpmUnClustered: [xpos ypos] for all UNclustered initiations
%
% Uses: none
%
% Daniel Nunez, January 9, 2009

%for each movie clustering result
for iexp = 1:length(clusterResults)
    
    mpm(iexp).mpmCentroids = clusterResults(iexp).clusterCentroids;
    mpm(iexp).mpmClusteredInitiations = clusterResults(iexp).clusterResults(clusterResults(iexp).clusterResults(:,3)~=0,1:2);
    mpm(iexp).mpmUnClusteredInitiations = clusterResults(iexp).clusterResults(clusterResults(iexp).clusterResults(:,3)==0,1:2);
    
end %of for each movie

end %of function