function [rateOffPerClust] = clusterOnOffRatesAndDensity(sizeHist,timeStep)
%clusterOnOffRatesAndDensity calculates cluster on and off rates and densities
%
%   INPUT:   
%       sizeHist:   A 2D array with each row corresponding to an
%                   association or dissociation event.  Column 1 gives the
%                   segment number (as in compTracks.seqOfEvents), column 2
%                   gives the cluster size, columns 3 and 4 give the
%                   cluster's starting and ending frames, column
%                   5 is the lifetime (in frames) and column 6 indicates how the
%                   cluster ended - 1 = dissociation and 2 = association.
%       timeStep:   the time step between frames
%
%   OUTPUT:
%       rateAssocPerClust:   A 1D array of calculated association rates
%                            with each row corresponding to a cluster size.
%       rateDissocPerClust:  A 1D array of calculated dissociation rates 
%                            with each row corresponding to a cluster size.
%       eventTable:          A 2D array with dimensions of rows equal to
%                            the maximum cluster size and 8 columns.
%                            For each cluster size (row), columns give:
%                            1: total number of occurrences
%                            2: mean lifetime
%                            3: total number of dissociation events
%                            4/5: calculated probabilities/rates of dissociation
%                            6: total number of association events
%                            7/8: calculated probabilities/rates of association
%
%       eventTable_mono:     A 2D array with number of rows equal to the
%                            maximum cluster size and 3 columns containing
%                            information for monomer association events. 
%                            1: number of associations for monomers to
%                            cluster of size row+1
%                            2: probability of the transition
%                            3: calculated association rate
%
%   Khuloud Jaqaman, May 2015


%get maximum cluster size
maxClusterSize = max(sizeHist(:,2));
rateOffPerClust = NaN(maxClusterSize,1);

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size, and whether they ended by
    %association or dissociation
    %only look at clusters with known start and end time
    indxClust = find(sizeHist(:,2)==iSize&~isnan(sizeHist(:,5)));
    clustLft = sizeHist(indxClust,5);
    clustEndType = sizeHist(indxClust,6);
    
    %calculate dissociation rate
    rateOffPerClust(iSize) = ...
        (length(find(clustEndType==1))/length(clustEndType)) / ...
        (mean(clustLft)*timeStep);
    
end
    


