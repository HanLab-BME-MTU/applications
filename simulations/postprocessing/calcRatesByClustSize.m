function [rateAssocPerClust,rateDissocPerClust,eventTable] = calcRatesByClustSize(sizeHist,timeStep)
%CALCRATESBYCLUSTSIZE calculates association and dissociation rates per
%cluster size.
%
%   INPUT:   
%       sizeHist:   A 2D array with each row corresponding to an
%                   association or dissociation event.  Column 1 gives the
%                   track number (as in compTracks.seqOfEvents), columns 2
%                   gives the cluster size, columns 3 and 4 give the
%                   cluster's starting and ending iteration points, column
%                   5 is the lifetime and column 6 indicates how the
%                   cluster ended - 1 = dissociation and 2 = association.
%       timeStep:   the time step used in the simulation
%
%   OUTPUT:
%       rateAssocPerClust:   A 1D array of calculated association rates
%                            with each row corresponding to a cluster size.
%       rateDissocPerClust:  A 1D array of calculated dissociation rates 
%                            with each row corresponding to a cluster size.
%       eventTable:          A 2D array with dimensions of rows equal to
%                            the maximum cluster size and 8 columns.
%                            For each cluster size (row), column 1 gives
%                            the total number of occurrences, column 2
%                            gives the mean lifetime, column 3 gives the
%                            total number of dissociation events with
%                            columns 4 and 5 giving the calculated
%                            probabilities and rates of dissociation,
%                            column 6 gives the total number of association
%                            events and columns 7 and 8 give the calculated
%                            probabilities and rates of association.
%
%   Robel Yirdaw, 08/29/13
%

    maxClusterSize = max(sizeHist(:,2));
    eventTable = NaN(maxClusterSize,8);       

    %Cluster size = 1
    currClustIndx_All = (sizeHist(:,2) == 1);
    %Total # of events
    eventTable(1,1) = sum(currClustIndx_All);
    %Overall lifetime
    eventTable(1,2) = mean(sizeHist(currClustIndx_All,5))*timeStep;
    %Association events
    currClustIndx_Assoc = ((sizeHist(:,2) == 1) & (sizeHist(:,6) == 2));
    %Total # of association events
    eventTable(1,6) = sum(currClustIndx_Assoc);            
    %Prob. of association
    eventTable(1,7) = eventTable(1,6)/eventTable(1,1);
    %Association rate
    eventTable(1,8) = eventTable(1,7)/eventTable(1,2);

    %Cluster sizes 2 and higher
    for clustIndx=2:maxClusterSize

        currClustIndx_All = (sizeHist(:,2) == clustIndx);
        %Total # of events
        eventTable(clustIndx,1) = sum(currClustIndx_All);
        %Overall lifetime
        eventTable(clustIndx,2) = mean(sizeHist(currClustIndx_All,5))*timeStep;

        %Dissociation events
        currClustIndx_Dissoc = ((sizeHist(:,2) == clustIndx) & (sizeHist(:,6) == 1));
        %Total # of dissociation events
        eventTable(clustIndx,3) = sum(currClustIndx_Dissoc);            
        %Prob. of dissociation
        eventTable(clustIndx,4) = eventTable(clustIndx,3)/eventTable(clustIndx,1);
        %Dissociation rate
        eventTable(clustIndx,5) = eventTable(clustIndx,4)/eventTable(clustIndx,2);

        %Association events
        currClustIndx_Assoc = ((sizeHist(:,2) == clustIndx) & (sizeHist(:,6) == 2));
        %Total # of association events
        eventTable(clustIndx,6) = sum(currClustIndx_Assoc);            
        %Prob. of association
        eventTable(clustIndx,7) = eventTable(clustIndx,6)/eventTable(clustIndx,1);
        %Association rate
        eventTable(clustIndx,8) = eventTable(clustIndx,7)/eventTable(clustIndx,2);
    end
    
    rateDissocPerClust = eventTable(:,5);
    rateAssocPerClust = eventTable(:,8);
    
    