function [rateAssocPerClust,rateDissocPerClust,eventTable,eventTable_mono] =...
    calcRatesByClustSize(sizeHist,timeStep)
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
%       eventTable_mono:     A 2D array wiht number rows equal to the
%                            maximum cluster size and 3 columns containing
%                            information for monomer association events. 
%                            Column 1 gives the number of associations for 
%                            monomers to cluster of size row+1. Column 2
%                            gives the probability of the transition and
%                            column 3 gives the calculated association rate.
%
%   Robel Yirdaw, 08/29/13
%       Modified, 02/20/14
%

    maxClusterSize = max(sizeHist(:,2));
    eventTable = NaN(maxClusterSize,8);       

    %Index of all cluster size = 1
    %currClustIndx_All = (sizeHist(:,2) == 1);     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %021714
    %Get index of dimerization events - if both monomers listed, the event
    %should be counted as one. Column 8 in sizeHist must be set by
    %cluserHistoryFromCompTracks_aggregState.
    %indx_monoDuplEvents = ~isnan(sizeHist(:,8));
    %Total # of events
    %eventTable(1,1) = sum(currClustIndx_All) - sum(indx_monoDuplEvents);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %022014
    %Third way is to count only those events where both monomers have
    %start time. Column 8 in sizeHist must be set by
    %cluserHistoryFromCompTracks_aggregState.
    currClustIndx_All = ((sizeHist(:,2) == 1) & (sizeHist(:,8) == 1)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Total # of events
    eventTable(1,1) = sum(currClustIndx_All);
    %Overall lifetime
    eventTable(1,2) = mean(sizeHist(currClustIndx_All,5))*timeStep;
    
    %The following must correspond with what is done above.
    %Association events
    %currClustIndx_Assoc = ((sizeHist(:,2) == 1) & (sizeHist(:,6) == 2));
    %... counting dimerization as a single event
    %currClustIndx_Assoc = ((sizeHist(:,2) == 1) & (sizeHist(:,6) == 2) &...
    %    isnan(sizeHist(:,8)));
    %... counting dimerzation events where both monomers have start time -
    %both monomers are counted.
    currClustIndx_Assoc = ((sizeHist(:,2) == 1) & (sizeHist(:,6) == 2) &...
        (sizeHist(:,8) == 1));
    
    %Total # of association events
    eventTable(1,6) = sum(currClustIndx_Assoc);            
    %Prob. of association
    eventTable(1,7) = eventTable(1,6)/eventTable(1,1);
    %Association rate
    eventTable(1,8) = eventTable(1,7)/eventTable(1,2);

    %Cluster sizes 2 and higher
    for clustIndx=2:maxClusterSize
        
        %The following must correspond with what is done above.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Indx of all events for current cluster size
        % = (sizeHist(:,2) == clustIndx);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Indx of all events for current cluster size where associations are
        %only for those where the start time exists for both segments.
        %Need to enable the same below for association part.
        currClustIndx_All = ( (sizeHist(:,2) == clustIndx) &...
            ( (sizeHist(:,6) == 1) | ( (sizeHist(:,6) == 2) & (sizeHist(:,8) == 1) ) ) );
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Association events
        %currClustIndx_Assoc = ((sizeHist(:,2) == clustIndx) & (sizeHist(:,6) == 2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Association events for current cluste size  where the start time 
        %exists for both segments.
        currClustIndx_Assoc = ((sizeHist(:,2) == clustIndx) & (sizeHist(:,6) == 2)...
            & (sizeHist(:,8) == 1));
        %Total # of association events
        eventTable(clustIndx,6) = sum(currClustIndx_Assoc);            
        %Prob. of association
        eventTable(clustIndx,7) = eventTable(clustIndx,6)/eventTable(clustIndx,1);
        %Association rate
        eventTable(clustIndx,8) = eventTable(clustIndx,7)/eventTable(clustIndx,2);
    end
    
    rateDissocPerClust = eventTable(:,5);
    rateAssocPerClust = eventTable(:,8);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate association rates for cluster size 1 (monomer) per size of
    %associating partner
    
    %Cluster size = 1 events
    %clust1Indx_All = (sizeHist(:,2) == 1);    
    %... counting dimerization as a single event
    %clust1Indx_All = ((sizeHist(:,2) == 1) & isnan(sizeHist(:,8)));    
    %... counting dimerzation events where both monomers have start time -
    %both monomers are counted.    
    clust1Indx_All = ((sizeHist(:,2) == 1) & (sizeHist(:,6) == 2) &...
            (sizeHist(:,8) == 1));    
        
    %Sizes of resulting clusters in the above events
    minFinalSize = min(sizeHist(clust1Indx_All,7));
    maxFinalSize = max(sizeHist(clust1Indx_All,7));    
    if (minFinalSize == 1)
        fprintf('\n================================================================');
        fprintf('\nFound %d monomer association events with final size of 1.',...
            sum(sizeHist(clust1Indx_All,7) == 1));
        fprintf('\nStarting with size 2.');
        fprintf('\n================================================================\n');   
        minFinalSize = 2;
    end
    if (maxFinalSize == 1)
        %In this case the for-loop below does not execute.
        fprintf('\n=====================================================================');
        fprintf('\nMax final size for monomer association events (%d) is 1.',...
            sum(sizeHist(clust1Indx_All,7) == 1));        
        fprintf('\n=====================================================================\n');   
    end
    
    
    %Initialize eventTable for monomers
    eventTable_mono = NaN(maxFinalSize-1,3);
    
    for finalSizeIndx=minFinalSize:maxFinalSize
        %Number of associations with cluster of size finalSizeIndx-1
        eventTable_mono(finalSizeIndx-1,1) =...
            sum(sizeHist(clust1Indx_All,7) == finalSizeIndx);
        %Probability of these events calculated with respect to the total
        %number of events for monomers determined above (eventTable(1,1))
        eventTable_mono(finalSizeIndx-1,2) =...
            eventTable_mono(finalSizeIndx-1,1)/eventTable(1,1);
        %The rate of association with cluster of size finalSizeIndx-1
        %calculated using the overall lifetime of monomers determined above
        %(eventTable(1,2)
        eventTable_mono(finalSizeIndx-1,3) =...
            eventTable_mono(finalSizeIndx-1,2)/eventTable(1,2);
    end
    
    clear clust1Indx_All minFinalSize maxFinalSize finalSizeIndx


    