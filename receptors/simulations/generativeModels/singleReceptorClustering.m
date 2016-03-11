function aggregFlag = singleReceptorClustering(aggregFlagIN,clusterMembers,receptor2clusterPrev)
%SINGLERECEPTORCLUSTERING prevents the clustering of more than two
%receptors with each other or with an exisiting to cluster. This allows
%clusters to grow by one receptor only.  Also prevents dissociation and
%association events from occuring in the same iteration.
%   Note that in the overall process of clustering (aggregation), a cluster 
%   is created first and then, based on the values in aggregFlag, those
%   receptors who are not supposed to aggregate will be removed from the
%   cluster that was just formed.  The goal here is to impose more
%   restriction on the aggregation in that a cluster (including of size 1)
%   can only increase by size of 1, by one receptor.
%   Modification on 09/12/13 handles the case where a dissociation and
%   association within the same iteration is being attempted.  In this
%   case, the aggregFlag for all receptors that were involved in the
%   dissociation event will be set to 0 - they are not allowed to undergo
%   association in this iteration.  However, if there is a cluster leftover
%   from the dissociation event, the flag needs to be reset back to 1 so
%   that it is not dissociated in receptorAggregationAlg, upon return from
%   the call here.
%   
%   INPUT:   
%       aggregFlagIN:   1D array of booleans for each receptor in the cluster
%                       indicating whether the receptor can aggregate
%       clusterMembers: 1D array of receptors that make up the cluster
%       receptor2clusterPrev:   1D array of complete receptor to cluster
%                               associations for the previous iteration
%
%   OUTPUT:
%       aggregFlag: the modified (or not) array of receptor aggregation
%                   boolean values
%
%   Robel Yirdaw, 09/05/13
%       Modified, 09/12/13
%
      
        
    %Copy the aggregation flag
    aggregFlag = aggregFlagIN;

    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Determine the previous clustering state of all the members            
        prevSize = NaN(length(aggregFlag),1);
        %For each receptor in clusterMembers with aggregFlag values
        for memIter=1:length(aggregFlag)
            %index of receptors in receptor2clusterPrev that make the
            %current cluster
            memberIndx = clusterMembers(memIter);            
            %Size, in previous iteration point, of the current
            %receptor's cluster
            prevSize(memIter,1) = ...
                numel(find(receptor2clusterPrev == receptor2clusterPrev(memberIndx,1)));
        end            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %If this new cluster is composed of a previous cluster with flag 0,
        %reset the flag and return - else continue to process.
        if (any(aggregFlag(prevSize > 1) == 0))
            %If a particle has size > 1 with flag = 0, it is a cluster
            %forbidden to associate because it participated in a
            %dissociation event in the current iteration - set the flag
            %back to 1 so that receptorAggregatioAlg does not dissociate it
            %- receptors can be kept at 0.
            aggregFlag(prevSize > 1) = 1;
            aggregFlag(prevSize == 1) = 0;
        else
            %Get indicies and number of those receptors able to aggregate
            aggregRecepIndx = find(aggregFlag == 1);
            numMembersAggreg = numel(aggregRecepIndx);
            %Previous sizes for those with aggregation flags of 1 only
            prevSizeAggreg = prevSize(aggregRecepIndx);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This is the change in size of the clusters of the *relevant* members
            sizeChange = numMembersAggreg - prevSizeAggreg;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( sum(sizeChange > 1) > 1 ) 
                %Index for free receptors that are associating
                %freeRecepIndx = find(prevSize == 1);
                freeRecepIndx = aggregRecepIndx(prevSizeAggreg == 1);
                %Currently, receptors are able to dissociate from a cluster
                %and associate back all within the same iteration - the
                %receptors are transiently free but receptor2clusterPrev
                %will show them clustered.  Addressed this in main code -
                %testing.
                if (~isempty(freeRecepIndx))
                    %Reset flag for free receptors that are associating
                    aggregFlag(freeRecepIndx,1) = 0;            
                    %Pick only one free receptor to aggregate
                    recept2associate(1) = randi(numel(freeRecepIndx),1);
                    %Set the flag for picked receptor to 1.
                    %aggregFlag(aggregRecepIndx(freeRecepIndx(recept2associate(1)))) = 1;
                    aggregFlag(freeRecepIndx(recept2associate(1))) = 1;
                    if (all(prevSizeAggreg == 1))
                        %This is an initial cluster formation - allow two
                        %receptors to associate          
                        %Will randomly pick a receptor from those remaining so
                        %first exclude the one selected above.
                        freeRecepIndx(recept2associate(1)) = [];
                        %Pick the second receptor
                        recept2associate(2) = randi(numMembersAggreg-1,1); 
                        %Set the flag for picked receptor to 1.
                        aggregFlag(freeRecepIndx(recept2associate(2))) = 1;
                    end

                end %If freeRecepIndx not empty

            end %If sum of >1 size changes
            
        end %If clusters with aggregFlag = 0
        
    catch exception
        fprintf('\n\nError: singleReceptorClustering\n');
        exception.message
        for k=1:length(exception.stack)
            exception.stack(k)
        end
        
    end %try/catch
    
end
