function [clusterState,merge2splitTime] = DissRateTestWithAggregVarRates(associationRate,dissociationRate,numReceptors,simTime,timeStep)
%DISSRATETESTWITHAGGREG simulates the aggregation into and dissociation
%from a single cluster of a specified number of receptors.
%
%   The simulation starts with all receptors free.  The user-defined
%   aggregation probability will be used to aggregate free receptors *one
%   at a time* to form a single cluster. Similarly, a dissociation
%   probability, calculated from user-defined dissociation rate, will be
%   used to dissociate receptors from the cluster, *one at a time*.
%   Simulation will run for a fixed number of iterations as determined by
%   simTime and timeStep.
%
%   INPUT:   
%       associationRate:    1D array of rates of a receptor merging with 
%                           cluster depending on cluster size
%       associationRate:    1D array of rates of a receptor dissociating
%                           from a cluster depending on cluster size
%       numReceptors:       total number of receptors
%       simTime:            total simulation time
%       timeStep:           simulation time step
%
%   OUTPUT:
%      clusterState:    a structure with the following fields
%         .members:     a logical 2D array indicating the receptors that
%                       make up the cluster. Rows correspond to receptors.
%                       The first column is the initial cluster state where
%                       all of the receptors are free (all 0s).
%                       Subsequent columns will only show dissociation and
%                       aggregation events.  The iteration values for these
%                       events are saved in .mergeIters and .dissIters.
%         .mergeIters:  a 2D array - column 1 gives the iteration
%                       points in .members where a merge occured and
%                       column 2 gives the absolute iteration point for the
%                       merge event in .members 
%         .dissIters:   a 2D array - column 1 gives the iteration
%                       points in .members where dissociation occured and
%                       column 2 gives the absolute iteration point for the
%                       dissociation event in .members 
%         .sizeHistory: a 2D array with 5 columns with each giving cluster
%                       information.  Column 1 is cluster size, columns 2
%                       and 3 are the iteration points where the cluster
%                       was formed and lost respectively, column 4 is the
%                       cluster's lifetime and column 5 indicates whether
%                       the cluster was lost via aggregation with another
%                       receptor (indicated by 1) or via dissociation into
%                       a smaller cluster and receptor (indicated by 0).
%   
%      merge2splitTime: a 2D array with 5 columns with each row
%                       representing a merge-to-split event.  Column 1
%                       gives the ID of the receptor dissociating, column 2
%                       gives the merge iteration point, column 3 gives the 
%                       splitting iteration point (current iter), column 4
%                       gives time the receptor spent as part of the
%                       cluster and column 5 is the current size of the
%                       cluster.
%
%   Robel Yirdaw, 08/16/13
%
    %The association probability
    associationProb = associationRate*timeStep;
    %The dissociaton probability
    dissociationProb = dissociationRate*timeStep;
    %Maximum number of simulation steps
    iterMax = floor(simTime/timeStep);

    %Set up return variables
    clusterState = struct('members',NaN(numReceptors,iterMax),...
        'mergeIter',zeros(iterMax,2),'dissIter',zeros(iterMax,2),...
        'sizeHistory',NaN(iterMax,5));
    %Record of merge-to-split events
    merge2splitTime = zeros(iterMax,5);
    %An array to track the aggregation iteration points for each receptor
    mergeTime = zeros(numReceptors,1);
    
    %Initialize    
    %Initially, all receptors are free - no cluster
    clusterState.members(:,1) = 0;        
    clusterSize = 1;
    %Or can be partially/all clustered    
    %clusterState.members(:,1) = 1;        
    %clusterState.members(1:floor(end/2),1) = 1;
    %clusterSize = sum(clusterState.members(:,1));
    
    
    
    %Initialize size history and index accordingly
    sizeHistIndx = 1;
    clusterState.sizeHistory(1,1) = clusterSize;
    clusterState.sizeHistory(1,2) = sizeHistIndx;
    %clusterState.sizeHistory(1,5) = 0;
    
    mergeCounter = 0;
    dissCounter = 0;
    m2sCounter = 0;
        
    %Begin
    for iterIndx=2:iterMax
                
        %Try clustering.
        if (rand(1) < associationProb(clusterSize))
            %Aggregation occurs if there are free receptors.  No free receptors
            %if cluster size = number of receptors.
            if (clusterSize < numReceptors)
                %Free receptors exist. 
                
                %The current cluster state
                currClustState = clusterState.members(:,mergeCounter+dissCounter+1);
                
                %If cluster does not exist then the first two receptors
                %will be clustered. Otherwise, randomly pick one receptor
                %to add to the existing cluster
                if (clusterSize == 1)
                    %All receptors free.  Pick the first two to cluster.
                    aggregRecepIndx = [1;2];
                else
                    %Get index of free receptors
                    freeRecep = find(currClustState == 0);
                    %Pick one of the free receptors to add to cluster.
                    aggregRecepIndx = freeRecep(randi(numel(freeRecep),1));
                end
                
                %Cluster selected receptor
                currClustState(aggregRecepIndx,1) = 1;

                %Update merge iter (time) value for aggregated receptor
                mergeTime(aggregRecepIndx) = iterIndx;

                %Update merge event counter
                mergeCounter = mergeCounter + 1;                 
                
                %Update cluster state and iteration values
                clusterState.members(:,mergeCounter+dissCounter+1) = currClustState;
                %The location of the current merge event in .members
                clusterState.mergeIter(mergeCounter,1) = mergeCounter+dissCounter+1;                
                %The iteration point of the current merge event
                clusterState.mergeIter(mergeCounter,2) = iterIndx;                

                %Update cluster size
                clusterSize = sum(currClustState);  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Update cluster size history
                %This is an end of current cluster...
                clusterState.sizeHistory(sizeHistIndx,3) = iterIndx;
                clusterState.sizeHistory(sizeHistIndx,4) = ...
                    iterIndx - clusterState.sizeHistory(sizeHistIndx,2);
                %cluster ends by association
                clusterState.sizeHistory(sizeHistIndx,5) = 1;
                %Increment counter
                sizeHistIndx = sizeHistIndx + 1;
                %...but also the beginning of a new cluster
                clusterState.sizeHistory(sizeHistIndx,1) = clusterSize;
                clusterState.sizeHistory(sizeHistIndx,2) = iterIndx;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
                
            end
            
        %Else try dissociating the cluster.
        elseif (rand(1) < dissociationProb(clusterSize))                    
            %Dissociation will be attempted only if there is a cluster.
            if (clusterSize > 1)   
            
                %Dissociating. Will remove a single receptor.

                %The current cluster state
                currClustState = clusterState.members(:,mergeCounter+dissCounter+1);
                %Index of members
                currClustMembers = find(currClustState ~= 0);    
                
                if (clusterSize == 2)
                    %Cluster is a dimer.  Will reset state value for both.
                    dissRecepIndx = currClustMembers;
                else                                       
                    %recept2dissociate = randi(clusterSize,1);
                    %Determine which receptor will dissociate - 1st, 2nd,...
                    %and get index of a receptor
                    dissRecepIndx = currClustMembers(randi(clusterSize,1));
                end
                
                %Remove receptor
                currClustState(dissRecepIndx) = 0;

                %Update dissociation counter
                dissCounter = dissCounter + 1;            
                
                %Update cluster state and iteration values
                clusterState.members(:,mergeCounter+dissCounter+1) = currClustState;
                %The location of the current dissociation event in .members
                clusterState.dissIter(dissCounter,1) = mergeCounter+dissCounter+1;
                %The iteration point of the current dissociation event
                clusterState.dissIter(dissCounter,2) = iterIndx;
                    
                %Update merge-to-split event counter
                m2sCounter = m2sCounter + 1;
                
                %Update merge-to-split values
                merge2splitTime(m2sCounter,:) = [dissRecepIndx(1) ...
                        mergeTime(dissRecepIndx(1)) iterIndx ...
                        (iterIndx - mergeTime(dissRecepIndx(1))) clusterSize];                
                
                if(numel(dissRecepIndx) == 2 && ...
                        (mergeTime(dissRecepIndx(1)) ~= mergeTime(dissRecepIndx(2))) )
                    %This dimer resulted from receptors merging with the 
                    %cluster at different time points - save merge-to-split
                    %time values for both and also increment counter                    
                    m2sCounter = m2sCounter + 1;
                    merge2splitTime(m2sCounter,:) = [dissRecepIndx(2) ...
                        mergeTime(dissRecepIndx(2)) iterIndx ...
                        (iterIndx - mergeTime(dissRecepIndx(2))) clusterSize];
                end
                        
                %Reset mergeTime for dissociated receptor
                mergeTime(dissRecepIndx) = 0;
                
                %Update cluster size
                clusterSize = sum(currClustState);                  
                %081613
                if (clusterSize == 0) 
                    clusterSize = 1;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Update cluster size history
                %This is an end of current cluster...
                clusterState.sizeHistory(sizeHistIndx,3) = iterIndx;
                clusterState.sizeHistory(sizeHistIndx,4) = ...
                    iterIndx - clusterState.sizeHistory(sizeHistIndx,2);
                %cluster ends by dissociation
                clusterState.sizeHistory(sizeHistIndx,5) = 0;
                %Increment counter
                sizeHistIndx = sizeHistIndx + 1;
                %but also the beginning of a new cluster
                clusterState.sizeHistory(sizeHistIndx,1) = clusterSize;
                clusterState.sizeHistory(sizeHistIndx,2) = iterIndx;                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            end
            
        end %If there is a cluster to dissociate               
        
    end %for every iteration step
    
    %Remove empty array rows and columns from return variables
    merge2splitTime(m2sCounter+1:end,:) = [];
    clusterState.members(:,mergeCounter+dissCounter+2:end) = [];
    clusterState.mergeIter(mergeCounter+1:end,:) = [];
    clusterState.dissIter(dissCounter+1:end,:) = [];
    clusterState.sizeHistory(sizeHistIndx+1:end,:) = [];

    
    
    