function [clusterState,dissEvents] = DissRateTestII(dissociationRate,numReceptors,simTime,timeStep)
%DISSRATETEST simulates the dissociation of a single cluster with an
%   initial size equal to the number of receptors. The dissociation occurs
%   along two paths, i.e. a parallel reaction. 
%   The simulation starts with all receptors forming a single cluster.  The
%   simulation will run until all receptors have dissociated or the maximum
%   iteration has been reached.  A dissociation probability for the cluster
%   will be calculated from the user-defined dissociation rate.  If the
%   cluster is dissociating, a receptor will be randomly selected and
%   removed. Therefore, receptors are removed one at a time. At every
%   iteration, a cluster can dissociate via one of the two paths.
%
%   INPUT:      
%       dissociationRate:   2-by-1 array of dissociation rates
%       numReceptors:       total number of receptors
%       simTime:            total simulation time
%       timeStep:           simulation time step
%
%   OUTPUT:
%       clusterState:   a structure with the following fields
%           .members:   a logical 2D array indicating the receptors that
%                       make up the cluster. Rows correspond to receptors.
%                       The first column is the initial cluster state where
%                       all of the receptors are clustered (all 1s).
%                       Subsequent columns will only show dissociation
%                       events, each time a receptor is removed.
%
%           .iters:     iteration points corresponding to the columns in
%                      .members.
%
%       dissEvents:     a 2D array with each row corresponding to a
%                       dissociaiton event.  Column 1 gives the size of the
%                       dissociating cluster, column 2 gives the iteration
%                       point where the cluster was formed, column 3 gives
%                       the iteration point where the cluster ends
%                       (a receptor dissociates), column 4 is the lifetime
%                       of the cluster and column 5 indicates how the 
%                       cluster dissociated, the "reaction" type, where 1
%                       is type 1 via dissociation rate 1 and 2 is type 2
%                       via dissociation rate 2.
%                       
%   Robel Yirdaw, 08/08/13
%
    
    %Set up return variables
    clusterState = struct('members',zeros(numReceptors,numReceptors),...
        'iter',zeros(numReceptors,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dissEvents = zeros(numReceptors-1,5);       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The dissociaton probability
    dissociationProb = dissociationRate*timeStep;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Maximum number of simulation steps
    iterMax = floor(simTime/timeStep);
    
    %Initialize
    %All receptors form a single cluster.
    clusterState.members(:,1) = 1;
    clusterSize = sum(clusterState.members(:,1));
    
    dissEvents(1,1) = numReceptors;
    dissEvents(1,2) = 1;
    
    %Dissociation coutner
    dissCounter = 1;
    %Sim counter
    iterIndx = 2;
    
    %Begin
    while ((clusterSize > 1) && (iterIndx <= iterMax))
        %Try dissociating the cluster. Two possible paths for dissociation.        
        if (rand(1) < dissociationProb(1))
            eventID = 1;
        elseif (rand(1) < dissociationProb(2))
            eventID = 2;
        else
            eventID = 0;
        end
        
        if (eventID ~= 0)
            %Dissociating. Will remove a single receptor.
            
            %First determine which receptor will dissociate - 1st, 2nd,...
            recept2dissociate = randi(clusterSize,1); 
            
            %The current cluster state
            currClustState = clusterState.members(:,dissCounter);
            %Index of members
            currClustMembers = find(currClustState ~= 0);    
            %Get index of receptor to be removed
            dissRecepIndx = currClustMembers(recept2dissociate);
            %Remove receptor
            currClustState(dissRecepIndx) = 0;            
            
            %Update cluster state 
            clusterState.members(:,sum(dissCounter)+1) = currClustState;
            clusterState.iter(sum(dissCounter)+1,1) = iterIndx;
            
            %%%END OF A CLUSTER%%%
            %The current cluster ends
            dissEvents(dissCounter,3) = iterIndx;
            %Lifetime of the cluster
            dissEvents(dissCounter,4) = iterIndx - dissEvents(dissCounter,2);            
            %Type of dissociation event that ended the cluster
            dissEvents(dissCounter,5) = eventID;            

            %Update cluster size
            clusterSize = sum(currClustState);  
            %Update dissociation event counter
            dissCounter = dissCounter + 1;               
            
            %%%BEGINNING OF A NEW CLUSTER%%%
            %The current cluster size
            dissEvents(dissCounter,1) = clusterSize;
            %The beginning of a new cluster
            dissEvents(dissCounter,2) = iterIndx;
            
             
        end
        
        %Increment sim step.
        iterIndx = iterIndx + 1;
        
    end %while cluster exists and not end of sim
    
    