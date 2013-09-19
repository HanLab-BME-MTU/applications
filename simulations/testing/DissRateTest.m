function [clusterState,dissEvents] = DissRateTest(dissociationRate,numReceptors,simTime,timeStep)
%DISSRATETEST simulates the dissociation of a single cluster with an
%   initial size equal to the number of receptors.
%   The simulation starts with all receptors forming a single cluster.  The
%   simulation will run until all receptors have dissociated or the maximum
%   iteration has been reached.  A dissociation probability for the cluster
%   will be calculated from the user-defined dissociation rate.  If the
%   cluster is dissociating, a receptor will be randomly selected and
%   removed. Therefore, receptors are removed one at a time.
%
%   INPUT:      
%       dissociationRate:   rate of dissociating a receptor from cluster
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
%                       dissociaiton event.  Column 1 gives the ID of the
%                       receptor currently dissociating, column 2 gives the
%                       current cluster size and column 3 is the time the
%                       receptor spent as part of the clsuster (elapsed
%                       time since cluster is formed at beginning of 
%                       simulation).
%
%   Robel Yirdaw, 07/16/13
%
    
    %Set up return variables
    clusterState = struct('members',zeros(numReceptors,numReceptors),...
        'iter',zeros(numReceptors,1));
    dissEvents = zeros(numReceptors-1,3);    
    %The dissociaton probability
    dissociationProb = dissociationRate*timeStep;
    %Maximum number of simulation steps
    iterMax = floor(simTime/timeStep);
    
    %Initialize
    %All receptors form a single cluster.
    clusterState.members(:,1) = 1;
    clusterSize = sum(clusterState.members(:,1));
    
    %Dissociation coutner
    dissCounter = 1;
    %Sim counter
    iterIndx = 1;
    
    %Begin
    while ((clusterSize > 1) && (iterIndx <= iterMax))
        %Try dissociating the cluster.
        if (rand(1) < dissociationProb)
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
            clusterState.members(:,dissCounter+1) = currClustState;
            clusterState.iter(dissCounter+1,1) = iterIndx;
            
            %The receptor dissociating
            dissEvents(dissCounter,1) = dissRecepIndx;
            %The current cluster size
            dissEvents(dissCounter,2) = clusterSize;
            %Time the receptor spent in the cluster
            dissEvents(dissCounter,3) = iterIndx;

            %Update cluster size
            clusterSize = sum(currClustState);                 
            %Update dissociation event counter
            dissCounter = dissCounter + 1;                 
        end
        
        %Increment sim step.
        iterIndx = iterIndx + 1;
        
    end %while cluster exists and not end of sim
    
    