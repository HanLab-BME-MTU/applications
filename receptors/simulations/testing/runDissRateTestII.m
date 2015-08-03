%   Script to run DissRateTestII.
%   RBY, 08/08/13
    
    rng(100);
    
    %Dissociation occurs two ways
    dissociationRate = [2; 2];
    %Total number of receptors (starting cluster size)
    numReceptors = 5;
    %Total simulation time and time-step
    simTime = 1000;
    timeStep = 0.01;

    %Total number of simulations
    simNum = 1000;

    %Cluster lifetimes will be saved
    clusterLTAll = zeros(numReceptors-1,simNum);
    clusterLTAll1 = zeros(numReceptors-1,simNum);
    clusterLTAll2 = zeros(numReceptors-1,simNum);
    
    for simIter=1:simNum
        [clusterState,dissEvents] = DissRateTestII(dissociationRate,numReceptors,simTime,timeStep);
        
        %Calculate the lifetime of each cluster size, in descending order
        %of cluster size.  Three types of values: overall, rate 1 (event 1)
        %only and rate 2 (event 2) only.
        
        %Overall
        clusterLTAll(1:numReceptors-1,simIter) = dissEvents(1:end-1,4).*timeStep;     
        %Event 1
        event1Indx = dissEvents(1:end-1,5) == 1;
        clusterLTAll1(event1Indx,simIter) = dissEvents(event1Indx,4).*timeStep;                       
        %Event 2
        event2Indx = dissEvents(1:end-1,5) == 2;
        clusterLTAll2(event2Indx,simIter) = dissEvents(event2Indx,4).*timeStep;                      

    end
    
    %Calculate mean lifetime values
    avgClusterLT = mean(clusterLTAll,2);
    avgClusterLT1 = zeros(numReceptors-1,1);
    avgClusterLT2 = zeros(numReceptors-1,1);
    avgClusterLT12 = zeros(numReceptors-1,1);

    for indx=1:(numReceptors-1)
        avgClusterLT1(indx,1) = mean(clusterLTAll1(indx,clusterLTAll1(indx,:)~=0),2);
        avgClusterLT2(indx,1) = mean(clusterLTAll2(indx,clusterLTAll2(indx,:)~=0),2);
    end    
    %Mean lifetime taking into account that there are two events possible.
    %The factor of 2 is for the special case of k1 = k2.
    if (dissociationRate(1) == dissociationRate(2))
        avgClusterLT12 = mean(2*clusterLTAll,2);
    end
    
    %Dissociation rates for the above
    [1./avgClusterLT 1./avgClusterLT1 1./avgClusterLT2 1./avgClusterLT12]
    
    
   