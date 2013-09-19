%   Script to run DissRateTest.
%   Modified on 07/30/13 to calculate cluster lifetimes. The dissociation
%   rate can be calculated as the inverse of the mean cluster lifetimes.
%   RBY, 07/16/13
    
    
    dissociationRate = 0.25;
    numReceptors = 100;
    simTime = 1000;
    timeStep = 0.01;

    %Total number of simulations
    simNum = 10;
    %dissEvents returned for each simulation
    dissEventsAll = zeros(numReceptors-1,simNum);
    %Cluster lifetimes
    clusterLTAll = zeros(numReceptors-1,simNum);
    
    for simIter=1:simNum
        [~,dissEvents] = DissRateTest(dissociationRate,numReceptors,simTime,timeStep);
        %Saving only the time each receptor spent in the cluster
        dissEventsAll(:,simIter) = dissEvents(:,3);
        %Calculate the lifetime of each cluster size. Descending order.
        %Since initially clustered, the lifetime for the first cluster,
        %largest cluster, is the number of iterations since simulation
        %started. For the rest, it is the difference between consecutive
        %iterations.
        clusterLTAll(1,simIter) = dissEvents(1,3).*timeStep;
        clusterLTAll(2:numReceptors-1,simIter) = diff(dissEvents(:,3)).*timeStep;        
    end
    
    %The average time spent in the cluster, by descending cluster size
    avgTimeByClust = mean(dissEventsAll(:,1:simNum),2).*timeStep;
    %The overall average time spent in the cluster
    avgTime = mean(dissEventsAll(1:simNum*(numReceptors-1)))*timeStep;
    %Average lifetime of each cluster size
    avgClusterLT = mean(clusterLTAll,2);
    
   