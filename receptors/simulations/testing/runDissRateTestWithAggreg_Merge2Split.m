%   Script to run DissRateTestWithAggreg_Merge2Split.
%   RBY, 07/24/13
    

    aggregationProb = 0.5;
    dissociationRate = 5;
    numReceptors = 10;
    simTime = 1000;
    timeStep = 0.01;

    %Total number of simulations
    simNum = 1;
    %Will save merge-to-split time values and corresponding dissociation 
    %rates by cluster for each simulation
    m2sTimeByClustSizeAll = zeros(numReceptors,simNum);
    calcDissRateByClust = zeros(numReceptors,simNum);
    %Will save average merge-to-split time values for each simulation
    avgMerge2splitTime = zeros(simNum,1);    
    
    %rng(100);
    
    for simIter=1:simNum
        
        [clusterState,m2sTime] = DissRateTestWithAggreg_Merge2Split(aggregationProb,dissociationRate,numReceptors,simTime,timeStep);
        
        %Calculate average merge-to-split time
        avgMerge2splitTime(simIter,1) = mean(m2sTime(:,4).*timeStep); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Merge-to-split times by cluster size
        maxClusterSize = max(m2sTime(:,5));
        %Will save the cluster size in column 1, number of events in column
        %2 and average merge-to-split time in column 3, for each cluster
        %size (rows)
        m2sTimesByClust = zeros(maxClusterSize,3);
        %Set up first column, the cluster sizes
        m2sTimesByClust(:,1) = 1:maxClusterSize;
        %Process all cluster sizes
        for clustIndx=2:maxClusterSize
            tempM2STimeIndx = m2sTime(:,5) == clustIndx;
            m2sTimesByClust(clustIndx,2) = sum(tempM2STimeIndx);
            m2sTimesByClust(clustIndx,3) = mean(m2sTime(tempM2STimeIndx,4).*timeStep);
        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        %Acculmulate merge-to-split times for repeated simulations
        m2sTimeByClustSizeAll(1:maxClusterSize,simIter) = m2sTimesByClust(:,3);    
        %Calculate the corresponding dissociation rate
        calcDissRateByClust(:,simIter) = 1./m2sTimeByClustSizeAll(:,simIter);
                
    end
    
    %Calculate dissociation rate
    calcDissRate = 1./avgMerge2splitTime;
        
    

    