%   Script to run DissRateTestWithAggreg.
%   Modified 07/30/13 to calculate dissociation rates from cluster
%   lifetime vaues.
%   RBY, 07/24/13
    

    aggregationProb = 0.1;
    dissociationRate = 5;
    numReceptors = 5;
    simTime = 1000;
    timeStep = 0.01;

    %Total number of simulations
    simNum = 10;
    %Will save merge-to-split time values and corresponding dissociation 
    %rates by cluster for each simulation
    m2sTimeByClustSizeAll = zeros(numReceptors,simNum);
    calcDissRateByClust = zeros(numReceptors,simNum);
    %Will save average merge-to-split time values for each simulation
    avgMerge2splitTime = zeros(simNum,1);    
    %Will save lifetime of each cluster size for each simulation
    clusterLifeTimeAll = zeros(numReceptors,simNum);
    
    %rng(100);
    %allRN = zeros(simNum,1);
    
    for simIter=1:simNum
        
        %rN = randi(1000,1);
        %allRN(simIter,1) = rN;
        %rN = allRN(simIter);
        %rng(rN);
        
        [clusterState,m2sTime] = DissRateTestWithAggreg(aggregationProb,dissociationRate,numReceptors,simTime,timeStep);
        
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
        %Acculmulate merge-to-split times for repeated simulations
        m2sTimeByClustSizeAll(1:maxClusterSize,simIter) = m2sTimesByClust(:,3);    
        %Calculate the corresponding dissociation rate
        calcDissRateByClust(:,simIter) = 1./m2sTimeByClustSizeAll(:,simIter);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        %Cluster lifetime by cluster size
        %The last cluster will not be used since its end is unknown
        sizeHist = clusterState.sizeHistory(1:end-1,:);
        maxClusterSize = max(sizeHist(:,1));
        clusterLifeTime = zeros(maxClusterSize,10);
        %Set up first column, the cluster sizes - first one is size 0.
        clusterLifeTime(2:end,1) = 2:maxClusterSize;
        %Process all cluster sizes - starting with size 2. In the following
        %the lifetime values are calculated three different ways.  Each has
        %three columns: Col.1 = the number of times the cluster existed,
        %Col.2 = total time cluster existed, Col.3 = cluster lifetime.
        for clustIndx=2:maxClusterSize              
            %Lifetime of clusters ending in aggregation and dissociation
            currClustIndx = (sizeHist(:,1) == clustIndx);
            clusterLifeTime(clustIndx,1) = sum(currClustIndx);
            clusterLifeTime(clustIndx,2) = sum(sizeHist(currClustIndx,4).*timeStep);            
            clusterLifeTime(clustIndx,3) = clusterLifeTime(clustIndx,2)/clusterLifeTime(clustIndx,1);    
            
            %Lifetime of clusters ending in dissociation only
            currClustIndx_Diss = (sizeHist(:,1) == clustIndx) & (sizeHist(:,5) == 0);
            clusterLifeTime(clustIndx,4) = sum(currClustIndx_Diss);
            clusterLifeTime(clustIndx,5) = sum(sizeHist(currClustIndx_Diss,4).*timeStep);            
            clusterLifeTime(clustIndx,6) = clusterLifeTime(clustIndx,5)/clusterLifeTime(clustIndx,4);   
            
            %Lifetime of clusters ending in aggregation only
            currClustIndx_Merge = (sizeHist(:,1) == clustIndx) & (sizeHist(:,5) == 1);
            clusterLifeTime(clustIndx,7) = sum(currClustIndx_Merge);
            clusterLifeTime(clustIndx,8) = sum(sizeHist(currClustIndx_Merge,4).*timeStep);            
            clusterLifeTime(clustIndx,9) = clusterLifeTime(clustIndx,8)/clusterLifeTime(clustIndx,7);              

            %Lifetime values are modified depending on the number of
            %reaction paths available for a current state - i.e.
            %association and/or dissociation, for a cluster of given size. 
            %Both are possible for all cluster sizes except the largest.            
            %See summary from 080113.
            if (clustIndx < maxClusterSize)                
                %The current cluster size can undergo association and
                %dissociation
                clusterLifeTime(clustIndx,10) = 2*clusterLifeTime(clustIndx,3);  
            else
                %The current cluster size can undergo dissociation only
                clusterLifeTime(clustIndx,10) = clusterLifeTime(clustIndx,3);  
            end

        end
        %Acculmulate lifetimes by cluster size, for repeated simulations
        clusterLifeTimeAll(1:maxClusterSize,simIter) = clusterLifeTime(:,6);
        
    end
    
    %Calculate dissociation rate from merge-to-split time 
    calcDissRate = 1./avgMerge2splitTime;           
    
    %Calculate dissociation rate from lifetime values
    calcDissRateFromLT = 1./mean(clusterLifeTimeAll(2:end,:),2);
    
    