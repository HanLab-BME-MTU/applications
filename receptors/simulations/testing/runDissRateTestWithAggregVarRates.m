%   Script to run DissRateTestWithAggregVarRates and calculate association and
%   dissociation rates by first determining the probability of an event and
%   then using the lifetime corresponding to that event to obtain the rate
%   value, i.e., the rate of event r is given by 
%       kr = pr*sum(ki), where sum(ki) is inverse of lifetime
%
%   Note: association and dissociation rates for individual cluster sizes
%   can be provided here.
%
%   RBY, 08/19/13
    
    %Total number of receptors
    numReceptors = 5;
    %Rates are vectors with each row giving a value for cluster size = row
    %aggregationRate(1:numReceptors,1) = 8;
    %dissociationRate(1:numReceptors,1) = 4;
    aggregationRate = [3;3;0;0;0];
    dissociationRate = [0;2;4;0;0];
    
    %Simultation time and time step
    simTime = 100;
    timeStep = 0.001;

    %Total number of simulations
    simNum = 50;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %08/12/13    
    calcDissocRateAll = zeros(simNum,1);
    calcAssocRateAll = zeros(simNum,1);
    eventTableAll = cell(simNum,1);
    
    calcDissocRateByClust = zeros(numReceptors,simNum);
    calcAssocRateByClust = zeros(numReceptors,simNum);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
    %allRN = zeros(simNum,1);
    
    for simIter=1:simNum
        
        %rN = randi(1000,1);
        %allRN(simIter,1) = rN;
        rN = allRN(simIter);
        rng(rN);
        
        [clusterState,~] = DissRateTestWithAggregVarRates(aggregationRate,dissociationRate,numReceptors,simTime,timeStep);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %08/12/13
        %Overall rates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sizeHist = clusterState.sizeHistory(1:end-1,:);
        
        lifeTimeAny = mean(sizeHist(:,4))*timeStep;
        
        numAssocEvents = sum(sizeHist(:,5) == 1);       
        numDissocEvents = sum(sizeHist(:,5) == 0);

        probAssocAll = numAssocEvents/(numAssocEvents+numDissocEvents); 
        probDissocAll = numDissocEvents/(numAssocEvents+numDissocEvents);
        
        calcAssocRateAll(simIter,:) = probAssocAll./lifeTimeAny;
        calcDissocRateAll(simIter,:) = probDissocAll./lifeTimeAny;        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        %Rates by cluster size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        maxClusterSize = max(sizeHist(:,1));
        eventTable = zeros(maxClusterSize,8);       

        %Cluster size = 1
        currClustIndx_All = (sizeHist(:,1) == 1);
        %Total # of events
        eventTable(1,1) = sum(currClustIndx_All);
        %Overall lifetime
        eventTable(1,2) = mean(sizeHist(currClustIndx_All,4))*timeStep;
        %Association events
        currClustIndx_Assoc = ((sizeHist(:,1) == 1) & (sizeHist(:,5) == 1));
        %Total # of association events
        eventTable(1,6) = sum(currClustIndx_Assoc);            
        %Prob. of association
        eventTable(1,7) = eventTable(1,6)/eventTable(1,1);
        %Association rate
        eventTable(1,8) = eventTable(1,7)/eventTable(1,2);
        
        %Cluster sizes 2 and higher
        for clustIndx=2:maxClusterSize
            
            currClustIndx_All = (sizeHist(:,1) == clustIndx);
            %Total # of events
            eventTable(clustIndx,1) = sum(currClustIndx_All);
            %Overall lifetime
            eventTable(clustIndx,2) = mean(sizeHist(currClustIndx_All,4))*timeStep;
            
            %Dissociation events
            currClustIndx_Dissoc = ((sizeHist(:,1) == clustIndx) & (sizeHist(:,5) == 0));
            %Total # of dissociation events
            eventTable(clustIndx,3) = sum(currClustIndx_Dissoc);            
            %Prob. of dissociation
            eventTable(clustIndx,4) = eventTable(clustIndx,3)/eventTable(clustIndx,1);
            %Dissociation rate
            eventTable(clustIndx,5) = eventTable(clustIndx,4)/eventTable(clustIndx,2);

            %Association events
            currClustIndx_Assoc = ((sizeHist(:,1) == clustIndx) & (sizeHist(:,5) == 1));
            %Total # of association events
            eventTable(clustIndx,6) = sum(currClustIndx_Assoc);            
            %Prob. of association
            eventTable(clustIndx,7) = eventTable(clustIndx,6)/eventTable(clustIndx,1);
            %Association rate
            eventTable(clustIndx,8) = eventTable(clustIndx,7)/eventTable(clustIndx,2);
        end
        
        %Save rates for repeated simulations
        calcDissocRateByClust(1:maxClusterSize,simIter) = eventTable(:,5);
        calcAssocRateByClust(1:maxClusterSize,simIter) = eventTable(:,8);
        %Save eventTable for repeated simulations
        eventTableAll{simIter,1} = eventTable;
        
        clear eventTable

    end %Num sims
    
    
    