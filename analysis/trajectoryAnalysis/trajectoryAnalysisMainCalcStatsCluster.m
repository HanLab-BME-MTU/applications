function [bestClusters,numClusters,indClusters] = trajectoryAnalysisMainCalcStatsCluster(dataListG,gosIdx,constants)
%TAMCS-CLUSTER clusters the speeds in dataListG indexed by gosIdx according to the settings in constants
%
% OUT: clusterStruct.bestClusters: #of clusters, means, std, weight
%                   .numClusters : #of clusters, #of sig Custers every
%                   iteration ([1stIter(clustertryX2),2ndIter,3rdIter etc])
%                   
%
% 3/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% Defaults
%=============
CLUSTERTOL = 1e-7;
clustermax = constants.CLUSTERMAX;
clustermin = constants.CLUSTERMIN;
%============

%============
% init
%============

bestClusters = [];         
numClusters  = [];          
indClusters(1:(constants.CLUSTERIND-1)*constants.CLUSTERMAX+1)  = deal([]);

%============


%========================
% prepare clusterall-loop
%========================

% because it is not possible (yet) to weigh measurements in
% mixtures4, we repeat speeds with longer times (numIntervals = 1,2,3 compared to shortest)
indSpeeds       = 60*dataListG(gosIdx,4);
numIntervals    = dataListG(gosIdx,7);
numIntervals    = round(numIntervals/min(numIntervals));
indSpeeds       = repeatEntries(indSpeeds,numIntervals);

%========================

%========================
% clusterall-loop
%=======================

% loop until we found n clusters with a minimum weight of 5%
% each.
done = 0;
clusterNumber = 0;
muInit = [];
ppInit = [];

while ~done
    
    % init array, because you never know how many clusters you'll get
    bestClusters = zeros(constants.CLUSTERTRY,3*clustermax+2); 
    sumWeight = zeros(constants.CLUSTERTRY,1);
    
    % repeat CLUSTERTRY times to avoid the (hopefully) lone outlier
    for iTry = 1:constants.CLUSTERTRY
        
        % initialize new cluster with 'reasonable' values
        if ~isempty(muInit)
            mu0 = muInit(1,:) + randn(1,clusterNumber) .* muInit(2,:);
            switch sign(indSpeeds(1))
                case 1 % growth
                    mu0 = max(mu0, zeros(1,clusterNumber));
                    mu0 = min(mu0, repmat(max(indSpeeds),[1,clusterNumber]));
                case -1 % shrinkage
                    mu0 = min(mu0, zeros(1,clusterNumber));
                    mu0 = max(mu0, repmat(min(indSpeeds),[1,clusterNumber]));
            end
            % init weights. Make sure their sum is equal to 1 and none is 0
            pp0 = (ppInit(1,:) + randn(1,clusterNumber) .* ppInit(2,:));
            pp0 = max(pp0,repmat(0.01,[1,clusterNumber]));
            pp0 = pp0/sum(pp0);
        else
            mu0 = [];
            pp0 = [];
        end
        
        % #of clusters, relative weights, positions, covariances = 
        %       m4(speeds,min#ofClusters,max#ofClusters,regularize,threshold,some option,pp0,mu0,verbose)
        [bestk,bestpp,bestmu,bestcov] = mixtures4(indSpeeds',...
            clustermin,clustermax,0,CLUSTERTOL,1,pp0,mu0,0);
        % and store: number of clusters, number of significant clusters,
        % means, variance, weight
        [means,muIdx] = sort(bestmu);
        weight = bestpp(muIdx);
        sumWeight(iTry) = sum(weight);
        numSigClusters = nnz(weight > constants.CLUSTERMINWEIGHT);
        variance = sqrt(bestcov(:,:,muIdx));
        bestClusters(iTry,1:bestk*3+2) = [bestk, numSigClusters, means, squeeze(variance)', weight]; 
    end % for iTry = 1:x
    
    % mixtures4 seems to hang itself sometimes and return sum of weights > 1 
    % particularly if it is initiated with weights > 1. Check anyway!
    bigWeightIdx= find(abs(sumWeight-1)>4*eps);
    if ~isempty(bigWeightIdx)
        warning('TRAJECTORYANALYSIS:clusterweight','weight of all clusters is > 1!')
        bestClusters(bigWeightIdx,:)
        bestClusters(bigWeightIdx,:)=[];
    end
    
    
    % now look for the correct number of significant clusters
    bestNumSigClusters = round(mean(bestClusters(:,2)));
    overallBestK       = round(mean(bestClusters(:,1)));
    
    
    
    if bestNumSigClusters == clusterNumber % the first time, this will not be satisfied
        
        % take the clusters whose number of clusters is equal to the number
        % of best significant clusters, and where all clusters are
        % significant
        bestIdx = find(all(bestClusters(:,1:2)==bestNumSigClusters,2));
        
        % that's it
        done = 1;
    else
        % we force the software to return the right number of clusters
        
        % clusterNumber = findClosest(bestClusters(:,2),bestNumSigClusters,'max',1);
        sigVec = bestClusters(:,2);
        sigVec = unique(sigVec);
        difference = abs(sigVec-bestNumSigClusters);
        [dummy,closeIdx] = min(difference);
        clusterNumber = max(sigVec(closeIdx));
        % for the moment: just try to cap the max number, allow all other
        % mins - I do not know how the algorithm behaves.
        clustermax = clusterNumber;
        
        % get info for starting points for clusters
        bestIdx = find(all(bestClusters(:,1:2)==clusterNumber,2));
        if ~isempty(bestIdx)
            % mu
            muInit = [];
            muInit(1,:) = mean(bestClusters(bestIdx, 3:3+clusterNumber-1),1);
            if length(bestIdx) > 1
                muInit(2,:) = std(bestClusters(bestIdx, 3:3+clusterNumber-1),0,1);
            else 
                muInit(2,:) = bestClusters(bestIdx, 3+clusterNumber:3+2*clusterNumber-1);
            end
            % weight (pp)
            ppInit = [];
            ppInit(1,:) = mean(bestClusters(bestIdx, 3+2*clusterNumber:3+3*clusterNumber-1),1);
            if length(bestIdx) > 1
                ppInit(2,:) =  std(bestClusters(bestIdx, 3+2*clusterNumber:3+3*clusterNumber-1),0,1);
            else
                ppInit(2,:) =  bestClusters(bestIdx, 3+2*clusterNumber:3+3*clusterNumber-1)*0.1;
            end
        else
            muInit = [];
            ppInit = [];
        end
    end
    
    % remember clusterNumbers. Put side by side, because number never
    % changes
    numClusters = [numClusters,bestClusters(:,1:2)];
    
end %======= while ~done =============

if isempty(bestIdx) % just for safety
    bestClusters = []; % avoid calling mean for nothing and creating warning
else
    % since the results are sorted, we take the correct means
    bestClusters = [bestNumSigClusters,mean(bestClusters(bestIdx,3:3*bestNumSigClusters+2),1)];
end



%========================


%===================================
% indCluster - just for completeness
%===================================

if constants.CLUSTERIND
    
    %find results for every cluster
    for kCluster = constants.CLUSTERMIN:constants.CLUSTERMAX
        % force EM to return data only for selected k of means
        [bestk,bestpp,bestmu,bestcov] = mixtures4(indSpeeds',...
            kCluster,kCluster,0,CLUSTERTOL,1,[],[],any(verbose==4));
        
        [means,muIdx] = sort(bestmu);
        weight = bestpp(muIdx);
        variance = sqrt(bestcov(:,:,muIdx));
        
        clusterStruct(kCluster).indClusters = [means',squeeze(variance),weight'];
    end
    
end

%========================