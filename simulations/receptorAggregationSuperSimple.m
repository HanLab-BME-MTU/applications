function time2aggreg = receptorAggregationSuperSimple(modelParam,simParam)
%receptorAggregationSuperSimple simulates the aggregation of receptors diffusing in a preset area; code stops when first aggregation happens
%
%SYNOPSIS time2aggreg = receptorAggregationSuperSimple(modelParam,simParam)
%
%INPUT  modelParam: Structure with the fields:
%           diffCoef        : Diffusion coefficient (microns^2/s).
%           confDim         : Confinement dimension (microns). Either one
%                             value, used for all probDims, or a row
%                             vector with a value for each probDim.
%                             NaN means free diffusion.
%           receptorDensity : Receptor density (#/microns^probDim).
%           aggregationDist : Distance between 2 receptors to consider them
%                             bumping into each other, thus potentially
%                             aggregating (microns).
%       simParam: Structure with the fields:
%           probDim         : System dimensionality (1, 2 or 3). Default: 2.
%           areaSideLen     : Area side length (microns). Either one
%                             value, used for all probDims, or a row
%                             vector with a value for each probDim.
%                             Default: 1 in all probDims.
%           timeStep        : Simulation time step (s).
%                             Default: 0.01/max(diffCoef,dissociationRate).
%           maxTime         : Maximum simulation time to stop simulaiton if
%                             nothing aggregates (s).
%                             Default: 100 * timeStep.
%           randNumGenSeeds : Row vector storing the seeds for the "rand"
%                             and "randn" random number generators.
%                             Default: [100 100].
%           numRep          : Number of times to repeat the simulation to
%                             get a distribution of time2aggreg.
%                             Default: 100.
%                 Whole structure optional. Individual fields also optional.
%
%OUTPUT time2aggreg: Vector of length numRep with distribution of times to
%                    aggregate.

%Khuloud Jaqaman, November 2013

%% Output

time2aggreg = [];
errFlag = 0;

%% Input

%check if correct number of arguments were used when function was called
if nargin < 1
    disp('--receptorAggregationSuperSimple: Too few input arguments');
    errFlag  = 1;
    return
end

%extract model parameters from modelParam

%diffusion coefficient
if isfield(modelParam,'diffCoef')
    diffCoef = modelParam.diffCoef;
else
    disp('--receptorAggregationSuperSimple: Please supply diffusion coefficient');
    errFlag = 1;
end

%confinement dimension
if isfield(modelParam,'confDim')
    confDim = modelParam.confDim;
else
    disp('--receptorAggregationSuperSimple: Please supply confinement dimension');
    errFlag = 1;
end

%receptor density
if isfield(modelParam,'receptorDensity')
    receptorDensity = modelParam.receptorDensity;
else
    disp('--receptorAggregationSuperSimple: Please supply receptor density');
    errFlag = 1;
end

%aggregation probability
aggregationProb = 1;

%aggregation distance
if isfield(modelParam,'aggregationDist')
    aggregationDist = modelParam.aggregationDist;
else
    disp('--receptorAggregationSuperSimple: Please supply aggregation distance');
    errFlag = 1;
end

%exit if there are missing model parameters
if errFlag == 1
    disp('--receptorAggregationSuperSimple: Please supply missing variables!');
    return
end

%extract simulation parameters from simParam

%if simParam wasn't supplied at all
if nargin < 2 || isempty(simParam)
    
    probDim = 2;
    areaSideLen = ones(1,probDim);
    timeStep = 0.01 / max(diffCoef,dissociationRate);
    maxTime = 100 * timeStep;
    randNumGenSeeds = [100 100];
    numRep = 100;
    
else
    
    %system probDimality
    if isfield(simParam,'probDim')
        probDim = simParam.probDim;
    else
        probDim = 2;
    end
    
    %observation side length
    if isfield(simParam,'areaSideLen')
        areaSideLen = simParam.areaSideLen;
        if length(areaSideLen) == 1
            areaSideLen = areaSideLen * ones(1,probDim);
        end
    else
        areaSideLen = ones(1,probDim);
    end
    
    %time step
    if isfield(simParam,'timeStep')
        timeStep = simParam.timeStep;
    else
        timeStep = 0.01 / max(diffCoef,dissociationRate);
    end
    
    %maximum simulation time
    if isfield(simParam,'maxTime')
        maxTime = simParam.maxTime;
    else
        maxTime = 100 * timeStep;
    end
    
    %random number generator seeds
    if isfield(simParam,'randNumGenSeeds')
        randNumGenSeeds = simParam.randNumGenSeeds;
    else
        randNumGenSeeds = [100 100];
    end
    
    %number of repetitions
    if isfield(simParam,'numRep')
        numRep = simParam.numRep;
    else
        numRep = 100;
    end
    
end

%expand confinement dimension if needed
if length(confDim) == 1
    confDim = confDim * ones(1,probDim);
end

%determine number of iterations to perform
numIterSim = ceil( maxTime / timeStep ) + 1;
numIterations = numIterSim;

%assuming that the displacement taken in a timeStep is normally distributed
%with mean zero, get standard deviation of step distribution from diffusion
%coefficient
stepStd = sqrt( 2 * diffCoef * timeStep );

%in case of confined diffusion, determine how many times simulation must be
%repeated in order to explore same number of receptors as in a free
%diffusion case
obsRegionSize = prod(areaSideLen);
if any(~isnan(confDim))
    confArea = prod(confDim);
    numRepConf = round(obsRegionSize/confArea);
    areaSideLen = confDim;
    obsRegionSize = confArea;
else
    numRepConf = 1;
end

%adjust aggregationDist to account for the finite simulation time step and
%the expected receptor displacement in that time step
aggregationDist = max(aggregationDist,2*sqrt(probDim)*stepStd);

%initialize random number generators
% rand('twister',randNumGenSeeds(1));
% randn('state',randNumGenSeeds(2));
rng(randNumGenSeeds(1),'twister')

%initialize output vector
time2aggreg = NaN(numRep,1);

progressText(0,'Simulation');

for iRep = 1 : numRep
    
    %initialize intermediate output vector
    time2aggregTmp = NaN(numRepConf,1);
    
    for iRepConf = 1 : numRepConf
        
        %% receptor initial positions and clustering
        
        %calculate number of receptors
        numReceptors = round(obsRegionSize * receptorDensity);
        
        %initialize receptor positions
        initPositions = rand(numReceptors,probDim) .* repmat(areaSideLen,numReceptors,1);
        
        %based on these positions, cluster receptors (and modify their positions)
        aggregationProbVec = aggregationProb*ones(numReceptors,1);
        [cluster2receptor,receptor2cluster,clusterSize,initPositions] = ...
            receptorAggregationAlg(initPositions,aggregationDist,aggregationProbVec);
        [numClusters,maxClustSize] = size(cluster2receptor);
        
        if maxClustSize > 1
            
            time2aggregTmp(iRepConf) = 0;
            
        else
            
            %% Trajectory generation
            
            %reserve memory for output vectors
            receptorTraj = zeros(numReceptors,probDim,numIterations);
            recept2clustAssign = zeros(numReceptors,numIterations);
            clust2receptAssign = zeros(numReceptors,maxClustSize,numIterations);
            
            %store initial information
            receptorTraj(:,:,1) = initPositions;
            recept2clustAssign(:,1) = receptor2cluster;
            clust2receptAssign(1:numClusters,1:maxClustSize,1) = cluster2receptor;
            
            %iterate in time
            iIter = 1;
            while maxClustSize == 1 && iIter <= numIterations
                
                %increase iteration index by 1
                iIter = iIter + 1;
                
                %get receptor positions at previous time point
                positionsOld = receptorTraj(:,:,iIter-1);
                
                %generate receptor displacements
                receptorDisp = stepStd*randn(numReceptors,probDim);
                
                %calculate the new receptor positions at the current time point
                positionsNew = positionsOld + receptorDisp;
                
                %make sure that receptors stay inside the region of interest
                correctionBoundaryLow = min(positionsNew,0);
                positionsNew = positionsNew - 2 * correctionBoundaryLow;
                correctionBoundaryUp = max(positionsNew - repmat(areaSideLen,numReceptors,1),0);
                positionsNew = positionsNew - 2 * correctionBoundaryUp;
                
                %for receptors that were clustered from before, make their aggregation
                %probability equal to 1, since they should stay clustered at this
                %point in the simulation
                aggregationProbVec = aggregationProb*ones(numReceptors,1);
                
                %based on these positions, cluster receptors (and modify their positions)
                [cluster2receptor,receptor2cluster,clusterSize,positionsNew] = ...
                    receptorAggregationAlg(positionsNew,aggregationDist,aggregationProbVec);
                [numClusters,maxClustSize] = size(cluster2receptor);
                
                %store new receptor information
                receptorTraj(:,:,iIter) = positionsNew;
                recept2clustAssign(:,iIter) = receptor2cluster;
                clust2receptAssign(1:numClusters,1:maxClustSize,iIter) = cluster2receptor;
                
            end %(while maxClustSize == 1 && iIter <= numIterations)
            
            if maxClustSize > 1
                time2aggregTmp(iRepConf) = (iIter-1)*timeStep;
            end
            
        end %(if maxClustSize > 1 ... else ...)
        
    end %(for iRepConf = 1 : numRepConf)
    
    %get time to aggregate for this repeat
    time2aggreg(iRep) = min(time2aggregTmp);
    
    progressText((iRep-1)/(numRep-1),'Simulation');
    
end %(for iRep = 1 : numRep)

%% ~~~ the end ~~~


%% subfunction 1

function [cluster2receptor,receptor2cluster,clusterSize,receptPositions] ...
    = receptorAggregationAlg(receptPositions,aggregationDist,aggregationProb)

%get number of receptors
numReceptors = size(receptPositions,1);

%determine whether each receptor would aggregate if given the chance
receptorAggregFlag = rand(numReceptors,1) < aggregationProb;

%calculate inter-receptor distances
interReceptDist = pdist(receptPositions);

%construct the hierarchical tree of receptors
receptTree = linkage(interReceptDist);

%construct clusters from the hierarchical tree
%receptor2cluster indicates the cluster to which each receptor belongs
receptor2cluster = cluster(receptTree,'cutoff',aggregationDist,'criterion','distance');

%get number of clusters
numClusters = max(receptor2cluster);

%put receptors that belong to the same cluster together
%first, allow every receptor to aggregate or not
cluster2receptor = zeros(numReceptors);
clusterSize = zeros(numReceptors,1);
numClustersInit = numClusters; %this is for storing receptors that refuse to aggregate
for iCluster = 1 : numClusters

    %find receptors belonging to this cluster
    clusterMembers = find(receptor2cluster==iCluster);
    numMembers = length(clusterMembers);

    %check aggregation status if there is more than one receptor in cluster
    if numMembers > 1

        %get receptor aggregation flags
        aggregFlag = receptorAggregFlag(clusterMembers);

        %if all receptors refuse to aggregate, assign the first one a flag
        %of 1 (just for algorithmic reasons)
        if all(aggregFlag == 0)
            aggregFlag(1) = 1;
        end

        %remove receptors with aggregation flag 0
        receptorsReject = clusterMembers(aggregFlag==0);
        clusterMembers = setdiff(clusterMembers,receptorsReject);
        numMembers = length(clusterMembers);

        %append to the end of the cluster vector those receptors not aggregating
        numClustersFinal = numClustersInit + length(receptorsReject);
        cluster2receptor(numClustersInit+1:numClustersFinal,1) = receptorsReject;
        clusterSize(numClustersInit+1:numClustersFinal) = 1;

        %update vector storing for every receptor its cluster number
        receptor2cluster(receptorsReject) = (numClustersInit+1:numClustersFinal)';
        numClustersInit = numClustersFinal;

    end

    %put (aggregating) receptors in their proper cluster
    cluster2receptor(iCluster,1:numMembers) = clusterMembers';
    clusterSize(iCluster) = numMembers;

end

%remove empty rows and columns and get final number of clusters
columnSum = sum(cluster2receptor);
cluster2receptor = cluster2receptor(:,columnSum~=0);
rowSum = sum(cluster2receptor,2);
cluster2receptor = cluster2receptor(rowSum~=0,:);
clusterSize = clusterSize(rowSum~=0);
numClusters = length(clusterSize);

%assign receptors the average position of the cluster they belong to
for iCluster = 1 : numClusters

    %get receptors belonging to this cluster
    clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));

    %calculate their average position
    averagePosition = mean(receptPositions(clusterMembers,:),1);

    %assign the average position to the receptors
    receptPositions(clusterMembers,:) = repmat(averagePosition,clusterSize(iCluster),1);

end

%% subfunction 2

function [cluster2receptor,receptor2cluster,clusterSize] = receptorDissociationAlg(...
    cluster2receptor,receptor2cluster,clusterSize,dissociationProb)

%get number of receptors and number of clusters
numReceptors = length(receptor2cluster);
numClusters = length(clusterSize);
maxClusterSize = max(clusterSize);

%copy some input variables for modification and output
cluster2receptorTmp = [cluster2receptor; zeros(numReceptors-numClusters,maxClusterSize)];
receptor2clusterTmp = receptor2cluster;
clusterSizeTmp = [clusterSize; zeros(numReceptors-numClusters,1)];

%determine whether each receptor would dissociate if given the chance
receptorDissociateFlag = rand(numReceptors,1) < dissociationProb;

%find clusters with more than one receptor
clustersBig = find(clusterSize>1);

%go over these clusters
numClustersInit = numClusters; %this is for storing receptors that dissociate
for iCluster = clustersBig'

    %get receptors belonging to this cluster
    clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));

    %get their dissociation flags
    dissociateFlag = receptorDissociateFlag(clusterMembers);

    %if all receptors want to dissociate, assign the first one a flag
    %of 0 (just for algorithmic reasons)
    if all(dissociateFlag == 1)
        dissociateFlag(1) = 0;
    end

    %remove receptors that want to dissocate from cluster
    receptorsReject = clusterMembers(dissociateFlag==1);
    clusterMembers = setdiff(clusterMembers,receptorsReject);
    numMembers = length(clusterMembers);

    %append to the end of the cluster vector those receptors that
    %dissociated
    numClustersFinal = numClustersInit + length(receptorsReject);
    cluster2receptorTmp(numClustersInit+1:numClustersFinal,1) = receptorsReject;
    clusterSizeTmp(numClustersInit+1:numClustersFinal) = 1;

    %keep the receptors that did not dissociate in their proper cluster
    cluster2receptorTmp(iCluster,:) = 0;
    cluster2receptorTmp(iCluster,1:numMembers) = clusterMembers';
    clusterSizeTmp(iCluster) = numMembers;

    %update vector storing for every receptor its cluster number
    receptor2clusterTmp(receptorsReject) = (numClustersInit+1:numClustersFinal)';
    numClustersInit = numClustersFinal;

end

%remove empty row and columns and get final number of clusters
columnSum = sum(cluster2receptorTmp);
cluster2receptorTmp = cluster2receptorTmp(:,columnSum~=0);
rowSum = sum(cluster2receptorTmp,2);
cluster2receptorTmp = cluster2receptorTmp(rowSum~=0,:);
clusterSizeTmp = clusterSizeTmp(rowSum~=0);

%copy temporary variable into output variables
cluster2receptor = cluster2receptorTmp;
receptor2cluster = receptor2clusterTmp;
clusterSize = clusterSizeTmp;


%% THIS WAS FOR receptorAggregationAlg

%%% OLD CODE THAT WORKED, BUT THAT I REALIZED WAS UNNECESSARY SINCE I COULD
%%% USE MATLAB'S LINKAGE AND CLUSTER FUNCTIONS, WHICH WERE MUCH FASTER.

% % % %calculate inter-receptor distances
% % % interReceptDist = createDistanceMatrix(receptPositions,receptPositions);
% % % for i = 1 : numReceptors
% % %     interReceptDist(i:numReceptors,i) = NaN;
% % % end
% % % 
% % % %find distances smaller than the "aggregation distance" - these are the
% % % %receptors that are considered to have bumped with each other
% % % receptorBumpingMat = interReceptDist < aggregationDist;
% % % 
% % % %determine whether the receptors that bumped into each other aggregated
% % % receptorAggregMat = (rand(numReceptors) < aggregationProb) .* ...
% % %     receptorBumpingMat;
% % % 
% % % %make an initial set of clusters where each receptors is on its own
% % % cluster2receptor = zeros(numReceptors);
% % % cluster2receptor(:,1) = (1:numReceptors)';
% % % 
% % % %go over all receptors and assign them to initial clusters
% % % for iReceptor = 1 : numReceptors
% % %     
% % %     %find receptors that this receptor aggregates with
% % %     associateReceptors = find(receptorAggregMat(iReceptor,:)==1);
% % %     
% % %     %if it has associate receptors, move it to their location and empty its
% % %     %location
% % %     if ~isempty(associateReceptors)
% % %         for iAssociate = associateReceptors
% % %             startCol = find(cluster2receptor(iAssociate,:)==0,1,'first');
% % %             cluster2receptor(iAssociate,startCol:end) = ...
% % %                 cluster2receptor(iReceptor,1:end-startCol+1);
% % %         end
% % %         cluster2receptor(iReceptor,:) = 0;
% % %     end
% % %     
% % % end
% % % 
% % % %remove all empty clusters
% % % cluster2receptor = cluster2receptor(cluster2receptor(:,1)~=0,:);
% % % 
% % % %get number of initial clusters
% % % numClusters = size(cluster2receptor,1);
% % % 
% % % %get clusters that have more than one member
% % % clusterSize = zeros(numClusters,1);
% % % for iCluster = 1 : numClusters
% % %     clusterSize(iCluster) = length(find(cluster2receptor(iCluster,:)~=0));
% % % end
% % % clustersBig = find(clusterSize>1)';
% % % numClustersBig = length(clustersBig);
% % % 
% % % %go over all clusters with more than one member, merge those which have
% % % %common members, and remove any member repetition in each cluster
% % % for iTmp = 1 : numClustersBig
% % %     
% % %     %get current cluster index
% % %     iCluster = clustersBig(iTmp);
% % % 
% % %     %get members in this cluster
% % %     clusterMembers = cluster2receptor(iCluster,cluster2receptor(iCluster,:)~=0);
% % %     
% % %     %check whether any of the following clusters have its members repeated
% % %     clustersBigTmp = clustersBig(iTmp+1:numClustersBig);
% % %     cluster2receptorTmp = cluster2receptor(clustersBigTmp,:);
% % %     memberRep = intersect(clusterMembers,cluster2receptorTmp(:));
% % %     
% % %     %proceed if any of the members are repeated in later clusters
% % %     if ~isempty(memberRep)
% % % 
% % %         %go over all later clusters
% % %         for jTmp = iTmp+1 : numClustersBig
% % % 
% % %             %get second cluster index
% % %             jCluster = clustersBig(jTmp);
% % % 
% % %             %get second cluster members
% % %             cluster2ndMembers = cluster2receptor(jCluster,cluster2receptor(jCluster,:)~=0);
% % % 
% % %             %if there are any common members, merge clusters into the second,
% % %             %erase first cluster
% % %             if ~isempty(intersect(clusterMembers,cluster2ndMembers))
% % %                 unionMembers = union(clusterMembers,cluster2ndMembers);
% % %                 cluster2receptor(jCluster,:) = 0;
% % %                 cluster2receptor(jCluster,1:length(unionMembers)) = unionMembers;
% % %                 cluster2receptor(iCluster,:) = 0;
% % %             end
% % % 
% % %         end
% % % 
% % %     end
% % %     
% % % end
% % %       
% % % %again remove all empty clusters
% % % cluster2receptor = cluster2receptor(cluster2receptor(:,1)~=0,:);
% % % 
% % % %get number of final clusters
% % % numClusters = size(cluster2receptor,1);
% % % 
% % % %calculate cluster sizes
% % % clusterSize = zeros(numClusters,1);
% % % for iCluster = 1 : numClusters
% % %     clusterSize(iCluster) = length(find(cluster2receptor(iCluster,:)~=0));
% % % end
% % % clustersBig = find(clusterSize>1)';
% % %   
% % % %remove any receptor repetition in any cluster
% % % for iCluster = clustersBig
% % % 
% % %     %get members in this cluster
% % %     clusterMembers = cluster2receptor(iCluster,cluster2receptor(iCluster,:)~=0);
% % %     
% % %     %remove any repetition and get new cluster size
% % %     clusterMembers = unique(clusterMembers);
% % %     clusterSize(iCluster) = length(clusterMembers);
% % %     
% % %     %put back in big matrix
% % %     cluster2receptor(iCluster,:) = 0;
% % %     cluster2receptor(iCluster,1:clusterSize(iCluster)) = clusterMembers;
% % %     
% % % end
% % % 
% % % %remove all unused columns
% % % columnSum = sum(cluster2receptor);
% % % cluster2receptor = cluster2receptor(:,columnSum~=0);
% % % 

