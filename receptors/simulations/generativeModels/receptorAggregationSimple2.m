function [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
    = receptorAggregationSimple(modelParam,simParam)
%RECEPTORAGGREGATIONSIMPLE simulates the aggregation of simply diffusing receptors
%
%SYNOPSIS [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
%    = receptorAggregationSimple(modelParam,simParam)
%
%INPUT  modelParam: Structure with the fields:
%           diffCoef        : Diffusion coefficient (microns^2/s).
%           receptorDensity : Receptor density (#/microns^dimension).
%           aggregationProb : Probability of aggregation if a receptor
%                             bumps into another receptor or receptor
%                             complex.
%           aggregationDist : Distance between 2 receptors to consider them
%                             bumping into each other, thus potentially
%                             aggregating.
%           dissociationRate: Rate of dissociation of a receptor from a
%                             receptor complex (#/s).
%           labelRatio      : Receptor labeling ratio.
%           intensityQuantum: Row vector with 2 values, the mean and std of
%                             receptor labeling intensities.
%       simParam: Structure with the fields:
%           dimension       : System dimensionality (1, 2 or 3). Default: 2.
%           observeSideLen  : Observation side length (microns). Either one
%                             value, used for all dimensions, or a row
%                             vector with a value for each dimension.
%                             Default: 1 in all dimensions.
%           timeStep        : Simulation time step (s).
%                             Default: 0.01/diffCoef.
%           totalTime       : Total time of simulation (s).
%                             Default:  100 * timeStep.
%           randNumGenSeeds : Row vector storing the seeds for the "rand"
%                             and "randn" random number generators.
%                             Default: [100 100].
%                 Whole structure optional. Individual fields also optional.
%
%OUTPUT receptorInfoAll     : Structure with fields:
%           .receptorTraj        : (Number of receptors) - by - (dimension)
%                                  - by - (number of iterations) array storing
%                                  receptor positions.
%           .recept2clustAssign  : (Number of receptors) - by - (number of
%                                  iterations) array storing the cluster
%                                  assignment of each receptor.
%           .clust2receptAssign  : (Maximum number of clusters) - by - (maximum
%                                  cluster size) - by - (number of iterations)
%                                  array storing the receptor membership to each
%                                  cluster.
%       receptorInfoLabeled : The same as receptorInfoAll, but for labaled
%                             receptors only. It also has the additional
%                             field:
%           .compTracks          : The receptor trajectories and
%                                  interactions over time, as converted by
%                                  convReceptClust2CompTracks. Same format
%                                  as output of trackCloseGapsKalman.
%       timeIterArray       : Vector storing time corresponding to each
%                             iteration.
%       errFlag : 0 if function executes normally, 1 otherwise.

%Khuloud Jaqaman, January 2009

%% Output

errFlag = 0;
receptorInfoAll = [];
receptorInfoLabeled = [];
timeIterArray = [];

%% Input

%check if correct number of arguments were used when function was called
if nargin < 1
    disp('--receptorAggregationSimple: Too few input arguments');
    errFlag  = 1;
    return
end

%extract model parameters from modelParam

%diffusion coefficient
if isfield(modelParam,'diffCoef')
    diffCoef = modelParam.diffCoef;
else
    disp('--receptorAggregationSimple: Please supply diffusion coefficient');
    errFlag = 1;
end

%receptor density
if isfield(modelParam,'receptorDensity')
    receptorDensity = modelParam.receptorDensity;
else
    disp('--receptorAggregationSimple: Please supply receptor density');
    errFlag = 1;
end

%aggregation probability
if isfield(modelParam,'aggregationProb')
    aggregationProb = modelParam.aggregationProb;
else
    disp('--receptorAggregationSimple: Please supply aggregation probability');
    errFlag = 1;
end

%aggregation distance
if isfield(modelParam,'aggregationDist')
    aggregationDist = modelParam.aggregationDist;
else
    disp('--receptorAggregationSimple: Please supply aggregation distance');
    errFlag = 1;
end

%dissociation rate
if isfield(modelParam,'dissociationRate')
    dissociationRate = modelParam.dissociationRate;
else
    disp('--receptorAggregationSimple: Please supply dissociation rate');
    errFlag = 1;
end

%labeling ratio
if isfield(modelParam,'labelRatio')
    labelRatio = modelParam.labelRatio;
else
    disp('--receptorAggregationSimple: Please supply labeling ratio');
    errFlag = 1;
end

%intensity quantum
if isfield(modelParam,'intensityQuantum')
    intensityQuantum = modelParam.intensityQuantum;
else
    disp('--receptorAggregationSimple: Please supply intensity quantum');
    errFlag = 1;
end

%exit if there are missing model parameters
if errFlag == 1
    disp('--receptorAggregationSimple: Please supply missing variables!');
    return
end

%check model parameters ...

%some must be positive
if any([receptorDensity aggregationDist labelRatio intensityQuantum] <= 0)
    disp('--receptorAggregationSimple: Receptor density, aggregation distance, labeling ratio and intensity quantum should be positive');
    errFlag = 1;
    return
end

%and some must be non-negative
if any([diffCoef aggregationProb dissociationRate] < 0)
    disp('--receptorAggregationSimple: Diffusion coefficient, aggregation probability and dissociation rate should be non-negative');
    errFlag = 1;
    return
end

%extract simulation parameters from simParam

%if simParam wasn't supplied at all
if nargin < 2 || isempty(simParam)

    dimension = 2;
    observeSideLen = ones(1,dimension);
    timeStep = 0.01 / diffCoef;
    totalTime = 100 * timeStep;
    randNumGenSeeds = [100 100];

else

    %system dimensionality
    if isfield(simParam,'dimension')
        dimension = simParam.dimension;
    else
        dimension = 2;
    end

    %observation side length
    if isfield(simParam,'observeSideLen')
        observeSideLen = simParam.observeSideLen;
        if length(observeSideLen) == 1
            observeSideLen = observeSideLen * ones(1,dimension);
        end
    else
        observeSideLen = ones(1,dimension);
    end

    %time step
    if isfield(simParam,'timeStep')
        timeStep = simParam.timeStep;
    else
        timeStep = 0.01 / diffCoef;
    end

    %total simulation time
    if isfield(simParam,'totalTime')
        totalTime = simParam.totalTime;
    else
        totalTime = 100 * diffCoef;
    end

    %random number generator seeds
    if isfield(simParam,'randNumGenSeeds')
        randNumGenSeeds = simParam.randNumGenSeeds;
    else
        randNumGenSeeds = [100 100];
    end

end

%determine number of iterations to perform
numIterations = ceil( totalTime / timeStep ) + 1;

%store time correponding to each iteration
timeIterArray = (0 : numIterations - 1)' * timeStep;

%assuming that the displacement taken in a timeStep is normally distributed
%with mean zero, get standard deviation of step distribution from diffusion
%coefficient
stepStd = sqrt( 2 * diffCoef * timeStep );

%calculate dissociation probability
dissociationProb = dissociationRate * timeStep;

%initialize random number generators
rand('twister',randNumGenSeeds(1));
randn('state',randNumGenSeeds(2));

%% receptor initial positions and clustering

%calculate observation region size
obsRegionSize = prod(observeSideLen);

%calculate number of receptors
numReceptors = round(obsRegionSize * receptorDensity);

%initialize receptor positions
initPositions = rand(numReceptors,dimension) .* repmat(observeSideLen,numReceptors,1);

%based on these positions, cluster receptors (and modify their positions)
aggregationProbVec = aggregationProb*ones(numReceptors,1);
[cluster2receptor,receptor2cluster,clusterSize,initPositions] = ...
    receptorAggregationAlg(initPositions,aggregationDist,aggregationProbVec);
[numClusters,maxClustSize] = size(cluster2receptor);

%% Trajectory generation

%reserve memory for output vectors
receptorTraj = zeros(numReceptors,dimension,numIterations);
recept2clustAssign = zeros(numReceptors,numIterations);
clust2receptAssign = zeros(numReceptors,maxClustSize,numIterations);

%store initial information
receptorTraj(:,:,1) = initPositions;
recept2clustAssign(:,1) = receptor2cluster;
clust2receptAssign(1:numClusters,1:maxClustSize,1) = cluster2receptor;

%iterate in time
for iIter = 2 : numIterations

    %allow receptors in clusters to dissociate in current time point
    [cluster2receptor,receptor2cluster,clusterSize] = receptorDissociationAlg(...
        cluster2receptor,receptor2cluster,clusterSize,dissociationProb);
    
    %get indices of clusters with more than one receptor
    clustersBig = find(clusterSize>1);

    %get receptor positions at previous time point
    positionsOld = receptorTraj(:,:,iIter-1);
    
    %generate receptor displacements
    %     receptorDisp = stepStd*randn(numReceptors,dimension);
    receptorDisp = randn(numReceptors,dimension) .* ...
        repmat(stepStd*(1.5-(abs(positionsOld(:,2)-0.1))/0.1),1,dimension);

    %assign receptors in a cluster the displacement of the receptor with
    %the smallest index
    for iCluster = clustersBig'

        %get receptors belonging to this cluster
        clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));

        %find the receptor with the smallest index
        indxSmall = min(clusterMembers);

        %assign all receptors in this cluster the displacement of that
        %receptor
        receptorDisp(clusterMembers,:) = repmat(receptorDisp(indxSmall,:),...
            clusterSize(iCluster),1);

    end
    
    %calculate the new receptor positions at the current time point
    positionsNew = positionsOld + receptorDisp;
    
    %make sure that receptors stay inside the region of interests
    correctionBoundaryLow = min(positionsNew,0);
    positionsNew = positionsNew - 2 * correctionBoundaryLow;
    correctionBoundaryUp = max(positionsNew - repmat(observeSideLen,numReceptors,1),0);
    positionsNew = positionsNew - 2 * correctionBoundaryUp;

    %for receptors that were clustered from before, make their aggregation
    %probability equal to 1, since they should stay clustered at this point
    %in the simulation
    aggregationProbVec = aggregationProb*ones(numReceptors,1);
    clusteredReceptorIndx = cluster2receptor(clustersBig,:);
    clusteredReceptorIndx = clusteredReceptorIndx(clusteredReceptorIndx~=0);
    aggregationProbVec(clusteredReceptorIndx) = 1;

    %based on these positions, cluster receptors (and modify their positions)
    [cluster2receptor,receptor2cluster,clusterSize,positionsNew] = ...
        receptorAggregationAlg(positionsNew,aggregationDist,aggregationProbVec);
    [numClusters,maxClustSize] = size(cluster2receptor);
    
    %store new receptor information
    receptorTraj(:,:,iIter) = positionsNew;
    recept2clustAssign(:,iIter) = receptor2cluster;
    clust2receptAssign(1:numClusters,1:maxClustSize,iIter) = cluster2receptor;

end %(for iIter = 2 : numIterations)

%remove empty rows and columns from clust2receptAssign
cluster2receptor = max(clust2receptAssign,[],3);
columnSum = sum(cluster2receptor);
clust2receptAssign = clust2receptAssign(:,columnSum~=0,:);
rowSum = sum(cluster2receptor,2);
clust2receptAssign = clust2receptAssign(rowSum~=0,:,:);

% %put receptor trajectories and clusters into the format of the output of
% %trackCloseGapsKalman
% compTracks = convReceptClust2CompTracks(clust2receptAssign,...
%     recept2clustAssign,receptorTraj);

%put information in receptorInfoAll
receptorInfoAll = struct('receptorTraj',receptorTraj,'recept2clustAssign',...
    recept2clustAssign,'clust2receptAssign',clust2receptAssign);

%% Receptor labeling and sub-sampling

%if labeling ratio is less than one ...
if labelRatio < 1

    %label some receptors based on the labeling ratio
    %1 = labeled, 0 = unlabeled
    labelFlag = rand(numReceptors,1) <= labelRatio;
    indxLabeled = find(labelFlag==1);
    indxNotLabeled = find(labelFlag==0);
    
    %extract the trajectories of the labeled receptors
    receptorTrajLabeled = receptorTraj(indxLabeled,:,:);
    
    %assign the intensities of the labeled receptors
    receptorIntensityLabeled = intensityQuantum(1) + ...
        randn(length(indxLabeled),numIterations) * intensityQuantum(2);
    
    %extract the cluster assignments of the labeled receptors
    recept2clustAssignLabeled = recept2clustAssign(indxLabeled,:);
    
    %modify the cluster-to-receptor assignments to include only labeled
    %receptors
    clust2receptAssignLabeled = clust2receptAssign;
    %convert receptors that are not labeled to zero
    for iNotLabeled = indxNotLabeled'
        clust2receptAssignLabeled(clust2receptAssignLabeled==iNotLabeled) = 0;
    end
    %update the indices of the labeled receptors
    for iLabeled = 1 : length(indxLabeled)
        clust2receptAssignLabeled(clust2receptAssignLabeled==indxLabeled(iLabeled)) = iLabeled;
    end
    %make sure that zeros come after receptor indices
    clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2,'descend');
    
    %remove empty clusters (which include unlabeled receptors)
    %modify receptor-to-cluster assignments accordingly
    for iIter = 1 : numIterations
        clustSize = sum(clust2receptAssignLabeled(:,:,iIter)~=0,2);
        indxFull = find(clustSize~=0);
        indxEmpty = find(clustSize==0);
        clust2receptAssignLabeled(:,:,iIter) = clust2receptAssignLabeled(...
            [indxFull;indxEmpty],:,iIter);
        for iFull = 1 : length(indxFull)
            recept2clustAssignLabeled(recept2clustAssignLabeled(:,iIter)...
                ==indxFull(iFull),iIter) = iFull;
        end
    end

    %remove empty rows and columns from clust2receptAssign
    cluster2receptor = max(clust2receptAssignLabeled,[],3);
    columnSum = sum(cluster2receptor);
    clust2receptAssignLabeled = clust2receptAssignLabeled(:,columnSum~=0,:);
    rowSum = sum(cluster2receptor,2);
    clust2receptAssignLabeled = clust2receptAssignLabeled(rowSum~=0,:,:);

    %put labeled receptor trajectories and clusters into the format of the
    %output of trackCloseGapsKalman
    compTracksLabeled = convReceptClust2CompTracks(clust2receptAssignLabeled,...
        recept2clustAssignLabeled,receptorTrajLabeled,receptorIntensityLabeled);

    %put information in receptorInfoLabeled
    receptorInfoLabeled = struct('receptorTraj',receptorTrajLabeled,...
        'recept2clustAssign',recept2clustAssignLabeled,...
        'clust2receptAssign',clust2receptAssignLabeled,...
        'compTracks',compTracksLabeled);

else %if all receptors are labeled
    
    %copy receptor information from all to labaled
    receptorInfoLabeled = receptorInfoAll;

    %assign receptor intensities
    receptorIntensity = intensityQuantum(1) + ...
        randn(numReceptors,numIterations) * intensityQuantum(2);
    
    %put labeled receptor trajectories and clusters into the format of the
    %output of trackCloseGapsKalman
    compTracksLabeled = convReceptClust2CompTracks(clust2receptAssign,...
        recept2clustAssign,receptorTraj,receptorIntensity);
    
    %store the new compound tracks with the labeling intensity
    receptorInfoLabeled.compTracks = compTracksLabeled;

end %(if labelRatio < 1 ... else ...)


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
%receptorMat2Clusters indicates the cluster to which each receptor belongs
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

%remove empty row and columns and get final number of clusters
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

