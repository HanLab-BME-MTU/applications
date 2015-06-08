function [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag,assocStats,collProbStats] ...
    = receptorAggregationSimple_new(modelParam,simParam)
%RECEPTORAGGREGATIONSIMPLE simulates the aggregation of simply diffusing receptors
%
%SYNOPSIS [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
%    = receptorAggregationSimple(modelParam,simParam)
%
%INPUT  modelParam: Structure with the fields:
%           diffCoef        : Diffusion coefficient (microns^2/s).
%           receptorDensity : Receptor density (#/microns^probDim).
%           aggregationProb : Probability of aggregation if a receptor
%                             bumps into another receptor or receptor
%                             complex.
%           aggregationDist : Distance between 2 receptors to consider them
%                             bumping into each other, thus potentially
%                             aggregating (microns).
%           dissociationRate: Rate of dissociation of a receptor from a
%                             receptor complex (#/s).
%           labelRatio      : Receptor labeling ratio.
%           intensityQuantum: Row vector with 2 values, the mean and std of
%                             receptor labeling intensities.
%           initPositions   : Receptor initial positions. If supplied, they
%                             will be used. If not, random positions will
%                             be chosen.
%       simParam: Structure with the fields:
%           probDim         : System dimensionality (1, 2 or 3). Default: 2.
%           observeSideLen  : Observation side length (microns). Either one
%                             value, used for all probDims, or a row
%                             vector with a value for each probDim.
%                             Default: 1 in all probDims.
%           timeStep        : Simulation time step (s).
%                             Default: 0.01/max(diffCoef,dissociationRate).
%           simTime         : Total time of simulation (s).
%                             Default: 100 * timeStep.
%           initTime        : Initialization time, to remove any
%                             initialization artifacts.
%                             Default: 100 * timeStep.
%           randNumGenSeeds : Row vector storing the seeds for the "rand"
%                             and "randn" random number generators.
%                             Default: [100 100].
%                 Whole structure optional. Individual fields also optional.
%
%OUTPUT receptorInfoAll     : Structure with fields:
%           .receptorTraj        : (Number of receptors) - by - (probDim)
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
%
%       assocStats: a struct containing details on sure, potential and
%                   probability based associations, at each iteraion.
%                   The fields are:
%                   numSureAssoc - number of sure associations
%                   numPotColl - number of potentail collisions
%                   numColl - number of collision from potential collisions
%                   numPotColl_Assoc - number of actual association that
%                                      resulted from potentail collisions
%                   sureAssocCountBySize - sure associations by cluster
%                                          size
%                   numCollProbPairs - number of pairs in collision via
%                                      probablity
%
%       collProbStats:  statistics specifically for collisions based on
%                       probablity at each iteration. This is a struct with
%                       fields:
%                       collisionProb - the probablity value
%                       pwDist - pairwise distance
%                       primaryNodeRadius - radius for primary node
%                       partnerNodeRadii - set of radii for partner nodes
%
%
%Khuloud Jaqaman, January 2009
%
%Modified June 2013 - December 2014, Robel Yirdaw
%

%% Output

errFlag = 0;
receptorInfoAll = [];
%09/05/14 (ryirdaw)
%need to block the following otherwise conversion from struct to double
%error at at the end
%receptorInfoLabeled = [];
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

%receptor initial positions
if isfield(modelParam,'initPositions')
    initPositions = modelParam.initPositions;
else
    initPositions = [];
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
%09/05/14 - modified to accomodate a vector labelRatio
if any([receptorDensity aggregationDist (labelRatio') intensityQuantum(1)] <= 0)
    disp('--receptorAggregationSimple: Receptor density, aggregation distance, labeling ratio and intensity quantum should be positive');
    errFlag = 1;
    return
end


%and some must be non-negative
%03/25/14 - modified to accomodate a vector aggregationProb
if any([diffCoef (aggregationProb') dissociationRate intensityQuantum(2)] < 0)
    disp('--receptorAggregationSimple: Diffusion coefficient, aggregation probability and dissociation rate should be non-negative');
    errFlag = 1;
    return
end

%extract simulation parameters from simParam

%if simParam wasn't supplied at all
if nargin < 2 || isempty(simParam)

    probDim = 2;
    observeSideLen = ones(1,probDim);
    timeStep = 0.01 / max(diffCoef,dissociationRate);
    simTime = 100 * timeStep;
    initTime = 100 * timeStep;
    randNumGenSeeds = [100 100];

else

    %system probDimality
    if isfield(simParam,'probDim')
        probDim = simParam.probDim;
    else
        probDim = 2;
    end

    %observation side length
    if isfield(simParam,'observeSideLen')
        observeSideLen = simParam.observeSideLen;
        if length(observeSideLen) == 1
            observeSideLen = observeSideLen * ones(1,probDim);
        end
    else
        observeSideLen = ones(1,probDim);
    end

    %time step
    if isfield(simParam,'timeStep')
        timeStep = simParam.timeStep;
    else
        timeStep = 0.01 / max(diffCoef,dissociationRate);
    end

    %simulation time
    if isfield(simParam,'simTime')
        simTime = simParam.simTime;
    else
        simTime = 100 * timeStep;
    end

    %initialization time
    if isfield(simParam,'initTime')
        initTime = simParam.initTime;
    else
        initTime = 100 * timeStep;
    end

    %random number generator seeds
    if isfield(simParam,'randNumGenSeeds')
        randNumGenSeeds = simParam.randNumGenSeeds;
    else
        randNumGenSeeds = [100 100];
    end

end

%determine number of iterations to perform
totalTime = simTime + initTime;
numIterSim = ceil( simTime / timeStep ) + 1;
numIterInit = ceil( initTime / timeStep );
numIterations = numIterSim + numIterInit;

%store time correponding to each iteration
timeIterArray = (0 : numIterSim - 1)' * timeStep;

%assuming that the displacement taken in a timeStep is normally distributed
%with mean zero, get standard deviation of step distribution from diffusion
%coefficient
stepStd = sqrt( 2 * diffCoef * timeStep );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%05/27/14 (ryirdaw)
%We are now introducing a collision probability modeled by a linear
%function with slope and interecept values calculated below. The
%probability is 1 at the user defined association distance
%(associationDist) and decreases linearly to reach 0 at the corrected
%association distance (associationDistCorr).
%Also, the corrected aggregationDist is now saved as aggregationDistCorr, 
%since both need to be passed to the association function which is also
%modified to accomodate these chagnes

%adjust aggregationDist to account for the finite simulation time step and
%the expected receptor displacement in that time step
%NOTE: 1/2 corrected to 2x (10/2013)
aggregationDistCorr = max(aggregationDist,sqrt(probDim)*stepStd*2);

%05/27/14 (ryirdaw) - use a struct for the associationDist values
associationDistVals = struct('associationDist',aggregationDist,...
    'associationDistCorr',aggregationDistCorr);
%{
%05/27/14 (ryirdaw) - calculate slope and intercept values of the collision probability
%model and place in a struct.
%Slope

%m = 1 / (aggregationDist - aggregationDistCorr);
m = 0;
%Intercept
%b = -1*m*aggregationDistCorr;
b = 1;
%Place values in a struct
collisionProbVals = struct('m',m,'b',b);
%}

%06/04/14 (ryirdaw)
%Implemented collision probability via the arc approach. To use this, set
%collisionProbVals to [], so that the linear model approach will be
%disabled. The above assignment should be blocked as well.
%collisionProbVals = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%06/24/14 (ryirdaw)
%Determine the number of displacements for closest distance/collision
%calculations of nodes
%{
%Calculate unit time of displacment from diffusion coefficient and
%aggregationDist via the rmsd formula. The actual value is rounded down to
%the nearest 1E-X, where X is the OM of the value determined at aggregDist.
%dispUnitTime = 1e-6;
dispUnitTime = 10^(floor(log10( (1/diffCoef)*(aggregationDist/4)^2 ) ) );
%Number of displacements
numDispPerUnitTime = timeStep/dispUnitTime;
%Time vector for each displacement
dispTimeVec = 0:dispUnitTime:timeStep;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate dissociation probability
dissociationProb = dissociationRate * timeStep;

%initialize random number generators
% rand('twister',randNumGenSeeds(1));
% randn('state',randNumGenSeeds(2));
rng(randNumGenSeeds(1),'twister')

%% receptor initial positions and clustering

%calculate observation region size
obsRegionSize = prod(observeSideLen);

if isempty(initPositions)
    
    %calculate number of receptors
    numReceptors = round(obsRegionSize * receptorDensity);
    
    %initialize receptor positions
    initPositions = rand(numReceptors,probDim) .* repmat(observeSideLen,numReceptors,1);
    
else
    
    %get number of receptors
    numReceptors = size(initPositions,1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%03/25/14 (ryirdaw)
%aggregationProb is now a vector giving prob of formation per size - first
%element is NaN. This changes the initial cluster formation as implemented
%in the original code (below).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
02/21/14 (ryirdaw)
Disabled the following call to generate clusters at the start of the
simulation - starting with all monomers instead. The 4 return vars are
needed for the subsequent lines so need to set here all except
initPositions. 
Additional Note:
One easy check for this would be to set aggregatinProbVec to
0 and let it continue as before. This will also start the sim. with all
free receptors but note that after return from receptorAggregationAlg the
cluster IDs will be scrambled and receptorDisp in the for-loop below will
produce a different set of displacement values than obtained if taking the
new approach (setting the three variables manually below). This is because
the sequence of rand values is different for the two since the call to
receptorAggregationAlg will result in rand calls.


%Use the following to start with random clusters when aggregationProb = 0.
aggregationProbVec = ones(numReceptors,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original. Modified 03/25/14 - aggregationProb can now be a vector in which
%case we set all receptor's aggregationProbVec with 1.

%based on these positions, cluster receptors (and modify their positions)
if (length(aggregationProb) == 1)
    aggregationProbVec = aggregationProb*ones(numReceptors,1);
else
    aggregationProbVec = ones(numReceptors,1);
end

[cluster2receptor,receptor2cluster,clusterSize,initPositions] = ...
    receptorAggregationAlg(initPositions,aggregationDist,aggregationProbVec,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%Starting with all receptors as monomers.  Set cluster2receptor and
%receptor2cluster as 1D vector of 1:numReceptors.  Then clusterSize will be
%numReceptors long with all values of 1. Note all are column vectors.
receptor2cluster = (1:numReceptors)';
cluster2receptor = (1:numReceptors)';
clusterSize = ones(numReceptors,1);

[numClusters,maxClustSize] = size(cluster2receptor);

%% Main simulation body

%reserve memory for output vectors
receptorTraj = zeros(numReceptors,probDim,numIterations);
recept2clustAssign = zeros(numReceptors,numIterations);
clust2receptAssign = zeros(numReceptors,maxClustSize,numIterations);

%store initial information
receptorTraj(:,:,1) = initPositions;
recept2clustAssign(:,1) = receptor2cluster;
clust2receptAssign(1:numClusters,1:maxClustSize,1) = cluster2receptor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnostic quantities (ryirdaw)
assocStats = struct('numSureAssoc',NaN(numIterations,1),...
    'numPotColl',NaN(numIterations,1),'numColl',NaN(numIterations,1),...
    'numPotColl_Assoc',NaN(numIterations,1),...
    'sureAssocCountBySize',NaN(numReceptors,numIterations),...
    'numCollProbPairs',NaN(numIterations,1));

collProbStatStruct = struct('collisionProb',NaN,'pwDist',NaN,...
    'primaryNodeRadius',NaN,'partnerNodeRadii',NaN);
collProbStats(numIterations,1) = collProbStatStruct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

progressText(0,'Simulation');

%iterate in time
for iIter = 2 : numIterations
    %fprintf('\niIter = %d\n',iIter);
    
    %% Dissociation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %06/24/13 (ryirdaw)
    %the dissociation probability is normalized to the number of receptors
    %in a given cluster - each receptor has an effective dissociation
    %probability.  Pass this to receptorDissociationAlg (below). 

    %numReceptorsPerCluster = sum( (cluster2receptor(receptor2cluster,:) ~= 0),2);
    %denom = factorial(numReceptorsPerCluster)./(factorial(numReceptorsPerCluster - 1));
    %effDissociationProb = dissociationProb./numReceptorsPerCluster;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if (any(numReceptorsPerCluster > 1))
    %allow receptors in clusters to dissociate in current time point
    [cluster2receptor,receptor2clusterDissAlg,clusterSize] = receptorDissociationAlg(...
        cluster2receptor,receptor2cluster,clusterSize,dissociationProb);
    %end
    
    %% Association flags based on outcome of dissociation
   
    %for receptors that were clustered from before, make their aggregation
    %probability equal to 1, since they should stay clustered at this
    %point in the simulation
    %aggregationProbVec = aggregationProb*ones(numReceptors,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %03/25/14 (ryirdaw)
    %This no longer depends on aggregationProb since aggregationProb is a
    %vector. Just set to 1 here and will be modified below for those
    %receptors that just underwent dissociaiton. The final vector is passed
    %to the associaition function.
    aggregationProbVec = ones(numReceptors,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %01/28/14 (ryirdaw)
    %Disable when using receptorAggregationAlg_maxWeightedMatching as there
    %is no need - clusters will not dissociate in the above function unlike
    %the original (must enable if using the original aggreg. function below)
    %{
    clusteredReceptorIndx = cluster2receptor(clustersBig,:);
    clusteredReceptorIndx = clusteredReceptorIndx(clusteredReceptorIndx~=0);
    aggregationProbVec(clusteredReceptorIndx) = 1;
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %09/05/13 (ryirdaw)
    %The call to receptorDissociationAlg above could have dissociated
    %receptors from existing clusters.  These free receptors should not be
    %allowed to associate in receptorAggregationAlg below. Set the
    %aggregation probability for these receptors to zero.
    %Note: receptor2cluster vector from receptorDissociationAlg is being
    %recieved now (above, receptor2clusterDissAlg).
    
    if (max(receptor2clusterDissAlg) > max(receptor2cluster))
        %A dissociation has occured. To confirm and identify the 
        %receptors invovled, determine the cluster size of each receptor
        %now and compare with size values from previous iterations. The
        %best way to do this is to tally the number of other receptors each
        %receptor is associated with.
        %Store the new and previous sizes on two columns for each receptor.        

        sizeNewPrev = NaN(length(receptor2clusterDissAlg),2);
        for recepID=1:length(receptor2clusterDissAlg)
            %Get cluster label for current receptor in this and previous
            %iteration
            clustIDNew = receptor2clusterDissAlg(recepID);
            clustIDPrev = receptor2cluster(recepID);
            %Get the total number of receptors associated
            sizeNewPrev(recepID,1:2) = [sum(receptor2clusterDissAlg == clustIDNew)...
                sum(receptor2cluster == clustIDPrev)];
        end
        %For those receptors who have dissociated set the
        %aggregationProbVec to 0.  NOTE: if the other receptors remain
        %clustered, their aggregationProbVec must stay as 1.        
        %{
        aggregationProbVec( (sizeNewPrev(:,1) - sizeNewPrev(:,2) < 0) & ...
            (sizeNewPrev(:,1) == 1) ) = 0;
        %}
        aggregationProbVec( (sizeNewPrev(:,1) - sizeNewPrev(:,2) < 0) ) = 0;
        
        %Reassign receptor2cluster to the new set reflecting the
        %dissociation.  NOTE: the position vector also shows the
        %dissociation that has occured.
        receptor2cluster = receptor2clusterDissAlg;
    end
    
    %% New receptor/cluster positions
    
    %get indices of clusters with more than one receptor
    clustersBig = find(clusterSize>1);

    %get receptor positions at previous time point
    positionsOld = receptorTraj(:,:,iIter-1);
    
    %generate receptor displacements
    receptorDisp = stepStd*randn(numReceptors,probDim);
    %     receptorDisp = randn(numReceptors,probDim) .* ...
    %         repmat(stepStd*(1.5-(abs(positionsOld(:,2)-0.1))/0.1),1,probDim);

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
    
    %make sure that receptors stay inside the region of interest
    correctionBoundaryLow = min(positionsNew,0);
    positionsNew = positionsNew - 2 * correctionBoundaryLow;
    correctionBoundaryUp = max(positionsNew - repmat(observeSideLen,numReceptors,1),0);
    positionsNew = positionsNew - 2 * correctionBoundaryUp;

    %% Association
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %07/03/14 (ryirdaw)
    %Perform association of nodes whose pairwise distance < aggregationDist
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %{
    %01/28/14 (ryirdaw)
    %If using this, must enable lines above that set clusters flags to 1
    
    %based on these positions, cluster receptors (and modify their positions)
    [cluster2receptor,receptor2cluster,clusterSize,positionsNew] = ...
        receptorAggregationAlg(positionsNew,aggregationDist,aggregationProbVec,receptor2cluster);    
    %}
    
    %{
    03/25/14 (ryirdaw)
    The association function has been modified to allow for variable 
    association probabilities, per cluster size.  This is the user defined
    quantity aggregationProb while aggregationProbVec identifies which
    receptors can undergo association, as updated above. Now passing both
    aggregationProb and aggregationProbVec. Use aggregationProbVec to
    prevent specific receptors from undergoing association. Note that
    giving an aggregatinProb value of 0.0 for a given cluster size will
    prevent that cluster and all larger clusters from being formed since
    aggregation happens one receptor at a time.
    %}
    
    try     
       %{
        %For associations via linear or arc collision probabilities, enable
        %the following. Must also enable and set associationDistVals above.
        %To not use either type of collision probability, call the original
        %function receptorAggregationAlg_maxWeightedMatching_noCollProb.
        %Alternatively, one can use the function below but set slope to 1
        %and intercept to 0 in the linear model.
        
        [cluster2receptor,receptor2cluster,clusterSize,positionsNew] = ...
            receptorAggregationAlg_maxWeightedMatching_CollProb(positionsNew,...
            associationDistVals,aggregationProbVec,aggregationProb,...
            collisionProbVals,receptor2cluster,cluster2receptor,clusterSize);
        %}
        
        numClustPre = length(cluster2receptor(:,1));
        
        [cluster2receptor,receptor2cluster,clusterSize,positionsNew,aggregationProbVec,sureAssocCount] = ...
            receptorAggregationAlg_maxWeightedMatching_sureColl(positionsNew,...
            aggregationDist,aggregationProbVec,aggregationProb,receptor2cluster,...
            cluster2receptor,clusterSize);        
        
        numClustPost = length(cluster2receptor(:,1));
        assocStats.numSureAssoc(iIter,1) = numClustPre - numClustPost;
        
        %assocStats.numSureAssoc(iIter,2) = nansum(sureAssocCount);
        assocStats.sureAssocCountBySize(1:length(sureAssocCount),iIter) = sureAssocCount;
        
        clear numClustPre numClustPost sureAssocCount
    
    catch newAssocFunExcep        
        fprintf('\nError at iIter = %d\n',iIter);
        disp(newAssocFunExcep.message);
        pause;
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %07/03/14 (ryirdaw)
    %Perform association of nodes not handled above and whose pairwise 
    %distance >= aggregationDist but < aggregationDistCorr    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %{
    %The following version of _potColl uses distance values calculated at
    %each displacement time points. The one below uses the approach where a
    %minimum for pairwise distances within the particles step is searched by
    %looking at the x and y components of the steps against time.
    %Calculation of individual dispalcement values was found to be time
    %consuming - the latter approach much faster.
    
    [aggregationProbVec,receptor2cluster,cluster2receptor,clusterSize,...
        positionsNew] = receptorAggregationAlg_maxWeightedMatching_potColl_070114(...
        cluster2receptor,receptor2cluster,positionsOld,positionsNew,...
        aggregationProbVec,aggregationProb,associationDistVals,...
        numDispPerUnitTime,dispTimeVec,clusterSize);     
    %}
    
    [~,receptor2cluster,cluster2receptor,clusterSize,...
        positionsNew,numPotColl,numColl,numPotColl_Assoc,numCollProbPairs,...
        collProbStatsTemp] =...
        receptorAggregationAlg_maxWeightedMatching_potColl(...
        cluster2receptor,receptor2cluster,positionsOld,positionsNew,...
        aggregationProbVec,aggregationProb,associationDistVals,clusterSize);
    
    assocStats.numPotColl(iIter) = numPotColl;
    assocStats.numColl(iIter) = numColl;
    assocStats.numPotColl_Assoc(iIter) = numPotColl_Assoc;
    assocStats.numCollProbPairs(iIter) = numCollProbPairs;

    collProbStats(iIter) = collProbStatsTemp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    [numClusters,maxClustSize] = size(cluster2receptor);
    
    %store new receptor information
    receptorTraj(:,:,iIter) = positionsNew;
    recept2clustAssign(:,iIter) = receptor2cluster;
    clust2receptAssign(1:numClusters,1:maxClustSize,iIter) = cluster2receptor;

    progressText((iIter-1)/(numIterations-1),'Simulation');
        
end %(for iIter = 2 : numIterations)

%% Post-processing

%remove the initialization period from simulation
receptorTraj = receptorTraj(:,:,numIterInit+1:end);
recept2clustAssign = recept2clustAssign(:,numIterInit+1:end);
clust2receptAssign = clust2receptAssign(:,:,numIterInit+1:end);

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

%KJ (150528): call function to label and sub-sample
receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
    labelRatio,intensityQuantum);


%% ~~~ the end ~~~


%% subfunction 1

%Function no longer used
%Outsourced to:
%receptorAggregationAlg_maxWeightedMatching_sureColl
%and
%receptorAggregationAlg_maxWeightedMatching_potColl

function [cluster2receptor,receptor2cluster,clusterSize,receptPositions] ...
    = receptorAggregationAlg(receptPositions,aggregationDist,aggregationProb,receptor2clusterPrev)

%get number of receptors
numReceptors = size(receptPositions,1);

%determine whether each receptor would aggregate if given the chance
receptorAggregFlag = rand(numReceptors,1) < aggregationProb;

%calculate inter-receptor distances
interReceptDist = pdist(receptPositions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%07/05/13 (ryirdaw)
%If avoiding aggregation of clusters, use interReceptDist determined above
%and call removeClusterAggregation along with previous receptor2cluster and
%association distance.
%Only necessary if there are two or more clusters currently. Determine the
%number of clusters without passing additional vars, using for loops or
%time consuming functions like find,hist, etc... - i.e.
%This is one compact option but takes time
%if ( sum(hist(receptor2clusterPrev,max(receptor2clusterPrev)) > 1) > 1 )  
%The following determines the presence of clusters without the total number
%of clusters - faster
if (sum(interReceptDist == 0) > 1)
    interReceptDist = removeClusterAggregation(interReceptDist,receptor2clusterPrev,aggregationDist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        

%construct the hierarchical tree of receptors
receptTree = linkage(interReceptDist);

%construct clusters from the hierarchical tree
%receptor2cluster indicates the cluster to which each receptor belongs
receptor2cluster = cluster(receptTree,'cutoff',aggregationDist,'criterion','distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%07/05/13 (ryirdaw)
%The cluster function above changes the cluster labels for receptors even
%if there has not been any change in the clusters from last iteration.
%Keep the same cluster labels as previous iteration if there has not been
%an aggregation event.
receptor2cluster = fixClusterLabels(receptor2cluster,receptor2clusterPrev);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %09/03/13 (ryirdaw)
        %The current cluster may have more than one receptor aggregating.
        %In this case, pick only one of these to aggregate. OR, a
        %dissociation event may have occured just prior to call to this
        %function, within the same iteration.
        %Needed only after initial step, for sizes larger than dimers and
        %when there is a change from previous iteration step. In the event
        %that a dissociation has occured, all receptors involved will have
        %flag = 0. The third condition for the if handles the case where a
        %dissociation has occured within this iteration and the leftover
        %cluster is a dimer - withouth the condition
        %the dimer will be forced to dissociate because numMembers == 2 and
        %so the function will not be invoked.
        %NOTE: keep the above aggregFlag line since the if right below
        %the call to singlerReceptorClustering may need it.
        if (~isempty(receptor2clusterPrev) &&...
                ((numMembers > 2)  ||  all(aggregFlag == 0)) )
            aggregFlag = singleReceptorClustering(aggregFlag,clusterMembers,receptor2clusterPrev);                       
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

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

%% subfunction 2 - still used

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine whether each receptor would dissociate if given the chance
%receptorDissociateFlag = rand(numReceptors,1) < dissociationProb;

%find clusters with more than one receptor
clustersBig = find(clusterSize > 1);

%06/27/13 (1 of 2)
%{
dissociationProb used on clusters instead of receptors (above).
Probability can be determined in two ways:  a single value for all current
clusters or each cluster has its own value. For a dissociating cluster, a
receptor will be randomly picked and removed from the cluster (below).
Dissociation happens one receptor at a time. 
%}
%Each cluster has its own value
clusterDissociateFlag = rand(numClusters,1) < dissociationProb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go over these clusters
numClustersInit = numClusters; %this is for storing receptors that dissociate
for iCluster = clustersBig'

    %get receptors belonging to this cluster
    clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get their dissociation flags
    %dissociateFlag = receptorDissociateFlag(clusterMembers);
    
    %06/27/13 (2 of 2)
    %For a cluster that is dissociating, determine which receptor to
    %dissociate by randomly picking one of the receptors.
    dissociateFlag = zeros(numel(clusterMembers),1);
    if (clusterDissociateFlag(iCluster))
        %Current cluster is dissociating. Pick a random integer between 1
        %and number of receptors in cluster.
        recept2Dissociate = randi(numel(clusterMembers),1);
        %Set the flag for picked receptor to 1.
        dissociateFlag(recept2Dissociate) = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %07/17/13
    %dissociateFlag = dissociateCluster(dissociationProb,clusterMembers);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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


%% OLD CODE WHEN LABELING AND SUBSAMPLING WAS DONE INSIDE FUNCTION

%got moved to outside function "genReceptorInfoLabeled" on 28 May 2015

% %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
% numLabelRatio = length(labelRatio);
% receptorInfoLabeled(numLabelRatio,1) = struct('receptorTraj',[],...
%     'recept2clustAssign',[],...
%     'clust2receptAssign',[],...
%     'compTracks',[],...
%     'labelRatio',[]);
%
% for lRindx=1:numLabelRatio
%
%     %if labeling ratio is less than one ...
%     if labelRatio(lRindx) < 1
%
%         %label some receptors based on the labeling ratio
%         %1 = labeled, 0 = unlabeled
%         %090514 - modified to allow a vector labelRatio
%         labelFlag = rand(numReceptors,1) <= labelRatio(lRindx);
%         indxLabeled = find(labelFlag==1);
%         indxNotLabeled = find(labelFlag==0);
%
%         %extract the trajectories of the labeled receptors
%         receptorTrajLabeled = receptorTraj(indxLabeled,:,:);
%
%         %assign the intensities of the labeled receptors
%         receptorIntensityLabeled = intensityQuantum(1) + ...
%             randn(length(indxLabeled),numIterSim) * intensityQuantum(2);
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %12/06/13 (ryirdaw)
%         %Reset intensity values that are < epsilon to epsilon
%         receptorIntensityLabeled(receptorIntensityLabeled < eps) = eps;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         %extract the cluster assignments of the labeled receptors
%         recept2clustAssignLabeled = recept2clustAssign(indxLabeled,:);
%
%         %modify the cluster-to-receptor assignments to include only labeled
%         %receptors
%         clust2receptAssignLabeled = clust2receptAssign;
%         %convert receptors that are not labeled to zero
%         for iNotLabeled = indxNotLabeled'
%             clust2receptAssignLabeled(clust2receptAssignLabeled==iNotLabeled) = 0;
%         end
%         %update the indices of the labeled receptors
%         for iLabeled = 1 : length(indxLabeled)
%             clust2receptAssignLabeled(clust2receptAssignLabeled==indxLabeled(iLabeled)) = iLabeled;
%         end
%         %make sure that zeros come after receptor indices
%         clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2,'descend');
%
%         %remove empty clusters (which include unlabeled receptors)
%         %modify receptor-to-cluster assignments accordingly
%         for iIter = 1 : numIterSim
%             clustSize = sum(clust2receptAssignLabeled(:,:,iIter)~=0,2);
%             indxFull = find(clustSize~=0);
%             indxEmpty = find(clustSize==0);
%             clust2receptAssignLabeled(:,:,iIter) = clust2receptAssignLabeled(...
%                 [indxFull;indxEmpty],:,iIter);
%             for iFull = 1 : length(indxFull)
%                 recept2clustAssignLabeled(recept2clustAssignLabeled(:,iIter)...
%                     ==indxFull(iFull),iIter) = iFull;
%             end
%         end
%
%         %remove empty rows and columns from clust2receptAssign
%         cluster2receptor = max(clust2receptAssignLabeled,[],3);
%         columnSum = sum(cluster2receptor);
%         clust2receptAssignLabeled = clust2receptAssignLabeled(:,columnSum~=0,:);
%         rowSum = sum(cluster2receptor,2);
%         clust2receptAssignLabeled = clust2receptAssignLabeled(rowSum~=0,:,:);
%
%         convTracksStartTime = tic;
%
%         %put labeled receptor trajectories and clusters into the format of the
%         %output of trackCloseGapsKalman
%         compTracksLabeled = convReceptClust2CompTracks(clust2receptAssignLabeled,...
%             recept2clustAssignLabeled,receptorTrajLabeled,receptorIntensityLabeled);
%
%         convTracksETime = toc(convTracksStartTime);
%         fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);
%
%         %put information in receptorInfoLabeled
%         %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
%         %{
%         %original
%         receptorInfoLabeled = struct('receptorTraj',receptorTrajLabeled,...
%             'recept2clustAssign',recept2clustAssignLabeled,...
%             'clust2receptAssign',clust2receptAssignLabeled,...
%             'compTracks',compTracksLabeled);
%         %}
%         receptorInfoLabeled(lRindx).receptorTraj = receptorTrajLabeled;
%         receptorInfoLabeled(lRindx).recept2clustAssign = recept2clustAssignLabeled;
%         receptorInfoLabeled(lRindx).clust2receptAssign = clust2receptAssignLabeled;
%         receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;
%         receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
%
%
%     else %if all receptors are labeled
%
%         %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
%         %{
%         %original
%         %copy receptor information from all to labaled
%         receptorInfoLabeled = receptorInfoAll;
%         %}
%
%         receptorInfoLabeled(lRindx).receptorTraj = receptorInfoAll.receptorTraj;
%         receptorInfoLabeled(lRindx).recept2clustAssign = receptorInfoAll.recept2clustAssign;
%         receptorInfoLabeled(lRindx).clust2receptAssign = receptorInfoAll.clust2receptAssign;
%         receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
%
%         %assign receptor intensities
%         receptorIntensity = intensityQuantum(1) + ...
%             randn(numReceptors,numIterSim) * intensityQuantum(2);
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %12/06/13 (ryirdaw)
%         %Reset intensity values that are < epsilon to epsilon
%         receptorIntensity(receptorIntensity < eps) = eps;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         convTracksStartTime = tic;
%
%         %put labeled receptor trajectories and clusters into the format of the
%         %output of trackCloseGapsKalman
%         compTracksLabeled = convReceptClust2CompTracks(clust2receptAssign,...
%             recept2clustAssign,receptorTraj,receptorIntensity);
%
%         convTracksETime = toc(convTracksStartTime);
%         fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);
%
%         %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
%         %{
%         %original
%         %store the new compound tracks with the labeling intensity
%         receptorInfoLabeled.compTracks = compTracksLabeled;
%         %}
%         receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;
%
%     end %(if labelRatio < 1 ... else ...)
%
% end %for each labelRatio
