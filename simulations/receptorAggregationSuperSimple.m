function [time2aggreg,fracAggreg] = receptorAggregationSuperSimple(modelParam,simParam)
%receptorAggregationSuperSimple simulates the aggregation of receptors diffusing in a preset area; code stops when first aggregation happens
%
%SYNOPSIS time2aggreg = receptorAggregationSuperSimple(modelParam,simParam)
%
%INPUT  modelParam: Structure with the fields:
%           diffCoef        : Diffusion coefficient (microns^2/s).
%           confDim         : Confinement dimension (microns). Either one
%                             value, used for all dimensions, or a row
%                             vector with a value for each dimension.
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
%       fracAggreg : Fraction of receptors that can potentially aggregate.
%
%Khuloud Jaqaman, November 2013

%% Output

time2aggreg = [];
errFlag = 0;

%% Input

%check if correct number of arguments were used when function was called
if nargin < 1
    disp('--receptorAggregationSuperSimple: Too few input arguments');
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

%calculate number of receptors
obsRegionSize = prod(areaSideLen);
numReceptors = round(obsRegionSize * receptorDensity);

%in case of confined diffusion, determine how many times simulation must be
%repeated in order to explore same number of receptors as in a free
%diffusion case
if any(~isnan(confDim))
    numRepConf = 2;
    for iDim = 1 : probDim
        gridDim{iDim} = 0:confDim(iDim):areaSideLen(iDim); %#ok<AGROW>
    end
else
    numRepConf = 1;
    receptInComp{1,1} = (1:numReceptors)';
    numReceptInComp = numReceptors;
    confDim = areaSideLen;
    indxRepConf = 1;
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

%default fraction of receptors that are able to aggregate
fracAggreg = 1;

% progressText(0,'Simulation');

for iRep = 1 : numRep
    
    %% receptor initial positions and clustering
    
    %initialize receptor positions
    initPositions = rand(numReceptors,probDim) .* repmat(areaSideLen,numReceptors,1);
    
    %check if any two receptors are within aggregationDist from each other, indicating clustering
    pairwiseDist = pdist(initPositions);
    minDist = min(pairwiseDist);
    
    %if receptors are already clustered, that's it, mission accomplished
    if minDist <= aggregationDist
        
        time2aggreg(iRep) = 0;
        
    else %otherwise
        
        %if receptors undergo confined diffusion ...
        if numRepConf > 1
            
            %put them in compartments and shift their initial positions so that
            %0 is the lower left corner of each compartment
            %let's just do 2D for now
            receptInComp = cell(length(gridDim{1})-1,length(gridDim{2})-1);
            numReceptInComp = ones(size(receptInComp));
            for iComp2 = 1 : length(gridDim{2})-1
                for iComp1 = 1 : length(gridDim{1})-1
                    receptInComp{iComp1,iComp2} = find( ...
                        initPositions(:,1)>=gridDim{1}(iComp1)&initPositions(:,1)<gridDim{1}(iComp1+1) & ...
                        initPositions(:,2)>=gridDim{2}(iComp2)&initPositions(:,2)<gridDim{2}(iComp2+1) );
                end
            end
            for iComp2 = 1 : length(gridDim{2})-1
                for iComp1 = 1 : length(gridDim{1})-1
                    indxIn = receptInComp{iComp1,iComp2};
                    initPositions(indxIn,:) = initPositions(indxIn,:) - repmat([gridDim{1}(iComp1) gridDim{2}(iComp2)],length(indxIn),1);
                    numReceptInComp(iComp1,iComp2) = length(indxIn);
                end
            end
            
            %find compartments that have more than 1 receptor
            indxRepConf = find(numReceptInComp>1);
            numRepConf = length(indxRepConf);
            
            %determine fraction of receptors able to aggregate
            fracAggreg = sum(numReceptInComp(indxRepConf))/numReceptors;
            
        end
        
        %initialize intermediate output vector
        time2aggregTmp = NaN(numRepConf,1);
        
        %if there is need to simulate
        for jRepConf = 1 : numRepConf
            
            iRepConf = indxRepConf(jRepConf);
            
            %% Trajectory generation
            
            %reserve memory for output vectors
            receptorTraj = zeros(numReceptInComp(iRepConf),probDim,numIterations);
            
            %store initial information
            indxIn = receptInComp{iRepConf};
            receptorTraj(:,:,1) = initPositions(indxIn,:);
            
            %iterate in time
            iIter = 1;
            minDist = 2*aggregationDist;
            while minDist > aggregationDist && iIter <= numIterations
                
                %increase iteration index by 1
                iIter = iIter + 1;
                
                %get receptor positions at previous time point
                positionsOld = receptorTraj(:,:,iIter-1);
                
                %generate receptor displacements
                receptorDisp = stepStd*randn(numReceptInComp(iRepConf),probDim);
                
                %calculate the new receptor positions at the current time point
                positionsNew = positionsOld + receptorDisp;
                
                %make sure that receptors stay inside the region of interest
                correctionBoundaryLow = min(positionsNew,0);
                positionsNew = positionsNew - 2 * correctionBoundaryLow;
                correctionBoundaryUp = max(positionsNew - repmat(confDim,numReceptInComp(iRepConf),1),0);
                positionsNew = positionsNew - 2 * correctionBoundaryUp;
                
                %check if any two receptors are within aggregationDist from each other, indicating clustering
                pairwiseDist = pdist(positionsNew);
                minDist = min(pairwiseDist);
                
                %store new receptor information
                receptorTraj(:,:,iIter) = positionsNew;
                
            end %(while ~clusterFlag && iIter <= numIterations)
            
            if minDist <= aggregationDist
                time2aggregTmp(jRepConf) = (iIter-1)*timeStep;
            end
            
        end %(for iRepConf = indxRepConf')
        
        %get time to aggregate for this repeat
        if numRepConf > 0
            time2aggreg(iRep) = min(time2aggregTmp);
        end
        
    end %(if minDist <= aggregationDist)
    
%     progressText((iRep-1)/(numRep-1),'Simulation');
    
end %(for iRep = 1 : numRep)

%% ~~~ the end ~~~
