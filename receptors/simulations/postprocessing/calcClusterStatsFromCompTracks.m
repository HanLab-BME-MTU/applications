function clusterStats = calcClusterStatsFromCompTracks(compTracks,observeSideLen,probDim)
%CALCCLUSTERSTATSFROMCOMPTRACKS calculates count, fraction and density of 
%each cluster size from compTracks.
%
%   Input:
%       1) compTracks: a cell containing compTracksALT structs for each
%                      simulation. aggregState from the defaultFormatTracks
%                      field will be used.
%       2) observeSideLen:  side length values used for the simulation
%                           which can be a single value or a value per side.
%       3) probDim: dimension of the simulation space (normally 2).
%
%   Output:
%       clusterStats: a struct with the following fields
%           a) clusterCount (max cluster x num iters x num sims): number of 
%              clusters per size per iteration per simulation
%           b) clusterCountMean (max cluster x num iters): mean of the
%              above over the simulations
%           c) clusterFrac (max cluster x num iters x num sims): fraction
%              of clusters per size per iteration per simulation
%           d) clusterFracMean (max cluster x num iters): mean of the
%              above over the simulations
%           e) clusterDensity (max cluster x num iters x num sims): density
%              of clusters per size per iteration per simulation
%           f) clusterDensityMean (max cluster x num iters): mean of the
%              above over the simulations
%           g) receptorDensity (num sims x num iters): density of receptors
%              per iteration per simulation
%           h) receptorDensityMean (1 x num iters): mean of the above over
%              the simulations
%           i) largestClusterSize: the largest cluster size from all sims
%           j) observeSideLen: the (input) value used for the calculations
%           k) probDim: the (input) value used for the calculations
%           l) simArea: calculated area of the simulation space used when
%              determining densities.
%
%   Robel Yirdaw, 02/14/14
%       Modified, 05/19/14
%       Modified, 09/11/14
%

    clusterStats = [];
    
    if (nargin == 1)
        observeSideLen = [];
        probDim = [];
        simArea = [];
        clusterDensity = [];
    elseif (nargin == 2)
        %Assign default value when not provided
        probDim = 2;
    end   

    numSideLenVals = length(observeSideLen);
    %Determine simulation area if possible
    if (numSideLenVals == 1)
        simArea = observeSideLen ^ probDim;
    elseif (numSideLenVals > 1)
        simArea = prod(observeSideLen);
    end     

    %compTracks from different simulations save in compTracks (cell vector)
    numSims = length(compTracks);
    
    %Get number of iterations - assuming all simulations have the same
    %number of iterations (i.e., columns in aggregState)
    numIters = length(compTracks{1}.defaultFormatTracks(1).aggregState(1,:));
    
    %Set a large cluster size for initializing matrices
    LRGST_CLUST_SIZE = 1000;
    
    %Cluster counts and fractions
    clusterCount = NaN(LRGST_CLUST_SIZE,numIters,numSims);
    clusterFrac = NaN(LRGST_CLUST_SIZE,numIters,numSims);    
    %If simulation area is defined, also initialize density matrix
    if (~isempty(simArea))
        clusterDensity = NaN(LRGST_CLUST_SIZE,numIters,numSims);
    end
    
    %Will keep track of the largest cluster size over all simulations
    largestClustSize = 0;
    
    %Begin
    progressText(0,'Processing');
    
    for simIndx=1:numSims
        %Pull out the current aggregState matrix, merged over all tracks.
        %When saving, it must be in full matrix form - histc below will
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %111814
        %currAggregState = vertcat(compTracks{simIndx}.defaultFormatTracks.aggregState);

        [~,~,~,~,currAggregState] =...
         convStruct2MatIgnoreMS(compTracks{simIndx}.defaultFormatTracks,0);
        
        for iterIndx=1:numIters
            %Determine the largest cluster size for the current iteration
            currLargestClustSize = full(max(currAggregState(:,iterIndx)));
            
            %Use histc to count how many of each size from 1 to
            %currLargestClustSize are present and save. When passing a
            %column from aggregState, pass the column in form since histc
            %does not take a sparce vector.
            clusterCount(1:currLargestClustSize,iterIndx,simIndx) =...
                histc(full(currAggregState(:,iterIndx)),(1:currLargestClustSize));
            
            %Track largest cluster size
            largestClustSize = max(largestClustSize,currLargestClustSize);
           
            %clear currLargestClustSize
            
        end % for each iteration
        
        clear currAggregState
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate total number of clusters at each iteration
        tempTotClusters = nansum(clusterCount(:,:,simIndx));
        %Use the totals to get fractions of clusters    
        clusterFrac(1:largestClustSize,:,simIndx) = clusterCount(1:largestClustSize,:,simIndx)./...
            repmat(tempTotClusters,largestClustSize,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate density of clusters if possible
        if (~isempty(simArea))    
            clusterDensity(1:largestClustSize,:,simIndx) =...
                clusterCount(1:largestClustSize,:,simIndx)./repmat(simArea,largestClustSize,numIters);
            %091114
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate receptor density.
            %First count the number of receptors for each cluster observed,
            %at, each iteration.
            receptorCountAll = repmat((1:largestClustSize)',1,numIters).*...
                clusterCount(1:largestClustSize,:,simIndx);
            %Calculate density. Here, rows correspond to simulations, columns
            %are iterations as above.
            receptorDensity(simIndx,:) = nansum(receptorCountAll)./simArea;            
        end
        
        clear tempTotClusters
        
        progressText((simIndx)/(numSims),'Processing');
        
    end %for each simulation
        
    
    %Trim extra rows
    clusterCount(largestClustSize+1:end,:,:) = [];
    clusterFrac(largestClustSize+1:end,:,:) = [];
    if (~isempty(simArea))    
        clusterDensity(largestClustSize+1:end,:,:) = [];
    end
    
    %Adjust NaN entries
    clusterCount(isnan(clusterCount)) = 0;
    clusterFrac(isnan(clusterFrac)) = 0;
    clusterDensity(isnan(clusterDensity)) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate mean count if more than one simulation
    if (numSims > 1)
        %Mean over simualtions, at each iteration
        clusterCountMean = mean(clusterCount,3);
        clusterFracMean = mean(clusterFrac,3);
        clusterDensityMean = mean(clusterDensity,3);
        receptorDensityMean = mean(receptorDensity);        
    else
        clusterCountMean = [];
        clusterFracMean = [];
        clusterDensityMean = [];
        receptorDensityMean = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    
    %Save values for return
    clusterStats.clusterCount = clusterCount;
    clusterStats.clusterCountMean = clusterCountMean;
    clusterStats.clusterFrac = clusterFrac;    
    clusterStats.clusterFracMean = clusterFracMean;
    clusterStats.clusterDensity = clusterDensity;
    clusterStats.clusterDensityMean = clusterDensityMean;
    clusterStats.receptorDensity = receptorDensity;
    clusterStats.receptorDensityMean = receptorDensityMean;    
    
    clusterStats.largestClustSize = largestClustSize;    
    clusterStats.observeSideLen = observeSideLen;
    clusterStats.probDim = probDim;
    clusterStats.simArea = simArea;
    
    
end
