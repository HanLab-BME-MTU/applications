function clusterStats = calcClusterStatsFromR2C(recept2clustAssign,observeSideLen,probDim)
%CALCCLUSTERSTATSFROMR2C calculates count, fraction and density of each
%cluster size from recept2clustAssign.
%
%   Input:
%       1) recept2clustAssign: matrix with receptor-to-cluster
%                              relationships as ouptut from 
%                              receptorAggregationSimple.
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
    
    %Get total number of simulations
    numSims = length(recept2clustAssign);
    
    %Get total number of receptors and iterations - assuming all
    %recept2clustAssign are for repeated simulations
    [numReceptors,numIters] = size(recept2clustAssign{1});
    
    largestClustSize = 0;
    clusterCount = NaN(numReceptors,numIters,numSims);
    clusterFrac = NaN(numReceptors,numIters,numSims);
    
    if (~isempty(simArea))
        clusterDensity = NaN(numReceptors,numIters,numSims);
    end
    
    progressText(0,'Processing');
    
    for simIndx=1:numSims
        %051914 - removed currR2C for memory optimization
        for iterIndx=1:numIters    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %051914 - a function call to histc instead of a third for-loop
            %to tally up the occurence of each cluster size is
            %significantly faster
            
            %Create a vector with all possible ID values
            clustIDvec = ( (1:max(recept2clustAssign{simIndx}(:,iterIndx)))' );
            
            %Use histc to get number of receptors with each ID
            currClustSizeVec = histc(recept2clustAssign{simIndx}(:,iterIndx),clustIDvec);
            
            %The current largest size
            currLargestClustSize = max(currClustSizeVec);
            %Use histc again to get count of each cluster size            
            tempClusterCount = histc(currClustSizeVec,1:currLargestClustSize);
            
            %Save
            clusterCount(1:currLargestClustSize,iterIndx,simIndx) = tempClusterCount;
                        
            %Track largest cluster size
            largestClustSize = max(largestClustSize,currLargestClustSize);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clear clustIDvec currClustSizeVec clustID currClustSize 
            clear currLargestClustSize tempClusterCount
            
        end %for each iteration
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate total number of clusters at each iteration
        tempTotClusters = nansum(clusterCount(:,:,simIndx));
        %Use the totals to get fractions of clusters    
        clusterFrac(1:largestClustSize,:,simIndx) = clusterCount(1:largestClustSize,:,simIndx)./...
            repmat(tempTotClusters,largestClustSize,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate density of clusters and receptors, if possible
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
    
            
        
    