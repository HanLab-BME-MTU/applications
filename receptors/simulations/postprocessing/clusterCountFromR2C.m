function clusterInfo = clusterCountFromR2C(recept2clustAssign)
%
%
    %Get total number of simulations
    numSims = length(recept2clustAssign);
    
    %Get total number of receptors and iterations - assuming all
    %recept2clustAssign are for repeated simulations
    [numReceptors,numIters] = size(recept2clustAssign{1});
    
    clusterCount = NaN(numReceptors,numIters,numSims);
    clusterCountMean = [];
    largestClustSize = 0;
    
    for simIndx=1:numSims
        currR2C = recept2clustAssign{simIndx};
        
        for iterIndx=1:numIters             
            for clustID=1:max(currR2C(:,iterIndx))
                currClustSize = sum(currR2C(:,iterIndx) == clustID);
                if (currClustSize > 0)
                    clusterCount(currClustSize,iterIndx,simIndx) =...
                        nansum([clusterCount(currClustSize,iterIndx,simIndx) 1]);
                    
                    if (currClustSize > largestClustSize)
                        largestClustSize = currClustSize;
                    end
                end
            end
            
        end        
       
        clear currR2C        
    end
        
    %calculate mean if more than one simulation
    if (numSims > 1)
        clusterCountMean = NaN(largestClustSize,1);
        
        for iterIndx=1:numIters
            for sizeIndx=1:largestClustSize
                clusterCountMean(sizeIndx,iterIndx) =...
                    nanmean(clusterCount(sizeIndx,iterIndx,:));
            end
        end
        
    end
    
    %Trim extra rows
    clusterCount(largestClustSize+1:end,:,:) = [];
    
    %Save values for return
    clusterInfo.clusterCount = clusterCount;
    clusterInfo.clusterCountMean = clusterCountMean;
    clusterInfo.largestClustSize = largestClustSize;
    
            
        
    