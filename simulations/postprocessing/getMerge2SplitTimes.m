function [merge2splitTimes,numRecepAgg,numRecepDiss] = getMerge2SplitTimes(recept2clustAssign)
%GETMERGE2SPLITTIMES Calculate receptor merge-to-split times from 
%           recept2clustAssign.
%
%   INPUT:  recept2clustAssign: A 2D array output by
%           receptorAggregationSimple as a struct field for either
%           receptorInfoAll and receptorInfoLabeled.
%
%   OUTPUT: merge2splitTimes: a 2D array of with the following column 
%           values
%              Col 1: Receptor ID
%              Col 2: Merge time
%              Col 3: Split time
%              Col 4: Merge-to-split time
%              Col 5: Cluster size (number of receptors) at dissociation
%
%   Robel Yirdaw, 07/09/13
%

    %Determine number of receptors and number of iterations in simulation
    [numReceptors, numIters] = size(recept2clustAssign);    
    %Vector to track merge times for each receptor
    mergeTime = nan(numReceptors,1);
    %Initialize the return using the number of iterations in simulation
    merge2splitTimes = zeros(numIters,5);
    
    %Number of receptors undergoing association (merging) at each iteration
    numRecepAgg = zeros(numIters,2);
    %First column is simply the iteration number - populate
    numRecepAgg(:,1) = 1:numIters;
    %Number of receptors undergoing dissociation (splitting) at each iteration
    numRecepDiss = zeros(numIters,2);
    %First column is simply the iteration number - populate here
    numRecepDiss(:,1) = 1:numIters;

    
    %Counter
    m2sTimeCounter = 1;
    
    for iterIndx=1:numIters
        for recepIndx=1:numReceptors           
            
            %Get the cluster id of current receptor
            clustNum = recept2clustAssign(recepIndx,iterIndx);
            %Determine if cluster is shared - i.e. receptors in cluster
            recepInClust = recept2clustAssign(:,iterIndx)==clustNum;
                                  
            %Count number of receptors in current receptor's cluster
            numRecepInClust = sum(recepInClust);
            if (numRecepInClust > 1)
                %Current receptor is part of a cluster.  Check if merging.
                %Receptor starting simulation merged is not a merging
                %event.
                if (mergeTime(recepIndx) == 0)
                    %Receptor merging. Save merge time.
                    mergeTime(recepIndx) = iterIndx;    
                    %===
                    %Tracking associations per receptor
                    numRecepAgg(iterIndx,2) = numRecepAgg(iterIndx,2) + 1;
                    %===
                end                
                
            else
                %Current receptor not part of a cluster.  
                %Prep to track merge-to-split if necessary
                if (isnan(mergeTime(recepIndx)))
                    mergeTime(recepIndx) = 0;
                else
                    %Check if splitting.
                    if (mergeTime(recepIndx) > 0)                        
                        %Splitting.  Calculate merge-to-split time.
                        merge2splitTimes(m2sTimeCounter,1) = recepIndx;
                        merge2splitTimes(m2sTimeCounter,2) = mergeTime(recepIndx);
                        merge2splitTimes(m2sTimeCounter,3) = iterIndx;
                        merge2splitTimes(m2sTimeCounter,4) = iterIndx - mergeTime(recepIndx);
                        %Determine number of receptors in cluster from
                        %which current receptor is dissociating
                        %Get the cluster id
                        clustIDPrev = recept2clustAssign(recepIndx,iterIndx-1);
                        %Number of receptors (cluster size)
                        merge2splitTimes(m2sTimeCounter,5) = sum(recept2clustAssign(:,iterIndx-1)==clustIDPrev);
                        %Increment merge-to-split event counter
                        m2sTimeCounter = m2sTimeCounter + 1;
                        %Reset merge time
                        mergeTime(recepIndx) = 0;
                        %===
                        %Tracking dissociations per receptor
                        numRecepDiss(iterIndx,2) = numRecepDiss(iterIndx,2) + 1;
                        %===
                    end %if splitting
                end
                
            end %If receptor is part of cluster
            
        end %For all receptors
        
    end %For each iteration
    
    merge2splitTimes(m2sTimeCounter:end,:) = [];
    
    
                    
            
        
    
    
    