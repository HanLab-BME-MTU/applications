function [merge2splitTimes,numRecepAgg,numRecepDiss] = getMerge2SplitTimes_Long(recept2clustAssign)
%GETMERGE2SPLITTIMES_Long calculates receptor merge-to-split times from 
%recept2clustAssign. This is a longer version of getMerge2SplitTimes.
%
%   INPUT:  recept2clustAssign: A 2D array output by
%           receptorAggregationSimple as a struct field for either
%           receptorInfoAll and receptorInfoLabeled.
%
%   OUTPUT: 
%           merge2splitTimes: a 2D array of with the following column 
%           values
%              Col 1: Receptor ID
%              Col 2: Merge time
%              Col 3: Split time
%              Col 4: Merge-to-split time
%              Col 5: Cluster size (number of receptors) at dissociation
%
%           numRecepAgg:    number of receptors undergoing association
%           numRecepDiss:   number of receptors undergoing dissociation
%
%   Robel Yirdaw, 07/09/13
%

    %Determine number of receptors and number of iterations in simulation
    [numReceptors, numIters] = size(recept2clustAssign);    
    %Vector to track merge times for each receptor
    mergeTime = nan(numReceptors,1);
    %Receptor association matrix
    recepAssociation = zeros(numReceptors,numReceptors);
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
                    
                    numRecepAgg(iterIndx,2) = numRecepAgg(iterIndx,2) + 1;
                end
                %Tracks receptor associations, even if a merge is not
                %occuring. This will be used when splitting.
                recepAssociation(recepIndx,:) = recepInClust;
                recepAssociation(:,recepIndx) = recepInClust;                     
                
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
                        clustSizePrev = recept2clustAssign(recepIndx,iterIndx-1);
                        %Number of receptors (cluster size)
                        merge2splitTimes(m2sTimeCounter,5) = sum(recept2clustAssign(:,iterIndx-1)==clustSizePrev);
                        
                        m2sTimeCounter = m2sTimeCounter + 1;

                        %Next determine how to reset mergeTime and
                        %recepAssociation.  Two general cases: 1) if
                        %current receptor is splitting from one other
                        %receptor *and* the two merged at the same time
                        %point, then reset entries for both receptors. 2)
                        %for all other cases, reset entries only for the
                        %current receptor.                        
                        
                        %Get all receptors with the same merge time.
                        %It is possible that more than two receptors merged
                        %at a single time point.
                        allReceptors = find(mergeTime == mergeTime(recepIndx));
                        
                        %From the association matrix determine how many
                        %receptors will be left when the current one is
                        %removed.  Note: this is separate from what
                        %allReceptors shows. There could have been
                        %successive merge events after an initial two
                        %particle merge.
                        numRecepLeft = sum(recepAssociation(recepIndx,:) ~= 0) - 1;
                                                
                        %This will be true two particle splitting only if
                        %allReceptors has two entries and the receptor
                        %association shows 1 receptor will be left.  It is
                        %possible to have > 1 receptors left even if
                        %allReceptors has two entries in which case only
                        %the current receptor is removed and respective
                        %variables reset. In the earlier case, variables for
                        %both receptors will be reset.
                        if ( (numel(allReceptors) == 2) && (numRecepLeft == 1) )
                            %Must check if
                            %second receptor is associated with other
                            %receptors. This is an issue when merging of
                            %clusters occurs or when a receptor merges with
                            %an existing cluster and now only two receptors
                            %left but with different merge times.
                            %NOTE: integerating the if block below with the
                            %above will not work because the check on the
                            %second recptor can not be done above - there
                            %maybe more than two receptors with the same
                            %merge time values.
                            
                            otherRecep = allReceptors(allReceptors ~= recepIndx);
                            if (sum(recepAssociation(otherRecep,:)) == 2)
                                %Split occuring between *two* receptors that
                                %merged at the same time.                                
                                mergeTime(allReceptors) = 0;
                                %Reset receptor association
                                recepAssociation(allReceptors,:) = 0;                          
                                recepAssociation(:,allReceptors) = 0;   
                                
                                numRecepDiss(iterIndx,2) = numRecepDiss(iterIndx,2) + 2;
                            else
                                %Split occuring between *two* receptors that
                                %came together at different time points.                                                            
                                mergeTime(recepIndx) = 0;
                                %Reset receptor association
                                recepAssociation(recepIndx,:) = 0;                            
                                recepAssociation(:,recepIndx) = 0;    
                                
                                numRecepDiss(iterIndx,2) = numRecepDiss(iterIndx,2) + 1;
                            end
                        else
                            %Split not occuring between two receptors -
                            %more than two receptors with the same merge
                            %time and/or more than one receptor will be
                            %left after dissociation occurs.
                            mergeTime(recepIndx) = 0;
                            %Reset receptor association
                            recepAssociation(recepIndx,:) = 0;                            
                            recepAssociation(:,recepIndx) = 0;    
                            
                            numRecepDiss(iterIndx,2) = numRecepDiss(iterIndx,2) + 1;
                        end
                        
                    end
                end
                
            end %If receptor is part of cluster
            
        end %For all receptors
        
    end %For each iteration
    
    merge2splitTimes(m2sTimeCounter:end,:) = [];
    
end %function
                    
            
        
    
    
    