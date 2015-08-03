function interReceptDist = removeClusterAggregation(interReceptDist,receptor2clusterPrev,aggregationDist)
%REMOVECLUSTERAGGREGATION   Check for and resolve cluster aggregation. 
%   At any given iteration point aggregation can occur involving the 
%   following
%       1. A pair of receptors
%       2. A receptor and a cluster
%       3. A pair of clusters
%   All of the above involving multiple entities can also occur - i.e.
%   multiple pairs of receptors aggregating, or multiple receptors
%   aggregating with a cluster or multiple clusters, or multiple clusters
%   aggregating with other clusters.
%
%   The goal of this function is to isolate and deal with case #3 without
%   affecting cases 1 and 2.  This is done by using the current
%   inter-receptor distances and identifying the clustered receptors
%   attempting to cluster.  In that case, the distance is increased beyond
%   the threshold in order to prevent the aggregation.  All other events
%   occur as before. The modified inter-receptor distance is returned.
%   Modifications on 09/10/13 handle the "bridging" case which was not
%   handled in the older (original) version from 07/05/13.  This scenario
%   involves two clusters and a free receptor where the clusters are beyond
%   the aggregationDist but the free receptor is within the aggregationDist
%   to both clusters.  In this case, the free receptor acts as a bridge and
%   results in the two clusters being aggregated along with the receptor.
%   This "bridging" effect can also invovle more than one free receptors
%   linking two receptors that are each within aggregationDist of two
%   clusters. This case is now handled here by preventing any
%   receptor/cluster to have only one interaction with another
%   receptor/cluster.
%
%   NOTE: the issue of dissociation followed by association within the same
%   iteration point, which results in the positions being updated but not
%   receptor2cluster, requires a specific check below. The case where a
%   receptor leaves its cluster (dissociation) and tries to aggregate back
%   to its original cluster is handled; receptor leaving its cluster and
%   trying to aggregate to another cluster not handled (not enough info
%   available in the passed args.).  
%   As of 09/12/13, dissociatons followed by associations described above
%   are being handled in the main code.  No changes made to this code.
%
%   INPUT:  interReceptDist:        inter-receptor distances from pdist
%           receptor2clusterPrev:   receptor2cluster array from previous
%                                   iteration with cluster assignments for
%                                   each receptor
%           aggregationDist:        user-defined receptor merge radius
%
%   OUTPUT: interReceptDist:        modified or original inter-receptor
%                                   distances
%
%   Robel Yirdaw, 07/05/13
%       Modified 09/12/13
%
%

    %Number of receptors
    numReceptors = length(receptor2clusterPrev);
    %First convert interReceptDist, which is a vector, into a matrix form
    interReceptDistMtrx = squareform(interReceptDist);
    %Counter for entries in interReceptDist (vector)
    vecCounter = 1;
    %Count number of clusters a receptor can aggregate with - column 1 is
    %the the count, column 2 is the cluster ID
    numClustInProximity = zeros(numReceptors,2);
    
    %Will traverse lower-diagonal of matrix - therefore no redundnacy
    for colIndx=1:(numReceptors-1)
        %The receptor corresponding to current column
        clustID = receptor2clusterPrev(colIndx,1);
        %Determine if cluster is shared - i.e. receptors in cluster
        recepInClust = receptor2clusterPrev(:,1)==clustID;
        %Iterate through all other receptors
        for rowIndx=(colIndx+1):numReceptors
            %Get the cluster id of the second receptor
            clustID2 = receptor2clusterPrev(rowIndx,1);
            %Determine if cluster is shared - i.e. receptors in cluster
            recepInClust2 = receptor2clusterPrev(:,1)==clustID2;
            
            %If a value in the square matrix is 0, then those receptors are
            %clustered already, skip.  Nonzero values are for *entities* not
            %clustered.
            if ( (clustID ~= clustID2) && (interReceptDistMtrx(rowIndx,colIndx) ~= 0) &&...
                    (interReceptDistMtrx(rowIndx,colIndx) < aggregationDist) )
                %Found new association event. Determine if it is cluster to
                %cluster, receptor to cluster or cluster to receptor...
                
                if ( (sum(recepInClust) > 1) && (sum(recepInClust2) > 1) )
                        %Found cluster pair that can aggregate.
                        interReceptDist(vecCounter) = 2*aggregationDist; 
                        
                else
                    %Now checking whether the free receptor can act as a
                    %"bridge" between two or more existing clusters.  Here,
                    %every particle (receptor or cluster) is allowed to
                    %have one interaction.  All others will be prevented by
                    %resetting the inter-receptor distance.
                    
                    %CASE1: The current *column* is the free receptor
                    if (numClustInProximity(colIndx,1) == 0)
                        numClustInProximity(colIndx,1) = numClustInProximity(colIndx,1) + 1;
                        numClustInProximity(colIndx,2) = clustID2;
                    elseif ( (numClustInProximity(colIndx,1) == 1) && ...
                            (numClustInProximity(colIndx,2) ~= clustID2) )
                        %Found "bridge" interaction (> 1 interaction)
                        interReceptDist(vecCounter) = 2*aggregationDist;                        
                    end

                    %CASE2: The current *column* is a receptor in cluster 
                    %clustID - the *row* is the free receptor interacting 
                    %with this cluster
                    if (numClustInProximity(rowIndx,1) == 0)
                        numClustInProximity(rowIndx,1) = numClustInProximity(rowIndx,1) + 1;
                        numClustInProximity(rowIndx,2) = clustID;
                    elseif ((numClustInProximity(rowIndx,1) == 1) && ...
                            (numClustInProximity(rowIndx,2) ~= clustID) )
                        %Found "bridge" interaction (> 1 interaction)
                        interReceptDist(vecCounter) = 2*aggregationDist;
                    end
                    
                end %If 2
                
            end %If 1

            %Update vector counter
            vecCounter = vecCounter + 1;

        end %For each row

    end %For each col           
    