function [cluster2receptor,receptor2cluster,clusterSize,receptPositions,...
    aggregationProbVec,sureAssocCount] = receptorAggregationAlg_maxWeightedMatching_sureColl(...
    receptPositionsPrev,aggregationDist,aggregationProbVec,aggregationProb,...
    receptor2clusterPrev,cluster2receptorPrev,clusterSizePrev)
%RECEPTAGGREGATIONALG_MAXWEIGHTEDMATCHING_SURECOLL performs receptor-to-receptor
%   and receptor-to-cluster associations for pairs whose pairwise distance
%   falls below the default value (aggregationDist). These are considered
%   sure collisions. This function replaces receptorAggregationAlg in the 
%   simulation.  maxWeightedMatching is used to perform associations, 
%   limiting all receptors/clusters to associate with only one other 
%   receptor. Cluster-to-cluster associations are not allowed.
%
%   INPUT:
%       receptPositonsPrev:     2D array of current receptor positions
%       aggregationDist:        the distance threshold for association
%       aggregationProbVec:     1D array of receptor aggregation flags
%       aggregationProb:        1D array giving domain from which 
%                               aggregation probabilities will be drawn,
%                               for each cluster size (formation). Element
%                               #1 is assumed to be a default value, used
%                               when aggregationProb is a scalar or when
%                               its size < a given cluster size.
%       receptor2clusterPrev:   1D array of current cluster labels for each
%                               receptor (rows = receptors)
%       cluster2receptorPrev:   2D array of clusters along the rows and
%                               member receptors along the columns
%       clusterSizePrev:        1D array of cluster sizes 
%
%   OUTPUT:
%       cluster2receptor:   2D array of clusters along the rows and
%                           member receptors along the columns
%       receptor2cluster:   1D array of cluster labels for each receptor
%       clusterSize:        1D array of cluster sizes
%       receptPositions:    2D array of recpetor positions 
%       aggregationProbVec: the updated array of receptor association flags
%       sureAssocCount:     number of associations that occured
%
%
%   Robel Yirdaw, 10/02/13
%       Modified, 01/28/14
%       Modified, 03/25/14
%       Modified, 07/03/14 (renamed sureColl)
%       Modified, 08/11/14 (sure assoc counts returned)
%

    %Total number of receptors
    numReceptors = length(aggregationProb);
    
    %Determine association flag for all receptors
    %01/28/14: moved below for nodes instead of receptors
    %aggregFlag = rand(numReceptors,1) < aggregationProb;
    
    %The new vectors initialized to previous
    receptor2cluster = receptor2clusterPrev;
    cluster2receptor = cluster2receptorPrev;
    receptPositions = receptPositionsPrev;
    clusterSize = clusterSizePrev;

    %081114
    sureAssocCount = NaN;
    
    %Current number of nodes
    numNodesInit = sum(cluster2receptorPrev(:,1) ~= 0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %01/28/14
    %Determine association flag for all nodes. A row in cluster2receptor is
    %a node. Using the aggregationProb value assigned to the first member
    %of the cluster (row) when generating the flag for the node. This is
    %necessary because a cluster or receptor may have aggregatinProb value
    %of 0 indicating that it shouldn't undergo association because it has
    %undergone dissociation in the current iteration. Also, setting
    %aggregationProb to 1 for all members of a cluster in the main
    %simulation has been disabled. Thus, the possible values found in
    %aggregationProb are 0 (for a receptor or a set of receptors forming a
    %cluster) or the input value; previously it also included the
    %value of 1 for all receptors of a cluster as well.
    %
    %nodeFlag = rand(numNodesInit,1) <...
    %    aggregationProb(cluster2receptorPrev(1:numNodesInit,1),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   03/25/14
    %   The function now inputs associationProb and associationProbVec. In
    %   the earlier version, associationProb was actually
    %   associationProbVec. Here, associationProb is a vector and gives the 
    %   probability of forming a cluster of a given size. On the other
    %   hand, associationProbVec is an association flag for each receptor.
    %
    %   Thus, we first set up nodeFlag using associationProbVec. Those with
    %   value of 0, will not undergo association, and this would be due to
    %   having participated in a dissociation event in the current
    %   iteration.  However, this quantity can be used for other purposes.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nodeFlag = aggregationProbVec(cluster2receptorPrev(1:numNodesInit,1));
    
    %If there are receptors that can associate, then continue.
    if (any(nodeFlag))
        %Current nodes
        nodesInit = cell(numNodesInit,1);
        %Node sizes
        nodeSizes = NaN(numNodesInit,1);        
        %Position vector for current nodes        
        nodePositionsInit = NaN(numNodesInit,2);
        
        %Association flag for nodes
        %01/28/14: moved above
        %nodeFlag = NaN(numNodesInit,1);
        
        %Current largest cluster size
        largestClusterInit = 0;
        
        %Build node information
        for nodeIter=1:numNodesInit
            %A row in cluster2receptor is a node
            nodeMembers = cluster2receptorPrev(nodeIter,cluster2receptorPrev(nodeIter,:) ~= 0);
            %Assign members
            nodesInit{nodeIter} = nodeMembers;
            %Determine size of node
            nodeSizes(nodeIter) = length(nodeMembers);            
            %Get positions
            nodePositionsInit(nodeIter,1:2) = receptPositionsPrev(nodeMembers(1,1),1:2);
            
            %Get flag - assuming all members of a cluster have the same
            %flag (1)
            %01/28/14: moved above for nodes instead of receptors
            %nodeFlag(nodeIter,1) = aggregFlag(nodeMembers(1,1),1);
            
            %Save largest cluster size
            if (length(nodeMembers) > largestClusterInit)
                largestClusterInit = length(nodeMembers);
            end
            
            clear nodeMembers
        end %for all initial nodes
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get id of those that can associate (flag = 1)        
        aggregNodeID_aggregFlag = find(nodeFlag == 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get pairwise distances for those receptors with flag 1
        aggregNodePWDist = squareform(pdist(nodePositionsInit(aggregNodeID_aggregFlag,:)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %In the pairwise distance matrix, find the nodes where distance
        %criteria is satisfied. The nodes are the columns in the matrix and
        %what is stored in aggreNodeIndx_aggregDist is the column index.
        %These are not (necessarily) node ids.
        aggregNodeIndx_aggregDist = find(any((aggregNodePWDist < aggregationDist) & (aggregNodePWDist ~= 0)));        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        if (~isempty(aggregNodeIndx_aggregDist))
            %Nodes associating
            aggregNodeID = aggregNodeID_aggregFlag(aggregNodeIndx_aggregDist);
            %Count of nodes associating
            numNodesAggreg = length(aggregNodeID);
            
            %The number of pairwise associations
            numNodePairsWithinDist = 0.5*(sum(sum((aggregNodePWDist < aggregationDist) & (aggregNodePWDist ~= 0)) ));            

            %Collect edges and weights.
            paramEdges = NaN(numNodePairsWithinDist,2);
            paramWeights = NaN(numNodePairsWithinDist,1);
            edgeCount = 0;
            
            %Iterate through each
            for aggregNodeIter=1:numNodesAggreg                
                %Get a receptor from list of aggregating nodes. This will
                %be used to access the columns of the pair-wise distance
                %matrix below.
                currNodeID = aggregNodeID(aggregNodeIter);

                %Get a list of potential partners
                partnerNodeIndx = find((aggregNodePWDist(:,aggregNodeIndx_aggregDist(aggregNodeIter)) < aggregationDist) & ...
                    (aggregNodePWDist(:,aggregNodeIndx_aggregDist(aggregNodeIter)) ~= 0) );                
                partnerNodeIndx(partnerNodeIndx < aggregNodeIndx_aggregDist(aggregNodeIter)) = [];
                %The node ids corresponding to the index values above -
                %note this must be taken from aggregNodeID_aggregFlag, not
                %aggregNodeID.
                partnerNodeID = aggregNodeID_aggregFlag(partnerNodeIndx);               
                
                %Cluster associations not allowed
                if ( (nodeSizes(currNodeID) > 1) &&...
                        any(nodeSizes(partnerNodeID) > 1) )
                    %Get flags of nodes in partnerNodeID with sizes > 1
                    cluster2clusterAssoc = (nodeSizes(partnerNodeID) > 1);
                    %Remove values for those nodes
                    partnerNodeID(cluster2clusterAssoc) = [];
                    partnerNodeIndx(cluster2clusterAssoc) = [];
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %03/25/14
                %Determine the sizes of clusters that can potentially be
                %formed and use these to draw probabilities from the given
                %range, in aggregationProb. If aggregationProb is a scalar,
                %or there are more sizes than values in aggregationProb,
                %the vector will be expanded for the maximum needed size 
                %and the default domain assigned for these sizes.
                %
                if (~isempty(partnerNodeID))
                    %Determine the cluster sizes that can be formed
                    finalSizes = nodeSizes(currNodeID) + nodeSizes(partnerNodeID);
                    %If all sizes are not accounted for in aggregationProb,
                    %expand it and assign default value which is element 1.
                    if (length(aggregationProb) < max(finalSizes))
                        aggregationProb(length(aggregationProb)+1:max(finalSizes),1) =...
                            aggregationProb(1,1);
                    end
                    %Draw a random value for each cluster to be formed,
                    %using its size, and get association flag.
                    eventFlag = (rand(length(partnerNodeID),1) < aggregationProb(finalSizes) );
                    %Remove from the list partner nodes not associating
                    partnerNodeID(eventFlag == 0) = [];
                    partnerNodeIndx(eventFlag == 0) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    
                %Number of possible partners for the current node
                numPartners = length(partnerNodeID);
                %Determine edges and weights of current node with all
                %potential partners
                if (numPartners > 0)
                    %Determine edges. Note that partnerNodeID can have more
                    %than one value.
                    tempEdges(1:numPartners,1) = currNodeID;
                    tempEdges(1:numPartners,2) =  partnerNodeID;
                    %Determine weights
                    tempWeights = 1./aggregNodePWDist(partnerNodeIndx,aggregNodeIndx_aggregDist(aggregNodeIter));

                    %Accumulate
                    paramEdges(edgeCount+1:(edgeCount+numPartners),1:2) = tempEdges;
                    paramWeights(edgeCount+1:(edgeCount+numPartners)) = tempWeights;
                    
                    %Increment total count of edges (pairs of nodes)
                    edgeCount = edgeCount + numPartners;
                     
                    clear tempEdges tempWeights partnerNodeID partnerNodeIndx numPartners
                end

            end %for each associating node

            %Next, depending on the number of edges determined, perform
            %selection - if more than one edge found, use
            %maxWeightedMatching to pick edges.  
            selectedEdges = [];
            if (edgeCount == 1)
                %Trim empty rows in edges only
                paramEdges(isnan(paramEdges(:,1)),:) = [];             
                selectedEdges = paramEdges;
            elseif (edgeCount > 1)
                %Trim blank rows in edges and weights
                paramEdges(isnan(paramEdges(:,1)),:) = [];
                paramWeights(isnan(paramWeights(:,1)),:) = [];
                %The node input parameter is the largest node number
                paramNumNodes = max(paramEdges(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Peform weighted matching to pick a single pairwise
                %interaction for every receptor                
                selectionFlags = maxWeightedMatching(paramNumNodes,paramEdges,paramWeights);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Per returned flags above, pick the edges
                selectedEdges = paramEdges(selectionFlags,:);
            end
               
            %If edges selected, then associations have occured. Construct
            %the new positions, cluster to receptor and receptor to cluster
            %relationships.
            if (~isempty(selectedEdges))               
                %The new largest cluster size
                largestClusterNew = largestClusterInit;
                %Copy initial nodes to be modified
                nodesNew = nodesInit;
                
                %081114
                %Initialize sure association count array. Rows are new node
                %sizes and the largest possible is the current largest + 1.
                sureAssocCount = zeros(max(nodeSizes) + 1,1);
                sureAssocCount(1) = NaN;
                
                %Build new nodes list
                for newClusterIter=1:length(selectedEdges(:,1))
                    %Get node pair from list
                    tempNode = selectedEdges(newClusterIter,1:2);
                    %Insert new node
                    nodesNew{tempNode(1)} = [nodesNew{tempNode(1)} nodesNew{tempNode(2)}];
                    %The merged node will be set to NaN
                    nodesNew{tempNode(2)}(:) = NaN;
                    
                    %Update positions - receptors in nodes that associated
                    %will get the average position of the two nodes
                    %Average of node positions
                    %newPos = mean(nodePositionsInit(tempNode(1:2),:));
                    %receptPositions(nodesNew{tempNode(1)},1) = newPos(1);
                    %receptPositions(nodesNew{tempNode(1)},2) = newPos(2);
                    
                    %Average of receptor positions
                    %newPos = mean(receptPositionsPrev([nodesInit{tempNode(1:2)}],:));
                    %receptPositions(nodesNew{tempNode(1)},1) = newPos(1);
                    %receptPositions(nodesNew{tempNode(1)},2) = newPos(2);
                    %Or equivalently, but with one less indexing than above
                    newPosX = sum(nodePositionsInit(tempNode(1:2),1).*nodeSizes(tempNode(1:2)))/sum(nodeSizes(tempNode(1:2)));
                    newPosY = sum(nodePositionsInit(tempNode(1:2),2).*nodeSizes(tempNode(1:2)))/sum(nodeSizes(tempNode(1:2)));
                    receptPositions(nodesNew{tempNode(1)},1) = newPosX;
                    receptPositions(nodesNew{tempNode(1)},2) = newPosY;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %070314
                    %Update aggregationProbVec since it is now returned and
                    %will be used by the second association function, for
                    %nodes whose pairwise distance falls between
                    %aggregationDist and aggregationDistCorr, i.e.,
                    %aggregationDist <= pairwise distance <
                    %aggreagtionDistCorr. Setting the flag to 0 here will
                    %prevent those that associated here from associating in
                    %the second function.
                    aggregationProbVec(nodesNew{tempNode(1)}) = 0;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %The new largest cluster size
                    if (length(nodesNew{tempNode(1)}) > largestClusterNew)
                        largestClusterNew = length(nodesNew{tempNode(1)});
                    end
                    
                    %081114
                    %Accumulate sure association counts.
                    currNewNodeSize = length(nodesNew{tempNode(1)});
                    sureAssocCount(currNewNodeSize) =...
                        sureAssocCount(currNewNodeSize) + 1;
                    
                    clear tempNode currNewNodeSize
                end %for each new cluster
                
                %Build new receptor2cluster, cluster2receptor and
                %clusterSize variables.  Initialize cluster2receptor with
                %zeros not NaN - unused columns of cluster rows are
                %supposed to have zero values, not NaN.
                %Note: combining the while loop below with the for above
                %is not a good choice specifically because of
                %receptor2cluster. The construction of cluster2receptor can
                %be handled in the for above easily but not
                %receptor2cluster.  A change to a row in c2r means the
                %cluster label for those receptors has to be changed.
                %However, the label for all receptors in clusters
                %subsequent to the changed one will also have to be
                %updated - the cluster labels are in ascending order.
                receptor2cluster = NaN(numReceptors,1);
                cluster2receptor = zeros(numNodesInit,largestClusterNew);
                clusterSize = NaN(numNodesInit,1);
                numNodesNew = 1;
                initNodesIter = 1;
                while (initNodesIter <= numNodesInit)
                    %NaN rows in the new set of nodes will be skipped
                    if (~isnan(nodesNew{initNodesIter}(:)))
                        clusterSize(numNodesNew) = length(nodesNew{initNodesIter});
                        cluster2receptor(numNodesNew,1:length(nodesNew{initNodesIter})) = ...
                            nodesNew{initNodesIter};
                        receptor2cluster(nodesNew{initNodesIter},1) = numNodesNew;
                        numNodesNew = numNodesNew + 1;
                    end
                    
                    initNodesIter = initNodesIter + 1;                    
                end %while
                
                %Trim empty array elements                
                cluster2receptor(numNodesNew:end,:) = [];
                clusterSize(numNodesNew:end,:) = [];
                
            end %if associations occured
                
        end %if any distance within aggregationDist exists in distance mtx
        
    end %if any node has association flag = 1
        
    
end %function
