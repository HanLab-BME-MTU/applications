function [receptor2cluster,cluster2receptor,clusterSize] = updateR2C2RwithNodeInfo(nodeList)
%UPDATER2C2RWITHNODEINFO creates receptor2cluster and cluster2receptor
%matrices from list of nodes.
%
%INPUT:     nodeList:   a 1D struct array of nodes.
%
%OUTPUT:    receptor2cluster:   a 2D array of receptor-to-cluster
%                               relationships same as in main simulation.
%                               Rows are receptor IDs and columns are
%                               iterations with entries indicating cluster.
%
%           cluster2receptor:   a 2D array of cluster-to-receptor
%                               relationships same as in main simulations.
%                               Rows are cluster IDs, columns are receptors
%                               with entries giving receptor IDs.
%
%   Robel Yirdaw, 07/18/14
%


    %receptor2cluster should have number of elements equal to number of
    %receptors which can be summed from the node sizes
    receptor2cluster = NaN(nansum([nodeList.size]),1);
    %cluster2receptor should have the same rows as number of nodes and
    %columns equal to the largest node size. Since nodeList contains blank
    %nodes from association events possibly, the length is longer than
    %actual number of nodes. These will be trimmed at the end
    cluster2receptor = zeros(length(nodeList),max([nodeList.size]));
    %Similarly for clusterSize
    clusterSize = zeros(length(nodeList),1);
    
    %Count the actual number of nodes present (minus blanks)
    nodeCount = 1;
    
    nodeIndx = 1;
    
    %Begin processing each element in nodeList
    while(nodeIndx <= length(nodeList))
        %Nodes with empty members are skipped
        if (~isempty(nodeList(nodeIndx).members))
            %Set cluster size
            clusterSize(nodeCount) = nodeList(nodeIndx).size;
            %Assign receptor IDs (members,columns) to the current cluster (rows)
            cluster2receptor(nodeCount,1:nodeList(nodeIndx).size) = ...
                nodeList(nodeIndx).members;
            %Set the cluster ID of the current member receptors. The
            %cluster ID is simply the node count.
            receptor2cluster(nodeList(nodeIndx).members,1) = nodeCount;            
            %Increment
            nodeCount = nodeCount + 1;
        end

        nodeIndx = nodeIndx + 1;                    
    end %while

    %Trim empty array elements                
    cluster2receptor(nodeCount:end,:) = [];
    clusterSize(nodeCount:end,:) = [];
    
end %function

