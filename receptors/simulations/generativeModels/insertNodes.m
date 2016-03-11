function [updatedNodeList,positionsNew] = insertNodes(nodePairsToMerge,currNodeList,positionsNew)
%INSERTNODES takes a set of indices corresponding to merging node pairs in
%the current node list, merges the pairs to form a set of single nodes and
%inserts these new nodes into the exsiting node list.
%
%INPUT:     nodePairsToMerge:   a 2D array of index values corresponding to
%                               the pairs of nodes to be merged. Number of
%                               elements equals the total number of merges
%                               to be performed (new nodes to be formed)
%                               and columns hold indices of the two nodes
%                               to be merged.
%
%           currNodeList:       a 1D struct array of nodes.
%
%           positionsNew:       2D coordinates of receptors compoing the
%                               nodes to be merged.
%
%OUTPUT:    updatedNodeList:    currNodeList updated with the associaitons.
%
%           positionsNew:       updated 2D coordinates of receptors after
%                               the merges.
%
%   Robel Yirdaw, 07/18/14
%


    %Initialize updated list
    updatedNodeList = currNodeList;
    
    %Build new nodes list
    for indx=1:length(nodePairsToMerge(:,1))
        
        %Get node pair ID
        tempNode = nodePairsToMerge(indx,1:2);
        
        %Merge node 2 with node 1 and insert new node as a new node 1
        updatedNodeList(tempNode(1)).members =...
            [currNodeList(tempNode(1)).members currNodeList(tempNode(2)).members];
        %Update its size
        updatedNodeList(tempNode(1)).size = length(updatedNodeList(tempNode(1)).members);
        
        %Update positions - receptors in nodes that associated
        %will get the average position of the two nodes
        newPosX = sum([currNodeList(tempNode(1:2)).xPos].*...
            [currNodeList(tempNode(1:2)).size])/sum([currNodeList(tempNode(1:2)).size]);
        newPosY = sum([currNodeList(tempNode(1:2)).yPos].*...
            [currNodeList(tempNode(1:2)).size])/sum([currNodeList(tempNode(1:2)).size]);
        
        %Update the new node list
        updatedNodeList(tempNode(1)).xPos = newPosX;
        updatedNodeList(tempNode(1)).yPos = newPosY;
        
        %Update the position vector with these values
        positionsNew(updatedNodeList(tempNode(1)).members,1) = newPosX;
        positionsNew(updatedNodeList(tempNode(1)).members,2) = newPosY;


        %The merged node, node 2, will be reset
        updatedNodeList(tempNode(2)).members = [];
        updatedNodeList(tempNode(2)).size = NaN;
        updatedNodeList(tempNode(2)).xPos = NaN;
        updatedNodeList(tempNode(2)).yPos = NaN;
        updatedNodeList(tempNode(2)).assocFlag = NaN;

        
        clear tempNode newPosX newPosY
        
    end %for each new cluster
    
    
end %function
