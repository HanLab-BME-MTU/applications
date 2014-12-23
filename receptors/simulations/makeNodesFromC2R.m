function nodesList = makeNodesFromC2R(cluster2receptor,receptPositions,receptAssocFlag)
%MAKENODESFROMC2R uses cluster2receptor information to construct nodes
%which are clusters of any size. A node is represented by a struct with
%fields as defined below.
%
%INPUT:     cluster2receptor: a table of cluster-receptor relationships
%                             as defined in the main simulation.
%           receptor2cluster: a vector of receptor-cluster relationships
%                             as defined in the main simulation.
%           receptAssocFlag:  flag indicating whether a given receptor can
%                             undergo associaton.
%   
%OUTPUT:    nodesList:  a struct array of nodes where a node is a struct 
%                       with the following fields:
%                       1) members:     receptors belonging to the node
%                       2) size:        number of receptors
%                       3) xPos:        x coordinate (in microns)
%                       4) yPos:        y coordinate (in microns)
%                       5) assocFlag:   a boolean indicating whether this
%                                       node can undergo association
%
%   Robel Yirdaw, 06/24/14
%

    %Define the basic node element
    nodeElement = struct('members',[],'size',NaN,'xPos',NaN,'yPos',NaN,'assocFlag',NaN);
    
    %Determine number of nodes that will be formed
    numNodes = sum(cluster2receptor(:,1) ~= 0);
    
    %Initialize list of nodes as a struct array of node elements
    nodesList(1:numNodes,1) = nodeElement;
    
    %Build node list
    for nodeIter=1:numNodes
        %Get member receptors - some columns can be 0
        tempNodeMembers = cluster2receptor(nodeIter,cluster2receptor(nodeIter,:) ~= 0);
        
        %Assign members
        nodesList(nodeIter).members = tempNodeMembers;
        
        %Assign size
        nodesList(nodeIter).size = length(tempNodeMembers);
        
        %Assign coordinates - all receptors in the same node (cluster) will
        %have the same coordinates.
        nodesList(nodeIter).xPos = receptPositions(tempNodeMembers(1,1),1);
        nodesList(nodeIter).yPos = receptPositions(tempNodeMembers(1,1),2);
        
        %Assign association flag - all receptors in the same node will have
        %the same flag
        nodesList(nodeIter).assocFlag = receptAssocFlag(tempNodeMembers(1,1));
        
        clear tempNodeMembers
        
    end
    
end %function


        
        
    