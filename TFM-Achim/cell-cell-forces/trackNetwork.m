function [constrForceField]=trackNetwork(constrForceField)
toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'network')
        toDoList=horzcat(toDoList,frame);
        constrForceField{frame}.network_tracked=[];
    end
end

% initialize the tracked network:
constrForceField{toDoList(1)}.network_tracked=constrForceField{toDoList(1)}.network;
% initialize the numbers of nodes and edges:
maxNodeId=length(constrForceField{toDoList(1)}.network.node);
maxEdgeId=length(constrForceField{toDoList(1)}.network.edge);
    
%track the nodes first:
for k=2:length(toDoList)
    prevFrame=toDoList(k-1);
    currFrame=toDoList(k);
    
    prevNumNodes=length(constrForceField{prevFrame}.network.node);
    currNumNodes=length(constrForceField{currFrame}.network.node);
    
    prevNumEdges=length(constrForceField{prevFrame}.network.edge);
    currNumEdges=length(constrForceField{currFrame}.network.edge);
    
    %**********************************************************************
    % Nodes:
    %**********************************************************************
    
    for nodeNum=1:currNumNodes
        pos=constrForceField{currFrame}.network.node{nodeNum}.pos;
        
        % find the corresponding node number in the previous frame:
        [corrNodeNum minD]=findNearestNode(pos,constrForceField{prevFrame}.network_tracked);
        reMapNode(nodeNum,:)=[corrNodeNum minD];
    end
    
    % It might be that cells divide, thus there might be more cells in the
    % current frame. Find the best match and set the other values to NaN.    
    for nodeNum=1:prevNumNodes
        ind=find(reMapNode(:,1)==nodeNum);
        if length(ind)>1
            % Then the map is ambiguous. Remove the correct point from the
            % remove list: 
            
            [~,id]=min(reMapNode(ind,2));
            ind(id)=[];
            
            reMapNode(ind,1)=NaN;
            reMapNode(ind,2)=NaN;

        end
    end
    
    % It might be that a cell breaks apart from the cluster (death
    % migration). No other cell should be within a certain radius. Find
    % those positions and warn the user, since this could be a potential
    % error source if at the same time a new cell joins the cluster...
    
    % Fill in the new node numbers:
    for nodeId=1:currNumNodes
        if isnan(reMapNode(nodeId,1))
            maxNodeId=maxNodeId+1;
            reMapNode(nodeId,1)=maxNodeId;
        end
        %constrForceField{currFrame}.network_tracked.node{reMapNode(nodeId,1)}=constrForceField{currFrame}.network.node{nodeId};
    end
    % Make sure that the number of node entries is: maxNodeId.
    constrForceField{currFrame}.network_tracked.node{maxNodeId}=[];
    % Now fill in the tracked values:
    for nodeId=1:currNumNodes
        constrForceField{currFrame}.network_tracked.node{reMapNode(nodeId,1)}=constrForceField{currFrame}.network.node{nodeId};
    end
       
    
    %**********************************************************************
    % Edges:
    %**********************************************************************
    
    for edgeNum=1:currNumEdges
        edgeStruc=constrForceField{currFrame}.network.edge{edgeNum};
        
        % finde the corresponding node number in the previous frame:
        [corrEdgeNum minD]=findNearestEdge(edgeStruc,constrForceField{prevFrame}.network_tracked);
        reMapEdge(edgeNum,:)=[corrEdgeNum minD];
    end
    
    % It might be that cells build up more interfaces (cell division), thus
    % there might be more edges in the current frame. Find the best match
    % and set the other values to NaN.
    for edgeNum=1:prevNumEdges
        ind=find(reMapEdge(:,1)==edgeNum);
        if length(ind)>1
            % Then the map is ambiguous. Remove the correct point from the
            % remove list: 
            
            [~,id]=min(reMapEdge(ind,2));
            ind(id)=[];
            
            reMapEdge(ind,1)=NaN;
            reMapEdge(ind,2)=NaN;
        end
    end    
    % It might be that an interface is lost. No other edge should be within
    % a certain radius. Find those positions and warn the user, since this
    % could be a potential error source if at the same time a new edge is
    % built up in the cluster.    
    
    % Fill in the new Edge numbers:
    for edgeId=1:currNumEdges
        if isnan(reMapEdge(edgeId,1))
            maxEdgeId=maxEdgeId+1;
            reMapEdge(edgeId,1)=maxEdgeId;
        end
        % constrForceField{currFrame}.network_tracked.edge{reMapEdge(edgeId,1)}=constrForceField{currFrame}.network.edge{edgeId};
    end
    % Make sure that the number of edge entries is: maxEdgeId.
    constrForceField{currFrame}.network_tracked.edge{maxEdgeId}=[];
    % Now fill in the tracked values:
    for edgeId=1:currNumEdges
        constrForceField{currFrame}.network_tracked.edge{reMapEdge(edgeId,1)}=constrForceField{currFrame}.network.edge{edgeId};
    end        
    
    % run through all the edges and adapt the field edge.nodes
    % appropriately:    
    for edgeID=1:maxEdgeId
        if ~isempty(constrForceField{currFrame}.network_tracked.edge{edgeID})
            nodesID=constrForceField{currFrame}.network_tracked.edge{edgeID}.nodes;
            newNodes=reMapNode(nodesID,1);
            constrForceField{currFrame}.network_tracked.edge{edgeID}.nodes=newNodes(:)';
        end
    end
    
    % run through all the edges and adapt the direction of the normal
    % vector to each interface since the direction might have changed from
    % frame to frame:   
    for edgeID=1:maxEdgeId
        if ~isempty(constrForceField{currFrame}.network_tracked.edge{edgeID})
            curr_nVec=constrForceField{currFrame}.network_tracked.edge{edgeID}.nVec_internal;
            % Only if the edge has existed previously:
            if edgeID<=length(constrForceField{prevFrame}.network_tracked.edge)
                prev_nVec=constrForceField{prevFrame}.network_tracked.edge{edgeID}.nVec_internal;
                if dot(curr_nVec,prev_nVec)<0
                    display('Change the direction of the interface normal vector, this part of the code needs to be debugged:')
                    % Then the product is most likely close to -1, that is, the
                    % dircetion has been flipped from frame to frame. In this
                    % case we have to invert the direction of the normal vector
                    % as well as the order of the interface positions:
                    constrForceField{currFrame}.network_tracked.edge{edgeID}.intf=flipud(constrForceField{currFrame}.network_tracked.edge{edgeID}.intf);
                    constrForceField{currFrame}.network_tracked.edge{edgeID}.intf_internal=flipud(constrForceField{currFrame}.network_tracked.edge{edgeID}.intf_internal);
                    constrForceField{currFrame}.network_tracked.edge{edgeID}.nVec_internal=-constrForceField{currFrame}.network_tracked.edge{edgeID}.nVec_internal;
                end
            end
        end
    end
    
    % run through all the nodes and adapt the field node.edges
    % appropriately:    
    for nodeID=1:maxNodeId
        if ~isempty(constrForceField{currFrame}.network_tracked.node{nodeID})
            edgesID=constrForceField{currFrame}.network_tracked.node{nodeID}.edges;
            newEdges=reMapEdge(edgesID,1);
            constrForceField{currFrame}.network_tracked.node{nodeID}.edges=newEdges(:)';
        end
    end
    
    
    constrForceField{currFrame}.network_tracked.stats=constrForceField{currFrame}.network.stats;
    
    clear reMapNode
    clear reMapEdge
    
end

return;


% to test this code execute:
constrForceField=trackNetwork(constrForceField);

toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'network')
        toDoList=horzcat(toDoList,frame);
    end
end

for frame=toDoList
    figure(frame)
    subplot(1,2,1)
    plotCellNetwork(constrForceField{frame}.network)
    subplot(1,2,2)
    plotCellNetwork(constrForceField{frame}.network_tracked)
end
