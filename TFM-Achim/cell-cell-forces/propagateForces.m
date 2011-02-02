function constrForceField=propagateForces(constrForceField,frame)
% This function has been cross validated with the old version
% propagateForcesOld. All good!
edge  =constrForceField{frame}.network.edge;
node  =constrForceField{frame}.network.node;
numNodes=length(node);
numEdges=length(edge);


% Run through all edges:
for k=1:numEdges
    % first, set all forces along the edges to []:
    edge{k}.f1 =[];
    edge{k}.f2 =[];
    
    % initialize the search:
    % the node sets are:
    n2visit(1).set=edge{k}.nodes(1);
    n2visit(2).set=edge{k}.nodes(2);
    
    nVisited(1).set=[];
    nVisited(2).set=[];
    
    % the edge sets are:
    e2visit(1).set=[];
    e2visit(2).set=[];
    
    eVisited(1).set=k;
    eVisited(2).set=k;
    
    for i=1:2
        while ~isempty(n2visit(i).set)
            currNode=n2visit(i).set(1);
            n2visit(i).set(1)=[];
            nVisited(i).set=union(nVisited(i).set,currNode);
            
            % new edges that have not been treated so far are:
            newEdges=setdiff(node{currNode}.edges,eVisited(i).set);
            
            % add those to the edge to-do-list:
            e2visit(i).set=union(e2visit(i).set,newEdges);
            
            while ~isempty(e2visit(i).set)
                currEdge=e2visit(i).set(1);
                e2visit(i).set(1)=[];
                eVisited(i).set=union(eVisited(i).set,currEdge);
                
                % new edges that have not been treated so far are:
                newNodes=setdiff(edge{currEdge}.nodes,nVisited(i).set);
                
                % add those to the edge to-do-list:
                n2visit(i).set=union(n2visit(i).set,newNodes);
            end
        end
    end
    
    % compare the two node lists:
    if isempty(intersect(nVisited(1).set,nVisited(2).set)) && length(union(nVisited(1).set,nVisited(2).set))==numNodes
        % Then forces can be calculated:
        f1=0;
        for nodeId=nVisited(1).set
            f1=f1+node{nodeId}.vec;
        end
        edge{k}.f1=f1;
        
        f2=0;
        for nodeId=nVisited(2).set
            f2=f2+node{nodeId}.vec;
        end
        edge{k}.f2=f2;
    elseif length(nVisited(1).set)==numNodes && length(nVisited(1).set)==numNodes
        display(['Edge: ',num2str(k),' is in a loop']);
    else
        display('Something went wrong')
    end
end

constrForceField{frame}.network.edge=edge;

% 
% 
% 
% 
% % run through all nodes and find the ones with connectivity=1.
% list_nodes_c1=findNodesC1(nodeRed);
% 
% while ~isempty(list_nodes_c1)
%     % take the first node:
%     currNode=list_nodes_c1(1);
%     % delete it from the list:
%     list_nodes_c1(1)=[];
%     
%     % calculate the communicated forces at this edge, that is the
%     % reduced force at this node (including already propagated forces)
%     f1=nodeRed{currNode}.vec;
%     nodeRed{currNode}.vec=[0 0];
%     
%     % and the sum over the resudual forces of all other cells (which
%     % are not yet zero):
%     f2=0;
%     for n=1:numNodes
%         f2=f2+nodeRed{n}.vec;
%     end
%     
%     % Now propagate this force to the neighboring node.
%     % get the edge number:
%     edgeNumb=nodeRed{currNode}.edges;
%     % find the number of the neighboring node:
%     neighNode = setdiff(edge{edgeNumb}.nodes, currNode);
%     % propagate the force:
%     nodeRed{neighNode}.vec=nodeRed{neighNode}.vec+f1;
%     
%     % assigne the forces to the edge in the real nodes!
%     [~, loc] = ismember(currNode, edge{edgeNumb}.nodes);
%     if loc==1
%         edge{edgeNumb}.f1=f1;
%         edge{edgeNumb}.f2=f2;
%     else
%         edge{edgeNumb}.f1=f2;
%         edge{edgeNumb}.f2=f1;
%     end
%     
%     % kill that edge at the two nodes and reduce the degree of both
%     % nodes:
%     nodeRed{currNode}.edges =setdiff(nodeRed{currNode} .edges,edgeNumb);
%     nodeRed{currNode}.deg   = nodeRed{currNode}.deg-1;
%     nodeRed{neighNode}.edges=setdiff(nodeRed{neighNode}.edges,edgeNumb);
%     nodeRed{neighNode}.deg  = nodeRed{neighNode}.deg-1;
%     
%     % There is no problem if the connected had only one edge left, too.
%     % Then, we had arrived at the final cell pair. In the next line,
%     % the list_nodes_c1 should be empty:
%     
%     % recalculate the list of nodes:
%     list_nodes_c1=findNodesC1(nodeRed);
% end
% 
% % prepare the output:
% constrForceField{frame}.network.edge=edge;
% 
%     function list_nodes_c1_nst=findNodesC1(node_nst)
%         list_nodes_c1_nst=[];
%         for n_nst=1:length(node_nst)
%             if node_nst{n_nst}.deg==1
%                 list_nodes_c1_nst=horzcat(list_nodes_c1_nst,n_nst);
%             end
%         end
%     end
% 
% end