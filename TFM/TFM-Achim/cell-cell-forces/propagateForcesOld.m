function constrForceField=propagateForcesOld(constrForceField,frame)
    edge  =constrForceField{frame}.network.edge;
    node  =constrForceField{frame}.network.node;
    numNodes=length(node);
    nodeRed=node;
    
    % go and set all forces along the edges to []:
    for k=1:length(edge)
        edge{k}.f1 =[];
        edge{k}.f2 =[];
    end

    % run through all nodes and find the ones with connectivity=1. 
    list_nodes_c1=findNodesC1(nodeRed);
    
    while ~isempty(list_nodes_c1)
        % take the first node:
        currNode=list_nodes_c1(1);
        % delete it from the list:
        list_nodes_c1(1)=[];
        
        % calculate the communicated forces at this edge, that is the
        % reduced force at this node (including already propagated forces)
        f1=nodeRed{currNode}.vec;
        nodeRed{currNode}.vec=[0 0];               
        
        % and the sum over the resudual forces of all other cells (which 
        % are not yet zero):
        f2=0;
        for n=1:numNodes
            f2=f2+nodeRed{n}.vec;
        end
                
        % Now propagate this force to the neighboring node.
        % get the edge number:
        edgeNumb=nodeRed{currNode}.edges;
        % find the number of the neighboring node:
        neighNode = setdiff(edge{edgeNumb}.nodes, currNode);
        % propagate the force:
        nodeRed{neighNode}.vec=nodeRed{neighNode}.vec+f1;
    
        % assigne the forces to the edge in the real nodes!
        [~, loc] = ismember(currNode, edge{edgeNumb}.nodes);
        if loc==1
            edge{edgeNumb}.f1=f1;
            edge{edgeNumb}.f2=f2;
        else
            edge{edgeNumb}.f1=f2;
            edge{edgeNumb}.f2=f1;
        end        
        
        % kill that edge at the two nodes and reduce the degree of both
        % nodes:
        nodeRed{currNode}.edges =setdiff(nodeRed{currNode} .edges,edgeNumb);
        nodeRed{currNode}.deg   = nodeRed{currNode}.deg-1;
        nodeRed{neighNode}.edges=setdiff(nodeRed{neighNode}.edges,edgeNumb);
        nodeRed{neighNode}.deg  = nodeRed{neighNode}.deg-1;        
        
        % There is no problem if the connected had only one edge left, too.
        % Then, we had arrived at the final cell pair. In the next line,
        % the list_nodes_c1 should be empty:

        % recalculate the list of nodes:
        list_nodes_c1=findNodesC1(nodeRed);
    end
    
    % prepare the output:
    constrForceField{frame}.network.edge=edge;

    function list_nodes_c1_nst=findNodesC1(node_nst)
        list_nodes_c1_nst=[];
        for n_nst=1:length(node_nst)
            if node_nst{n_nst}.deg==1
                list_nodes_c1_nst=horzcat(list_nodes_c1_nst,n_nst);
            end
        end
    end

end