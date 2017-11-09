function plotNodes(nodeList)
%PLOTNODES creates a plot showing the initial and final positions of nodes
%in the given node list (nodeList).
%
%   This function was primarily used when developing the path-based
%   association approach to visualize various possible paths a pair of
%   nodes can take. Initial positions are indicated by open circles and
%   final positions by closed circles. Paths are indicated by a dotted
%   line. See examples in the test script at directory 07/07/14.
%
%   INPUT:
%           nodeList:   a 2x1 array of node structs with initial and final
%                       positions of the nodes, i.e.,
%                       [node1;node2] where node1 and node2 are struct with
%                       initial and final x & y coordinates as fields:
%                           xi,yi - initial x & y coordinates
%                           xf,yf - final x & y coordinates
%
%   OUTPUT:
%           none.
%
%   Robel Yirdaw, July, 2014.
%
    numNodes = length(nodeList);
    %Create figure, get axes, and set to overlay plots
    figure();
    currAx = gca();
    set(currAx,'NextPlot','add');
    
    for nodeIndx=1:numNodes
        tempNode = nodeList(nodeIndx);
        %Generate a random color for each node
        randColor = rand(3,1);
        %Plot
        tempH1 = plot([tempNode.xi;tempNode.xf],[tempNode.yi;tempNode.yf],':o');
                
        %set(tempH1,'DisplayName',num2str(sizeIndx));
        set(tempH1,'Color',[randColor(1),randColor(2),randColor(3)]);
        
        %Update final position plots
        plot(tempNode.xf,tempNode.yf,'Marker','o',...
            'MarkerFaceColor',[randColor(1),randColor(2),randColor(3)],...
            'MarkerEdgeColor',[randColor(1),randColor(2),randColor(3)]);
       
        clear tempNode
    end
    
end %function