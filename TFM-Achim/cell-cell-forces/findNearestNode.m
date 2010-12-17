function [node minD]=findNearestNode(pos,network)

minD=Inf;

for nodeNum=1:length(network.node)
    if ~isempty(network.node{nodeNum})
        D=sqrt(sum((network.node{nodeNum}.pos-pos).^2));
    else
        % You can neve link to a non-existing node:
        D=Inf;
    end
    if D<minD
        node=nodeNum;
        minD=D;
    end
end