function [ cc ] = ccAddBackground( cc, edges )
%ccAddBackground Adds a background connected component to a connected
%components structure

if(isstruct(edges))
    edges = logical(labelmatrix(edges));
end

bg = ~imfill(edges,'holes');
cc = connectedComponents.ccAppend(cc,find(bg));


end

