function [edge minD]=findNearestEdge(edgeStruc,network)

minD=Inf;

for edgeNum=1:length(network.edge)
    if ~isempty(network.edge{edgeNum})
        % Distance to strPt:
        D1a=sqrt(sum((network.edge{edgeNum}.strPt-edgeStruc.strPt).^2));
        D1b=sqrt(sum((network.edge{edgeNum}.strPt-edgeStruc.endPt).^2));

        % Distance to endPt:
        D2a=sqrt(sum((network.edge{edgeNum}.endPt-edgeStruc.strPt).^2));
        D2b=sqrt(sum((network.edge{edgeNum}.endPt-edgeStruc.endPt).^2));
    else
        % You can never link to a non-exisiting edge:
        D1a=Inf;
        D1b=Inf;
        D2a=Inf;
        D2b=Inf;
    end
    
    % Take the pairing with the minimal additive sum:
    D12a=D1a+D2b;
    D12b=D1b+D2a;
    
    if D12a<=D12b
        D12=D12a;
    else
        D12=D12b;
    end
    
    if ~isempty(network.edge{edgeNum})
        % Take the center to center distance as an additional criteria:
        D3=sqrt(sum((mean(network.edge{edgeNum}.intf,1)-mean(edgeStruc.intf,1)).^2));
    else
        % You can never link to a non-exisiting edge:
        D3=Inf;
    end
        
    
    D=D12+D3;
    
    if D<minD
        edge=edgeNum;
        minD=D;
    end
end