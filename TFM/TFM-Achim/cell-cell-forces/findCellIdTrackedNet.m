function cellIdTracked=findCellIdTrackedNet(trackedNet,cellPos)
minD=Inf;
numNodes=length(trackedNet.node);
for iNode=1:numNodes
    [comp,diffVec]=compPts(cellPos,trackedNet.node{iNode}.pos);
    currD=sqrt(sum(diffVec(:).^2,1));
    if comp
        % we found exact match:
        cellIdTracked=iNode;
        return;
    elseif currD<minD
        % we take the best match if no perfect match is found:
        cellIdTracked=iNode;   
        minD=currD;
    end
end
display(['no perfect match achieved, go with min. deviation= ',num2str(minD)])
    