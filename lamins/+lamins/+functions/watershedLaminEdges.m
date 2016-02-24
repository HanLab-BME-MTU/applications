function [ w_masked ] = watershedLaminEdges( S )
%watershedLaminEdges Draw watershed boundaries around lamin edges

    neighbors = S.countNeighbors;
    justEdges = neighbors == 2 | neighbors == 1;
    distMap = bwdist(justEdges);
    w = watershed(distMap);
    w_masked = w;
    w_masked(~S.getMask) = 0;
end

