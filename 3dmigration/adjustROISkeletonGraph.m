function skelGraph = adjustROISkeletonGraph(skelGraph,roiInf,pixXY,pixZ)
%Adjusts the xyz coordinates to agree with the new origin of the ROI.

newOrig = [roiInf.cropX(1) roiInf.cropY(1) roiInf.cropZ(1)]-1;%Offset used to adjust for new image origin location
newOrig(3) = newOrig(3) * pixZ/pixXY;   %Adjust z for symmetric voxels

adjFun = @(x)(x - repmat(newOrig,[size(x,1),1]));

if ~isempty(skelGraph.vertices)
    skelGraph.vertices = adjFun(skelGraph.vertices);
end
if ~isempty(skelGraph.edgePaths)
    skelGraph.edgePaths = cellfun(adjFun,skelGraph.edgePaths,'Unif',0);
end

