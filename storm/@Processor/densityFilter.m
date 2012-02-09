function densityFilter(obj,nNeighborsThreshold,ballRadius)
indices = KDTreeBallQuery(obj.data.points,obj.data.points,repmat(ballRadius,obj.data.nPoints,1));
nNeighbors = cellfun(@numel,indices);
obj.data.points = obj.data.points(nNeighbors >= nNeighborsThreshold,:);
obj.data.intensity = obj.data.intensity(nNeighbors >= nNeighborsThreshold);
obj.data.frame = obj.data.frame(nNeighbors >= nNeighborsThreshold);

disp('Process: Points in low density areas removed!');
end

