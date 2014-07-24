function initialEdges(obj)

pointsStart = obj.data.points(obj.data.initialEdges(:,1),:);
pointsEnd = obj.data.points(obj.data.initialEdges(:,2),:);

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Initial Edges');

end

