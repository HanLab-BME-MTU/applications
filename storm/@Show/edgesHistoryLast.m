function edgesHistoryLast(obj)

pointsStart = obj.data.points(obj.data.edgesHistory{end}(:,1),:);
pointsEnd = obj.data.points(obj.data.edgesHistory{end}(:,2),:);

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Edges History');

end
