function edgesHistoryLast(obj,idx)

if nargin == 2
    pointsStart = obj.data.points(obj.data.edgesHistory{end-idx+1}(:,1),:);
    pointsEnd = obj.data.points(obj.data.edgesHistory{end-idx+1}(:,2),:);
else
    pointsStart = obj.data.points(obj.data.edgesHistory{end}(:,1),:);
    pointsEnd = obj.data.points(obj.data.edgesHistory{end}(:,2),:);
end

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Edges History');

end
