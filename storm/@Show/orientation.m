function orientation(obj,modelLength)

% Fixed segment length
pointsStart = obj.data.points+obj.data.orientation*modelLength/2;
pointsEnd = obj.data.points-obj.data.orientation*modelLength/2;

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Orientation');

end

