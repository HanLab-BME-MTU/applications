function orientationWithVariableLength(obj,maxModelLength)

% Variable segment length
maxMagnitude = max(obj.data.magnitude);
pointsStart = obj.data.points+obj.data.orientation.*repmat(obj.data.magnitude,1,3)/maxMagnitude*maxModelLength;
pointsEnd = obj.data.points-obj.data.orientation.*repmat(obj.data.magnitude,1,3)/maxMagnitude*maxModelLength;

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Orientation');

end

