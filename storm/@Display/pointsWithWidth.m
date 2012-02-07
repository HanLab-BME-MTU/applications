function pointsWithWidth(obj)

% Normalize the width
pointRadii = obj.data.width/mean(obj.data.width)*obj.pointSize;

obj.imaris.displayPoints(obj.data.points,pointRadii',obj.pointColor,'Display: Points');

obj.imaris.fitCamera();

end

