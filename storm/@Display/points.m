function points(obj,varargin)

if nargin == 2 % Display only a subset of points
    obj.imaris.displayPoints(obj.data.points(varargin{1},:),obj.pointSize,obj.pointColor,'Display: Points');
else % Display all the points
    obj.imaris.displayPoints(obj.data.points,obj.pointSize,obj.pointColor,'Display: Points');
end

obj.imaris.fitCamera();

end

