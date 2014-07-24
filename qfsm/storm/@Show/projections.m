function projections(obj)

points = [];
for c=1:obj.data.nClusters
    if obj.data.modelType(c) == 1;
        curvePoints = renderBezier(obj.data.modelBezCP{c},obj.data.modelProj{c});
        points(end+1:end+size(curvePoints,1),1:3) = curvePoints;
    else
        curvePoints = renderBezier(obj.data.modelBezCP{c},obj.data.modelProj{c});
        points(end+1:end+size(curvePoints,1),1:3) = curvePoints;
    end
    
end

obj.imaris.displayPoints(points,obj.projectionSize,obj.projectionColor,'Display: Projections');

end