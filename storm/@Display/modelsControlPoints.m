function modelsControlPoints(obj)

controlPoints = vertcat(obj.data.modelBezCP{:});

obj.imaris.displayPoints(controlPoints,obj.controlPointSize,obj.controlPointColor,'Display: Control Points');

end

