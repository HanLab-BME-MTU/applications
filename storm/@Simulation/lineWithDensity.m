function lineWithDensity(obj,startPoint,endPoint,density)
obj.line(startPoint,endPoint,ceil(density*norm(startPoint-endPoint)));
end