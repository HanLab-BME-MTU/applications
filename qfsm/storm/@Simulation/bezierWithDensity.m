function t = bezierWithDensity(obj,controlPoints,density)
nSamples = ceil(lengthBezier(controlPoints)*density);
t = obj.bezier(controlPoints,nSamples);
end