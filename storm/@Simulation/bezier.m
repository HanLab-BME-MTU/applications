function t = bezier(obj,controlPoints,nSamples)
if strcmp(obj.samplingType,'random')
    % Random sampling
    t = rand(nSamples,1);
    points = renderBezier(controlPoints,arcLengthToNativeBezierParametrization(controlPoints,t));
elseif strcmp(obj.samplingType,'regular')
    % Regular sampling
    t = linspace(0,1,nSamples)';
    points = renderBezier(controlPoints,arcLengthToNativeBezierParametrization(controlPoints,t));
end
obj.data.points = [obj.data.points;points];
end