function t = bezierWithOrthogonalGaussianNoise2D(obj,controlPoints,nSamples,sigmaNoise)

% Convert 2D control points to 3D control points
if size(controlPoints,2) == 2
    controlPoints(:,3) = zeros(size(controlPoints,1));
end

% Sample the curve
if strcmp(obj.samplingType,'random')
    % Random sampling
    t = rand(nSamples,1);
    points = renderBezier(controlPoints,arcLengthToNativeBezierParametrization(controlPoints,t));
elseif strcmp(obj.samplingType,'regular')
    % Regular sampling
    t = linspace(0,1,nSamples)';
    points = renderBezier(controlPoints,arcLengthToNativeBezierParametrization(controlPoints,t));
end

% Compute the vector orthogonal to the curve
[~,normTanVec] = tangentBezier(controlPoints,t);
normOrthoVec = [-normTanVec(:,2),normTanVec(:,1),normTanVec(:,3)];

% Project the noise on the orthogonal vector
noiseGaussian = randn(nSamples,1)*sigmaNoise;
noiseGaussian = repmat(noiseGaussian,1,3).*normOrthoVec;

% Add the noise to the points
points = points + noiseGaussian;
obj.data.points = [obj.data.points;points];

end


