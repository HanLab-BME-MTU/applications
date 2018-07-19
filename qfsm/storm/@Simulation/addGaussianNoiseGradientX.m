function addGaussianNoiseGradientX(obj,sigmaNoiseStart,sigmaNoiseEnd,xStart,xEnd)
sigmaNoise = (obj.data.points(:,1)-xStart)./(xEnd-xStart).*(sigmaNoiseEnd-sigmaNoiseStart)+sigmaNoiseStart;
obj.data.error = [sigmaNoise sigmaNoise sigmaNoise];

errors = repmat(sigmaNoise,1,3);

noiseGaussian = randn(obj.data.nPoints,3).*errors;
obj.data.points = obj.data.points + noiseGaussian;

obj.data.error = errors;

disp('Simulation: Gaussian position-noise gradient added!');
end