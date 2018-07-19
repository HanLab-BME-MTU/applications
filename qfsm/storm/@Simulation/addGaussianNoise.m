function addGaussianNoise(obj,sigmaNoise)
if numel(sigmaNoise) == 1
    noiseGaussian = randn(obj.data.nPoints,3)*sigmaNoise;
elseif numel(sigmaNoise) == 3
    noiseGaussian = randn(obj.data.nPoints,3).*repmat(sigmaNoise,obj.data.nPoints,1);
end

obj.data.points = obj.data.points + noiseGaussian;
disp('Simulation: Gaussian position-noise added!');

end