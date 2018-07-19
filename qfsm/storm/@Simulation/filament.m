function filament(obj,controlPoints,nSamples,sigmaNoise,sourceSpacing,fractionOfDamagedSources)

len = lengthBezier(controlPoints);
nSources = ceil(len/sourceSpacing)+1;
nDamagedSources = round(nSources*fractionOfDamagedSources);

% Create the sources
sourcesS = linspace(0,1,nSources)';

% Remove the damaged sources
sourcesIdx = randperm(nSources);
intactSourcesIdx = sourcesIdx(1:nSources-nDamagedSources);
sourcesS = sourcesS(intactSourcesIdx);

% Randomly pick sources
low = 1; high = nSources-nDamagedSources;
sourcesIdx = round(low + (high-low) * rand(nSamples,1));
sourcesS = sourcesS(sourcesIdx);

points = renderBezier(controlPoints,arcLengthToNativeBezierParametrization(controlPoints,sourcesS));

% Add the noise to the points
if numel(sigmaNoise) == 1
    noiseGaussian = randn(nSamples,3)*sigmaNoise;
elseif numel(sigmaNoise) == 3
    noiseGaussian = randn(nSamples,3).*repmat(sigmaNoise,nSamples,1);
end
points = points + noiseGaussian;
obj.data.points = [obj.data.points;points];

end
