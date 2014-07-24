function statistics(obj)
totalLength = sum(obj.data.modelLength);
nNonNoisePoints = sum(cellfun(@numel,obj.data.clusters));
nNoisePoints = obj.data.nPoints - nNonNoisePoints;
noisePercent = nNoisePoints/nNonNoisePoints*100;
pointDensity = nNonNoisePoints/totalLength;
averagePointSpacing = 1/pointDensity;
averageClusterSize = mean(cellfun(@numel,obj.data.clusters));
averageModelLength = mean(obj.data.modelLength);

% Display statistics
disp('================ STATISTICS ================');
fprintf('              Total length: %.2f nm\n',totalLength);
fprintf('Number of non-noise points: %d points\n',nNonNoisePoints);
fprintf('    Number of noise points: %d points\n',nNoisePoints);
fprintf('              Noise points: %.2f %%\n',noisePercent);
fprintf('           Average spacing: %.2f nm\n',averagePointSpacing);
fprintf('        Number of clusters: %d clusters\n',obj.data.nClusters);
fprintf('      Average cluster size: %.2f points\n',averageClusterSize);
fprintf('      Average model length: %.2f nm\n',averageModelLength);
disp('============================================');
end