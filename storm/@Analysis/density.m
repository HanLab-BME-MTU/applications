function density(obj)
nBins = round(sqrt(obj.data.nClusters));
pointDensity = obj.data.clusterSize./obj.data.modelLength;

% Weight the individual densities with the rounded length
histInput = arrayfun(@repmat,pointDensity,round(obj.data.modelLength),ones(size(obj.data.modelLength)),'UniformOutput',false);
histInput = vertcat(histInput{:});

% Display the weighted point density histogram
hist(histInput,20*nBins);
xlabel('Point density [1/nm]');
ylabel('Amount of filament length [nm]');
end