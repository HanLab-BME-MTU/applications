function density(obj)
nBins = round(sqrt(obj.data.nClusters));
bigClusterIdx = 1:obj.data.nClusters;
len = obj.data.modelLength(bigClusterIdx);
nPoints = cellfun(@numel,obj.data.clusters(bigClusterIdx));
pointDensity = num2cell(nPoints./len);

% Weight the individual densities with the rounded length
histInput = cellfun(@repmat,pointDensity,num2cell(round(len)),num2cell(ones(size(len))),'UniformOutput',false);
histInput = vertcat(histInput{:});

% Display the weighted point density histogram
hist(histInput,2*nBins);
xlabel('Point density [1/nm]');
ylabel('Amount of filament length [nm]');
end