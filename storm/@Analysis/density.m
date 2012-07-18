function density(obj)
nBins = round(sqrt(obj.data.nClusters));
pointDensity = obj.data.clusterSize./obj.data.modelLength;

% Weight the individual densities with the rounded length
histInput = arrayfun(@repmat,pointDensity,round(obj.data.modelLength),ones(size(obj.data.modelLength)),'UniformOutput',false);
% histInput = arrayfun(@repmat,pointDensity,1+0*round(obj.data.modelLength),ones(size(obj.data.modelLength)),'UniformOutput',false);
histInput = vertcat(histInput{:});

nBins = nBins * 4; % 2,4,8,16
binCenters = 0.01:0.01:0.85;

% [~,xout] = hist(histInput,nBins);
[~,xout] = hist(histInput,binCenters);
binSize = xout(2)-xout(1)

% Display the weighted point density histogram
% hist(histInput,nBins);
hist(histInput,binCenters);
% xlim([0,0.9]);
xlabel('Point density [1/nm]');
ylabel('Amount of filament length [nm]');
end