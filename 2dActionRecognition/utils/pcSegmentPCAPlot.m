

function [] = pcSegmentPCAPlot(distributionMaps,threshold)
sortedMap = sort(distributionMaps(:),'descend');
accumulatedProbabilities = cumsum(sortedMap);
ind = find(accumulatedProbabilities > threshold,1);
probTH = sortedMap(max(ind-1,1));

BIN = distributionMaps > probTH;

figure; imagesc(BIN);

end