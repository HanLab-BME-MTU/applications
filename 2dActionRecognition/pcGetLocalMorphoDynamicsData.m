%% Returns a list of size 1 x N for each measure, were N is the total accumulated number of cells in all experiments in all time points
function [allSizes,allEccentricities,allSpeeds,allMatchingScores,allLabels, out] = ...
    pcGetLocalMorphoDynamicsData(allExps, labels)

nExps = length(allExps);

allSizes = [];
allEccentricities = [];
allSpeeds = [];
allMatchingScores = [];
allLabels = [];
for e = 1 : nExps
    curExp = allExps{e};
    curLabel = labels(e);
    
    fname = [curExp.dir '/results/' curExp.name '_localMorphDynam.mat'];
    
    load(fname); % 'sizes','eccentricities','speeds','matchingScores'
    n = size(sizes);
    
    allSizes = [allSizes sizes];
    allEccentricities = [allEccentricities eccentricities];
    allSpeeds = [allSpeeds speeds];
    allMatchingScores = [allMatchingScores matchingScores];
    allLabels = [allLabels; curLabel.*ones(n,1)];    
end

params.percentile = 2;
params.nBins = 20;
[out.sizes] = histogramQuantification(allSizes,allLabels,params); % out.measure.hists, out.measure.diffHists

[out.eccentricities] = histogramQuantification(allEccentricities,allLabels,params);
    
[out.speeds] = histogramQuantification(allSpeeds,allLabels,params);
    
[out.matchingScores] = histogramQuantification(allMatchingScores,allLabels,params);

end