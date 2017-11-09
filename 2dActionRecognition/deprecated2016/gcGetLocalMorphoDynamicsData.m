%% Returns a list of size 1 x N for each measure, were N is the total accumulated number of cells in all experiments in all time points
function [allSizes,allEccentricities,allSpeeds,allMatchingScores,allLabels, out] = ...
    gcGetLocalMorphoDynamicsData(allExps, labels, uniqueLabelsStr, uniqueLabelsValues)

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
    
    if ~exist(fname,'file')
        continue;
    end
    
    load(fname); % 'sizes','eccentricities','speeds','matchingScores'
    n = size(accumulatedSizes,2);
    
    allSizes = [allSizes accumulatedSizes];
    allEccentricities = [allEccentricities accumulatedEccentricity];
    allSpeeds = [allSpeeds accumulatedLocalSpeed];
    allMatchingScores = [allMatchingScores accumulatedLocalMatchingScore];
    allLabels = [allLabels; curLabel.*ones(n,1)];    
end

params.percentile = 2;
params.nBins = 40;
[out.sizes] = histogramQuantification(allSizes,allLabels,uniqueLabelsStr,uniqueLabelsValues,params); % out.measure.hists, out.measure.diffHists

[out.eccentricities] = histogramQuantification(allEccentricities,allLabels,uniqueLabelsStr,uniqueLabelsValues,params);
    
[out.speeds] = histogramQuantification(allSpeeds,allLabels,uniqueLabelsStr,uniqueLabelsValues,params);
    
[out.matchingScores] = histogramQuantification(allMatchingScores,allLabels,uniqueLabelsStr,uniqueLabelsValues,params);

end