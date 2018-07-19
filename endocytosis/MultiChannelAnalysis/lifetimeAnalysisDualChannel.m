function [resCargo resNoCargo resGlobal] = lifetimeAnalysisDualChannel(data, kWeibull, cutoffIdx)

if nargin < 2
    kWeibull = [2 2 1];
end
if nargin < 3
    cutoffIdx = 1;
end


nExp = length(data);
for k = 1:nExp

    load([data(k).source 'LifetimeInfo' filesep 'lftInfo.mat']);
    lftMat = full(lftInfo.Mat_lifetime);
    statMat =  full(lftInfo.Mat_status);
        
    cargoIdx = find(data(k).cargoStatus == 1 & data(k).csValid == 1);
    nocargoIdx = find(data(k).cargoStatus == 0 & data(k).csValid == 1);
    
    cargoLifetimeVect = NaN(1, length(cargoIdx));
    nocargoLifetimeVect = NaN(1, length(nocargoIdx));
    
    lifetimeVect = NaN(size(lftMat,1),1);
    for t = 1:size(lftMat,1)
        cstat = nonzeros(statMat(t,:));
        if ( (min(cstat)==1) && (max(cstat)<5) )
            lifetimeVect(t) = max(lftMat(t,:));
        end
    end
    lifetimeVect = lifetimeVect(isfinite(lifetimeVect));
    data(k).lftHist = hist(lifetimeVect, 1:data(k).movieLength);
    
    for t = 1:length(cargoIdx)
        cstat = nonzeros(statMat(cargoIdx(t),:));
        if ( (min(cstat)==1) && (max(cstat)<5) )
            cargoLifetimeVect(t) = max(lftMat(cargoIdx(t),:));
        end
    end
    cargoLifetimeVect = cargoLifetimeVect(isfinite(cargoLifetimeVect));
    
    for t = 1:length(nocargoLifetimeVect)
        cstat = nonzeros(statMat(nocargoIdx(t),:));
        if ( (min(cstat)==1) && (max(cstat)<5) )
            nocargoLifetimeVect(t) = max(lftMat(nocargoIdx(t),:));
        end
    end
    nocargoLifetimeVect = nocargoLifetimeVect(isfinite(nocargoLifetimeVect));
    
    data(k).cargoLftHist = hist(cargoLifetimeVect, 1:data(k).movieLength);
    data(k).nocargoLftHist = hist(nocargoLifetimeVect, 1:data(k).movieLength);    
end

% create mean histograms
lftHistMean = stageFitLifetimes(data, 'lftHist');
cargoLftHistMean = stageFitLifetimes(data, 'cargoLftHist', lftHistMean.detectionCutoff);
nocargoLftHistMean = stageFitLifetimes(data, 'nocargoLftHist', lftHistMean.detectionCutoff);


resGlobal = lifetimeAnalysis(data, lftHistMean, kWeibull, cutoffIdx, 'all pits');
resCargo = lifetimeAnalysis(data, cargoLftHistMean, kWeibull, cutoffIdx, 'w/ cargo');
resNoCargo = lifetimeAnalysis(data, nocargoLftHistMean, kWeibull, cutoffIdx, 'w/o cargo');
