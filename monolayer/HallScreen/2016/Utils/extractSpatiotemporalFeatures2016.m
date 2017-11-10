function [allFeatures] = extractSpatiotemporalFeatures2016(mainDirname,metaData)
allDataDirname = [mainDirname '../../allData/'];

kymographDir = [allDataDirname 'kymographs/'];
healingRateDir = [allDataDirname 'healingRate/'];
measures = {'speed','directionality','coordination'};

allFeatures.speedFeats = extractFeatures(kymographDir,metaData,measures{1});

allFeatures.directionalityFeats = extractFeatures(kymographDir,metaData,measures{2});
% patch for fix directionality > 8
allFeatures.directionalityFeats.features(allFeatures.directionalityFeats.features > 10) = 10;

allFeatures.coordinationFeats = extractFeatures(kymographDir,metaData,measures{3});

allFeatures.healingRate = getHealingRate(healingRateDir,metaData);
end



%%
function [feats] = extractFeatures(kymographDir,metaData,measureStr)
feats.kymographs = cell(1,metaData.N);
feats.features = cell(1,metaData.N);
for i = 1 : metaData.N
    load([kymographDir measureStr '/' metaData.fnames{i} '_' measureStr 'Kymograph.mat']);
    if exist('speedKymograph','var')
        kymograph = speedKymograph;
    end
    if exist('directionalityKymograph','var')
        kymograph = directionalityKymograph;
    end
    if exist('strainRateKymograph','var')
        kymograph = strainRateKymograph;
    end
    if exist('coordinationKymograph','var')
        kymograph = coordinationKymograph;
    end
    
    assert(exist('kymograph','var') > 0);
    
    % TODO: deal with different kymographs
    feats.kymographs{i} = kymograph;    
end

[minDistFromWound,allDistancesFromWound] = getDistanceFromWound(feats.kymographs);

% Distance for calculations
feats.metaDists.minDistFromWound = minDistFromWound;
feats.metaDists.allDistancesFromWound = allDistancesFromWound;

[feats.features, feats.kymographs]= getFeatures(feats.kymographs,metaData,minDistFromWound);
end

%% Find maximal distance for which all values are calculated
function [minDistFromWound,allDistancesFromWound] = getDistanceFromWound(kymographs)
minDistFromWound = inf;
n = length(kymographs);
allDistancesFromWound = zeros(1,n);
for i = 1 : n
    curKymograph = kymographs{i};
    ind = find(isnan(curKymograph(:,1)),1,'first')-1;
    if isempty(ind)
        ind = size(curKymograph,1);
    end
    allDistancesFromWound(i) = ind;
    if ~isempty(ind)
        minDistFromWound = min(minDistFromWound,ind);
    end
end
end

%%
function [features kymographs] = getFeatures(kymographs,metaData,nDist)
n = length(kymographs);
nTime = floor(metaData.timeToAnalyze/metaData.timePerFrame);
nFeats = (metaData.timePartition-metaData.initialTimePartition) * metaData.spatialPartition; %metaData.timePartition * metaData.spatialPartition;
features = zeros(nFeats,n);

for i = 1 : n
    curKymograph = kymographs{i};    
    ys = 1 : floor(nDist/(metaData.spatialPartition)) : (nDist+1);
    xs = 1 : floor(nTime/(metaData.timePartition)) : (nTime+1);
    curFeatI = 0;
    for y = 1 : metaData.spatialPartition
        for x = (1+metaData.initialTimePartition) : metaData.timePartition%x = 1 : metaData.timePartition
            curFeatI = curFeatI + 1;
            values = curKymograph(ys(y):(ys(y+1)-1),xs(x):(xs(x+1)-1));
            values = values(~isinf(values));
            values = values(~isnan(values));
            features(curFeatI,i) = mean(values(:));
        end
    end
    kymographs{i} = curKymograph(1:nDist,1:nTime);
end

end

%% 
function [allHealingRates] = getHealingRate(healingRateDir,metaData)
allHealingRates = zeros(1,metaData.N);
nTime = floor(metaData.timeToAnalyze/metaData.timePerFrame);

for i = 1 : metaData.N
    clear averageHealingRate;
    load([healingRateDir '/' metaData.fnames{i} '_healingRate.mat'])
    allHealingRates(i) = averageHealingRate(nTime);    
end
end