function [] = pcCollectGlobalData(exps,allParams,allDirs)
nexp = length(exps);

allData.ncells = 0;
allData.texture = [];
allData.accumulatedSizes = [];
allData.accumulatedEccentricity = [];
allData.accumulatedLocalSpeed = [];
allData.accumulatedLocalMatchingScore = [];
allData.cellID = [];
allData.cellExp = [];
for i = 1 : nexp
     [allData] = collectAllData(allData,allParams{i},allDirs{i},i);
end

save([allDirs{1}.results 'allData.mat'],'allData');

end

%%
function [allData] = collectAllData(allData,params,dirs,exp)
% Texture
textureDataFname = [dirs.results dirs.expname '_lbpPerFrame.mat'];
load(textureDataFname); % 'accumulatedLBP','accumulatedCellsID'

% Local morphodynamics
localMorphDynamicsFname = [dirs.results dirs.expname '_localMorphDynam.mat'];
load(localMorphDynamicsFname); 
% 'accumulatedSizes','accumulatedEccentricity','accumulatedLocalSpeed','accumulatedLocalMatchingScore','accumulatedCellsID'


ncells = size(accumulatedLBP,2);

allData.texture = [allData.texture, accumulatedLBP];

allData.accumulatedSizes = [allData.accumulatedSizes accumulatedSizes];
allData.accumulatedEccentricity = [allData.accumulatedEccentricity accumulatedEccentricity];
allData.accumulatedLocalSpeed = [allData.accumulatedLocalSpeed accumulatedLocalSpeed];
allData.accumulatedLocalMatchingScore = [allData.accumulatedLocalMatchingScore accumulatedLocalMatchingScore];

allData.ncells = allData.ncells + ncells;
allData.cellID = [allData.cellID, accumulatedCellsID]; % (pairs of time, #cells)
allData.cellExp = [allData.cellExp repmat(exp,1,ncells)]; % index for experiment
end