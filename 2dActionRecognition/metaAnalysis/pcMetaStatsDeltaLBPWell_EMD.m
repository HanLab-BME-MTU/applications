%% At the well level: visualizing and quantifying differences between categories
% Assaf Zaritsky, June. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaStatsDeltaLBPWell_EMD()

addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
addpath(genpath('/home2/azaritsky/code/DeBiasUtils'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/dLBP/'];

outDname = [lbpDirname filesep 'FOV/'];

nScales = 4;


for iScale = 1 : nScales 
    outDnameScale = [outDname filesep num2str(iScale) filesep 'summary_emd' filesep];
    
    if ~exist(outDnameScale,'dir')
        unix(sprintf('mkdir %s',outDnameScale));
    end
    
    % Patch: load for dateAll,cellTypeAll
    patchFname = [analysisDirname 'metaAnalysis/LBPWell/FOV' filesep num2str(iScale) filesep 'allFeats.mat'];    
    load(patchFname);% 'mapsAll','featsAll','strAll','dateAll','cellTypeAll','cellTypeIndAll','sourceAll','metEffAll','n'    
    
    allFeatsFname = [outDnameScale filesep '../allFeats.mat'];
    load(allFeatsFname); % 'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allDeltaLBP','lowTH','highTH','n'
    
    
    featsAll = getDLBPFeats(distributionsAll);
    
    D = getSimilarityMatrixEMD(featsAll);
            
    % Check out if there are bugs in indices of similarity!
    %     pvalDayCellTypeSimilarity = pcMetaVarianceDayCellType(featsAll,dateAll,cellTypeAll,outDnameScale);
    pvalDayCellTypeSimilarity = pcMetaVarianceDayCellType(D,dateAll,cellTypeAll,outDnameScale);
    fprintf(sprintf('p(%d) = %.4f\n',iScale,pvalDayCellTypeSimilarity));
    
    %     [pvalCellTypeSource,pvalCellTypeSourceNoMelanocyte] = pcMetaVarianceCellTypeSource(featsAll,cellTypeAll,sourceAll,outDnameScale);
    [pvalCellTypeSource,pvalCellTypeSourceNoMelanocyte] = pcMetaVarianceCellTypeSource(D,cellTypeAll,sourceAll,outDnameScale);
    fprintf(sprintf('p(%d) = %.4f,%.4f\n',iScale,pvalCellTypeSource,pvalCellTypeSourceNoMelanocyte));
    
    
    pcMetaSimilarity(featsAll,D,sourceAll,outDnameScale,'Source');
    metEffAllCellArray = getMetEffCellArray(metEffAll);
    pcMetaSimilarity(featsAll,D,metEffAllCellArray,outDnameScale,'MetEff');
    pcMetaSimilarity(featsAll,D,cellTypeAll,outDnameScale,'CellType');
end
end


function [metEffAllCellArray] = getMetEffCellArray(metEffAll)
n = length(metEffAll);
metEffAllCellArray = cell(1,n);
for i = 1 : n
    if isnan(metEffAll{i})
        metEffAllCellArray{i} = '';
    else if metEffAll{i} == 0
            metEffAllCellArray{i} = 'Low';
        else if metEffAll{i} == 1
                metEffAllCellArray{i} = 'High';
            end
        end
    end    
end
end


function featsAll = getDLBPFeats(distributionsAll)
n = length(distributionsAll);
k = length(distributionsAll{1});
featsAll = nan(k,n);
for i = 1 : n
    featsAll(:,i) = distributionsAll{i}';
end
end


function D = getSimilarityMatrixEMD(featsAll)
n = size(featsAll,2);

D = zeros(n,n);

for i = 1 : n
    for j = 1 : n
        D(i,j) = EMD(featsAll(:,i),featsAll(:,j));        
    end
end

end