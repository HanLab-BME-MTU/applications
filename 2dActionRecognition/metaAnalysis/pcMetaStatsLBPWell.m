%% At the well level: visualizing and quantifying differences between categories
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaStatsLBPWell()

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/LBPWell/'];

outDnameFov = [lbpDirname filesep 'FOV/'];
outDnameBck = [lbpDirname filesep 'BCK/'];
outDnameFwd = [lbpDirname filesep 'FWD/'];

% fovStrs = {'fov','bck','fwd'};

nScales = 4; % 4 scales (1 to 1/8)
% scales = 1.0./2.^((1:nScales)-1);
% nScales = 1; % Debug 201604
% scales = 1; % Debug 201604

% outDirs = {outDnameFov,outDnameBck,outDnameFwd}; 
outDirs = {outDnameFov}; % Debug 201604

% testStrs = {'fov'}; % Debug 201604

nTest = length(outDirs);

for iTest = 1 : nTest 
    outDname = outDirs{iTest};    
    for iScale = 1 : nScales % resolution (1- maximal)
        outDnameScale = [outDname filesep num2str(iScale) filesep];
        
        allFeatsFname = [outDnameScale filesep 'allFeats.mat'];
        
        load(allFeatsFname); % 'mapsAll','featsAll','strAll','dateAll','cellTypeAll','cellTypeIndAll','sourceAll','metEffAll','n'
        
        %         Dvec = pdist(featsAll');
        %         D = squareform(Dvec);
        
        %% Check out if there are bugs in indices of similarity!
        %         pvalDayCellTypeSimilarity = pcMetaVarianceDayCellType(featsAll,dateAll,cellTypeAll,outDnameScale);
        %         fprintf(sprintf('p(%s,%d) = %.4f\n',fovStrs{iTest},iScale,pvalDayCellTypeSimilarity));
        
        %         [pvalCellTypeSource,pvalCellTypeSourceNoMelanocyte] = pcMetaVarianceCellTypeSource(featsAll,cellTypeAll,sourceAll,outDnameScale);
        %         fprintf(sprintf('p(%s,%d) = %.4f,%.4f\n',fovStrs{iTest},iScale,pvalCellTypeSource,pvalCellTypeSourceNoMelanocyte));
        
                
        pcMetaSimilarity(featsAll,sourceAll,outDnameScale,'Source');
        metEffAllCellArray = getMetEffCellArray(metEffAll);
        pcMetaSimilarity(featsAll,metEffAllCellArray,outDnameScale,'MetEff');
        pcMetaSimilarity(featsAll,cellTypeAll,outDnameScale,'CellType');
    end
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