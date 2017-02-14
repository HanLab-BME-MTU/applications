%% At the well level: combining features for all scales, visualizing and quantifying differences between categories
% Assaf Zaritsky, Oct. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaStatsDeltaLBPWell_allScales()

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
lbpDirname = [analysisDirname 'metaAnalysis/dLBP/'];

outDname = [lbpDirname filesep 'FOV/'];

nScales = 4;


featsAcrossScales = [];
for iScale = 1 : nScales 
    dnameScale = [outDname filesep num2str(iScale) filesep];      
    
    % Patch: load for dateAll,cellTypeAll
    patchFname = [analysisDirname 'metaAnalysis/LBPWell/FOV' filesep num2str(iScale) filesep 'allFeats.mat'];    
    load(patchFname);% 'mapsAll','featsAll','strAll','dateAll','cellTypeAll','cellTypeIndAll','sourceAll','metEffAll','n'    
    
    allFeatsFname = [outDname filesep num2str(iScale) filesep 'allFeats.mat'];
    load(allFeatsFname); % 'distributionsAll','distributionsLocations','meansAll','meansLocations','strLocations','allDeltaLBP','lowTH','highTH','n'
        
    featsAcrossScales = [featsAcrossScales; getDLBPFeats(distributionsAll)];
end
metEffAllCellArray = getMetEffCellArray(metEffAll);
% % pcMetaTSNE(featsAcrossScales,metEffAllCellArray,outDname,'MetEff');
% % pcMetaTSNE(featsAcrossScales,sourceAll,outDname,'Source');
% % pcMetaTSNE(featsAcrossScales,cellTypeAll,outDname,'CellType');
% % pcMetaTSNE(featsAcrossScales(:,~strcmp(sourceAll,'Melanocytes')),sourceAll(~strcmp(sourceAll,'Melanocytes')),outDname,'TumorCellLine');
pcMetaDeltaLbpNnCls(featsAcrossScales,metEffAllCellArray,cellTypeAll,dateAll,outDname,'MetEff');
indsNoMelano = ~strcmp(sourceAll,'Melanocytes');
pcKnnCls(featsAcrossScales(:,indsNoMelano),sourceAll(indsNoMelano),cellTypeAll(indsNoMelano),dateAll(indsNoMelano),outDname,'MetEff');
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


function feats = getDLBPFeats(distributionsAll)
n = length(distributionsAll);
k = length(distributionsAll{1});
feats = nan(k,n);
for i = 1 : n
    feats(:,i) = distributionsAll{i}';
end
end