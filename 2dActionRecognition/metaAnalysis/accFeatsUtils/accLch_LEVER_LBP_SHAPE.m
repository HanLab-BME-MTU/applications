%% 
% Pooling LBP and shape feature across experiments using LEVER segmentation
% 
function [] = accLch_LEVER_LBP_SHAPE(featsStrID,featsStrOut,metaDataFname)

if nargin < 3
    featsStrID = 'LEVER_LBP_SHAPE';
    featsStrOut = 'LEVER_LBP_SHAPE';
    metaDataFname = 'ExperimentsAll20171016a.mat';
end

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname 'Cells/LEVER_LBP_SHAPE/'];% input directory
params.featsStrIn = 'LEVER_LBP_SHAPE';
params.featsStrOut = featsStrOut;
params.featsStrID = featsStrID;
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/' metaDataFname]; % meta data (from excel DB)
params.fGetFeats = @getLchWrapper_LEVER_LBP_SHAPE;

pcAccumulateFeatsGeneric2018(params); % Pool all features to a single data-structure
pcPartitionFeatsGeneric2018(params); % Partition to individual features
end

function [feats,cellID,TXY] = getLchWrapper_LEVER_LBP_SHAPE(dataTask)
[feats,cellID,TXY] = getLch_LEVER_LBP_SHAPE(dataTask);
end