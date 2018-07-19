%% 
% LBP + deltaLBP
function [] = accLchSpaceTime120_All(featsStrID,featsStrOut,metaDataFname)

if nargin < 3
    featsStrID = 'LBP_dLBP';
    featsStrOut = 'LBP_dLBP_120_All201705';
    metaDataFname = 'ExperimentsAll20170514.mat';
end

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname 'Cells/dLBP/'];% input directory
params.featsStrIn = 'dLBP';
params.featsStrOut = featsStrOut;
params.featsStrID = featsStrID;
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/' metaDataFname]; % meta data (from excel DB)
params.fGetFeats = @getLchSpaceTimeAll120;

pcAccumulateFeatsGenericThirdGen(params);
end

function [feats,cellID,TXY] = getLchSpaceTimeAll120(cellData)
[feats,cellID,TXY] = getLchSpaceTimeT_10min(cellData,120);
end