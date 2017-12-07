%% 
% LBP + deltaLBP
function [] = accLchSpaceTime60_All_201711(featsStrID,featsStrOut,metaDataFname)

if nargin < 3
    featsStrID = 'LBP_dLBP';
    featsStrOut = 'LBP_dLBP_60_All201711';
    metaDataFname = 'ExperimentsAll20170810.mat';
end

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname 'Cells/dLBP/'];% input directory
params.featsStrIn = 'dLBP';
params.featsStrOut = featsStrOut;
params.featsStrID = featsStrID;
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/' metaDataFname]; % meta data (from excel DB)
params.fGetFeats = @getLchSpaceTimeAll60;

pcAccumulateFeatsGenericThirdGen(params);
end

function [feats,cellID,TXY] = getLchSpaceTimeAll60(cellData)
[feats,cellID,TXY] = getLchSpaceTimeT(cellData,60); % takes the whole trajectory
end