%% 
% LBP + deltaLBP
function [] = accLchSpaceTime120_ThirdGen(featsStrID,featsStrOut,metaDataFname)

if nargin < 3
    featsStrID = 'LBP_dLBP';
    featsStrOut = 'LBP_dLBP_120_ThirdGen201701';
    metaDataFname = 'ThirdGen20170118_noControl.mat';
end

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname 'Cells/dLBP/'];% input directory
params.featsStrIn = 'dLBP';
params.featsStrOut = featsStrOut;
params.featsStrID = featsStrID;
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/' metaDataFname]; % meta data (from excel DB)
params.fGetFeats = @getLchSpaceTimeThirdGen120;

pcAccumulateFeatsGenericThirdGen(params);
end

function [feats,cellID,TXY] = getLchSpaceTimeThirdGen120(cellData)
[feats,cellID,TXY] = getLchSpaceTimeT(cellData,120);
end