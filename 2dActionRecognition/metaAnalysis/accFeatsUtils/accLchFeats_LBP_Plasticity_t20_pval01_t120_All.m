%% 
% Plasticity of LBP
function [] = accLchFeats_LBP_Plasticity_t20_pval01_t120_All(featsStrID,featsStrOut,metaDataFname)

if nargin < 3
    featsStrID = 'dLBP_Plasticity';
    featsStrOut = 'Plasticity_LBP_t20_pval01_t120';
    metaDataFname = 'ExperimentsAll20170514.mat';
    featsInDname = 'Cells/dLBP_Plasticity_t20_pval01/';
end

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname featsInDname];% input directory
params.featsStrIn = 'dLBP_Plasticity';
params.featsStrOut = featsStrOut;
params.featsStrID = featsStrID;
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/' metaDataFname]; % meta data (from excel DB)
params.fGetFeats = @getLchSpaceTimeAll120;

pcAccumulateFeatsGeneric201706(params);
end

function [feats,cellID,TXY] = getLchSpaceTimeAll120(featsInFname)
load(featsInFname); % 1 x nCells
cellData = dLbpWell;
[feats,cellID,TXY] = getFeatsT(cellData,dLbpPlasticityWell,120);% dLbpPlasticityWell
end

%% 
% LBP + deltaLBP
% take all cells with trajectories with frame time = t (frame)
function [feats, cellID, TXY] = getFeatsT(cellData,dLbpPlasticityWell,t)
n = length(cellData);
feats = [];
cellID = [];
TXY = {};
ncells = 0;
for i = 1 : n
    if sum(cellData{i}.ts==t)
        ncells = ncells + 1;
        feats = [feats,dLbpPlasticityWell{i}.plasticity];
        cellID = [cellID, i];        
        TXY{ncells}.ts = cellData{i}.ts;
        TXY{ncells}.xs = cellData{i}.xs;
        TXY{ncells}.ys = cellData{i}.ys;
    end
end
end