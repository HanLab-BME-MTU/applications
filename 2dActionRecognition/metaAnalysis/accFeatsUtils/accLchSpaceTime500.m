%% 
% LBP + deltaLBP
function [] = accLchSpaceTime500()

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

params.featsInDname = [analysisDirname 'Cells/dLBP/'];% input directory
params.featsStrIn = 'dLBP';
params.featsStrOut = 'LBP_dLBP_500';
params.featsOutDname  = [analysisDirname 'metaAnalysis/' params.featsStrOut]; % output directory
params.metaDataFname  = [analysisDirname 'MetaData/Experiments20151023.mat']; % meta data (from excel DB)
params.fGetFeats = @getLchSpaceTime500;

pcAccumulateFeatsGeneric(params);
end