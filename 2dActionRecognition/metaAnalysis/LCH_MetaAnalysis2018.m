%% Meta analysis given features
function [] = LCH_MetaAnalysis2018(featsDname,featsStr)

close all; clc;
always = true;

if nargin < 2
    featsDname = 'C:\Assaf\Research\LCH\data\CellExplorerData\2018';
    featsStr = 'LEVER_LBP_SHAPE';
end


%% Flags for processing
% Pre-processing: normalization (& concatenation of features?)
flags.orginizeFeatreus = true; % orginize data structure and normalize

% Analysis
flags.clsFeatures = true;
flags.clsFeatsVis = true;


%% Get list of features to process
accFeatsFnameAll = [featsDname filesep featsStr '_all.mat'];
feats2process = getFeatsToProcess(accFeatsFnameAll);
nfeats = length(feats2process);

for ifeat = 1 : nfeats
    curFeatStr = feats2process{ifeat};
    doMetaAnalysis(featsDname,curFeatStr,flags);
end

end

%%
function [] = doMetaAnalysis(featsDname,curFeatStr,flags)

if flags.orginizeFeatreus
    LCH_orginizeFeatures2018(featsDname,curFeatStr);
end

if flags.clsFeatures
    LCH_clsFeatures(baseDname,dateStr,featsStr); % Transform to new features
end

if flags.clsFeatsVis
    LCH_clsFeaturesVisualization_Train(baseDname,dateStr,featsStr);
end

end

%% collect features to process
function feats2process = getFeatsToProcess(accFeatsFnameAll)
load(accFeatsFnameAll); % allCells

feats = allCells.accFeats{1};
fields = fieldnames(feats);

feats2process = {};
nfeats = 0;

for ifield = 1 : numel(fields)
    curField = fields{ifield};        
    
    if strcmp(curField,'ncells') || strcmp(curField(1:4),'corr')
        continue;
    else
        nfeats = nfeats + 1;
        feats2process{nfeats} = curField;
    end
end
end




