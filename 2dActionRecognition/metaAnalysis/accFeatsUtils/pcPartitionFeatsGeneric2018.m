%% Generic accumulation of *ONE* feature from cell type + condition (over all locations),
%% Assumes pcAccumulateFeatsGeneric2018 was 
%   Params:
%        featsOutDname  % input/output directory
%        featsStrOut    % feature string for the output files

function [] = pcPartitionFeatsGeneric2018(params)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

featsOutDname = params.featsOutDname;
featsStrOut = params.featsStrOut; % prefix for output of the features

if ~exist(featsOutDname,'dir')
    unix(sprintf('mkdir %s',featsOutDname));
end

% features
accFeatsFnameAll = [featsOutDname filesep featsStrOut '_all.mat'];

assert(logical(exist(accFeatsFnameAll,'file')));

load(accFeatsFnameAll); % allCells
allCellsOld = allCells;
clear allCells;

feats = allCellsOld.accFeats{1};

fields = fieldnames(feats);

for ifield = 1 : numel(fields)
    curField = fields{ifield};
    
    curFeatsFname = [featsOutDname filesep curField '_all.mat'];
    
    if strcmp(curField,'ncells') || strcmp(curField(1:4),'corr') 
        continue;
    else
        allCells = allCellsOld;
        allCells.featsStrID = curField;
        n = length(allCells.accFeats);
        for i = 1 : n
            allCells.accFeats{i} = allCells.accFeats{i}.(curField);
        end
    end
    save(curFeatsFname,'allCells');
end

end