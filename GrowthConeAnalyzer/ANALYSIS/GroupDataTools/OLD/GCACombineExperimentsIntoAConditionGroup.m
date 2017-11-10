function [ toPlotCondGroup ] = GCACombineExperimentsIntoAConditionGroup( expDirs,selectionFile,colorForGroup)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% expDirs: 
% cell array of exerimental directories 
% selectionFile : name - assumes you would like to pool all the selection
% files with the from each experiment with the same name 
% 
allToPlots = cellfun(@(x) load([x filesep 'SELECTIONS' filesep  selectionFile '.mat']),expDirs,'uniformoutput',0);
projLists = arrayfun(@(x) allToPlots{x}.toPlot.info.projList{1},1:numel(expDirs),'uniformoutput',0); 

projListCond  = vertcat(projLists{:});

toPlotCondGroup.info.names{1} = allToPlots{1}.toPlot.info.names{1}; 
toPlotCondGroup.info.projList{1} = projListCond; 
toPlotCondGroup.info.colors{1} = colorForGroup; 
toPlotCondGroup.info.grouping = ones(size(projListCond,1),1); 

end

