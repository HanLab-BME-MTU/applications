function [ toPlot ] = GCAAnalysisToolsRemoveAGroup( toPlot,saveDir,filename )
% GCAAnalysisToolsRemoveAGroup
if nargin<2 || isempty(saveDir) 
    saveDir = pwd ; 
end 

if nargin<3 || isempty(filename) 
    filename = 'toPlotGroupRemoved'; 
end 

gNames = toPlot.info.names;
groupToDelete = listSelectGUI(gNames);
nGroups = numel(toPlot.info.names); 


% fix info
fieldsInfo = fieldnames(toPlot.info);
% grouping field is set up differently
fieldsInfo = fieldsInfo(cellfun(@(x) ~strcmpi(x,'grouping'),fieldsInfo));

for iParam = 1:numel(fieldsInfo)
    toPlot.info.(fieldsInfo{iParam})(groupToDelete) = [];
end

grouping = arrayfun(@(x) repmat(x,size(toPlot.info.projList{x},1)),1:nGroups-length(groupToDelete),...
    'uniformoutput',0);
toPlot.info.grouping = grouping; 

fieldsRest = fieldnames(toPlot); 
fieldsRest = fieldsRest(cellfun(@(x) ~strcmpi(x,'info'),fieldsRest));

for iParam = 1:numel(fieldsRest) 
toPlot.(fieldsRest{iParam})(groupToDelete) = [];  
end 
save([saveDir filesep filename],'toPlot'); 
end

