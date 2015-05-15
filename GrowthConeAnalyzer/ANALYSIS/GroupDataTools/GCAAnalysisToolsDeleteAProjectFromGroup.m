function [ toPlot ] = GCAAnalysisToolsDeleteAProjectFromGroup( toPlot,saveDir,filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if nargin< 2 || isempty(saveDir) 
    saveDir  = pwd ; 
end 

if nargin < 3 || isempty(filename)
   filename = 'toPlotMod.mat'; 
end 
    
    groupToMod = listSelectGUI(toPlot.info.names); 
    
    projList = toPlot.info.projList{groupToMod}; 
    projectsToDelete = listSelectGUI(projList); 
    projList(projectsToDelete) = []; 
    toPlot.info.projList{groupToMod} =  projList;  
    % fix the grouping 
    
    newGrouping = arrayfun(@(iGroup) repmat(iGroup,size(toPlot.info.projList{iGroup},1),1),1:numel(toPlot.info.names),... 
        'uniformoutput',0); 
    
    toPlot.info.grouping = vertcat(newGrouping{:}); 

    save([saveDir filesep filename],'toPlot'); 


