function [ output_args ] = GCAAnalysisToolsAddAProjectTotoPlot(toPlot )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% find the group to which to add the data 
% get names of all the groups 
names = toPlot.info.names;  
userPick = menu('Please Select a Group To Add Project',names); 
path = uigetdir(pwd,'Please Select Project to Add'); 

  
projList = toPlot.info.projList{userPick}; 



projList{end+1,1} = path; 
[~,~,num] = upDirectory(path,1); 
[~,date] = upDirectory(path,2); 
[~,name] = upDirectory(path,3); 
% figure out where it belongs 
%idx = isempty(strfind(name,names)); 

projList{end,2} = [name '_' date '_' num]; 

toPlot.info.projList{userPick} = projList; 
% update groupinfo 
grouping = toPlot.info.grouping; 



save('toPlotUpdated','toPlot'); 

end

