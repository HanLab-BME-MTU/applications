function [ output_args ] = GCAmakeToPlotCondGroup( toSearch )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(toSearch) 
   toSearch =  uigetdir(pwd);
end 

directories = dir(toSearch);
idxDir = vertcat(directories(:).isdir); 
directories = directories(idxDir); 
forListSelect  = arrayfun(@(x) [toSearch filesep directories(x).name],1:length(directories),'uniformoutput',0);
idxDirectoryInclude  = listSelectGUI(forListSelect,[],'move'); 
 
 directories = directories(idxDirectoryInclude).name ; 




end

