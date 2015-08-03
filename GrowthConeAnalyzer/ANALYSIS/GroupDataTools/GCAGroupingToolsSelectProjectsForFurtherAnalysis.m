function [ output_args ] = GCASelectProjectsForFurtherAnalysis(toSearch,selectionName)

if isempty(toSearch) 
   toSearch =  uigetdir(pwd);
end 

directories = dir(toSearch);
idxDir = vertcat(directories(:).isdir); 
directories = directories(idxDir); 
forListSelect  = arrayfun(@(x) [toSearch filesep directories(x).name],1:length(directories),'uniformoutput',0);
idxDirectoryInclude  = listSelectGUI(forListSelect,[],'move'); 
 
 directories = directories(idxDirectoryInclude).name ; 
 
 
 for iDir = 1:length(directories) 
   toPlot =   GCAAnalysisToolsMakeToPlotFile(directories(iDir),1); 
   outDir = [directories(iDir) filesep 'Grouping']; 
   if ~isdir(outDir) 
       
   mkdir(outDir)
   end 
  save([outDir filesep 'toPlot' selectionName '.mat'],'toPlot');      
 end 
end

