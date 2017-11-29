function [ output_args ] = GCAGroupAnalysisPlotHMMTransitions(toPlot,varargin)
% GCAGroupAnalysisPlotHMMTransitions

% 
groupIDs = [9,10]; 
nGroups = length(groupIDs); 
setAxis('on')
hold on

% filter by the max natural magnitude of the control 



for iGroup = 1:nGroups
    nProjs = size(toPlot.info.projList{groupIDs(iGroup)},1); 
    cmap = toPlot.info.colorShades{groupIDs(iGroup)}; 
    cmapFinal = [cmap(1,:) ; cmap(5,:); cmap(10,:)]; 
    
    
    
    for iProj = 1:nProjs 
       x =  toPlot.HMM.transTime{groupIDs(iGroup)}{iProj};
       
       if ~isempty(x)
           idx = find(x > 61,1); 
           y = toPlot.HMM.transMag{groupIDs(iGroup)}{iProj}; 
           firstY = y(idx); 
           firstX = x(idx); 
           scatter(firstY,firstX,50,cmapFinal(iProj,:),'filled');
       end 
           
 
    end 
end

