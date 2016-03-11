function [ output_args ] = compareOutgrowthTwoGroups(toPlotGroups)
%

% setFigure 
    subplot(1,2,1); 

for iGroup = 1:numel(toPlot.info.names) 
    
    projListC = toPlot.info.projList{iGroup}; 
    color = toPlot.info.colors{iGroup}; 
   
    deltasC= GCAplotMultNetOutgrowth(projListC,color,[]); 
    
    
    deltas{iGroup} = deltasC; 
    hold on 

end

 subplot(1,2,2); 
 horzcat(deltas{:}); 
 notboxplot(

 
 


