function [ DBIdx,DunnIdx,SC ] = GCAAnalysisRunDiscrimMetrics( toPlotMedian )
% toPlotGroupMedianPerMovie 
% tranform toPlot from a structure array with param fields to

% featCell: a nGroupx1 cell of rxc double arrays
% where r is the total number of features and c is the the number of movies measured.  
% 
 

featCell = cell(numel(toPlotMedian.info.names),1);
fields = fieldnames(toPlotMedian); 
fields  = fields(cellfun(@(x) ~strcmpi(x,'info'),fields)); 

for iGroup = 1:numel(toPlotMedian.info.names)
    
   features =   arrayfun(@(x) toPlotMedian.(fields{x}){iGroup}, 1:numel(fields),'uniformoutput',0); 
   
   
   featCell{iGroup,1} = vertcat(features{:}); 
    
end 



[DBIdx,DunnIdx,SC] = cellfun(@(x) whDiscriminationMeasures(x,featCell{1}),featCell); 


end

