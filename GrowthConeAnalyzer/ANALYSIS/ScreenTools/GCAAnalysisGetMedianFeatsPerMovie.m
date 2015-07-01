function [ toPlotMedian] = GCAAnalysisGetMedianFeatsPerMovie(toPlotAllFrames)
% 
% INPUT: toPlotAllFrames: 
% the toPlot structure with feature fields structured as rxc double arrays such 
% that r (row) is the number of frames per movie and c (column) is the
% movie number 
% 
% simply tranform so that each param is now an 1xc double array such that c is the movie number
% and the row is now the median value for the movie 

 

fields = fieldnames(toPlotAllFrames); 
% take out info 
fields = fields(cellfun(@(x) ~strcmpi(x,'info'),fields)); 

for iParam = 1:numel(fields)

 toPlotMedian.(fields{iParam})  = cellfun(@(x) nanmedian(x,1),toPlotAllFrames.(fields{iParam}),'uniformoutput',0); 

end
toPlotMedian.info = toPlotAllFrames.info; 

