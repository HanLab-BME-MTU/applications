function [ output_args ] = helperCalcFreqOfProjectForProblem(analInfo)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for iFrame = 1:length(analInfo)
    filoInfo = analInfo(iFrame).filoInfo; 
    
  commonPix = arrayfun(@(x) intersect(x.Ext_pixIndicesFor,x.Ext_pixIndicesBack),filoInfo,'uniformoutput',0);
  incorrectPerFrame(iFrame) = sum(cellfun(@(x) ~isempty(x),commonPix));
  
end

