function [ analInfo ] = performMedianFiltering(analInfo,neuriteLengthInfoMovie,k)
%GCAfixBodyEstimateBasedOnOutlierInfo : find detection likely outlier frames based on
%neurite outgrowth median filter fitting currently we implement this as an optional
% step at the end 
% 
% INPUT: 
% analInfoC: with the body information  
% neuriteLengthInfoMovie: scalar 1 if need to calculate 
% 
% reset the outlier flag 

%analInfo = rmfield(analInfo,'flagOutlier'); 
% 
[outgrowthFilt , outlierIdx] = findOutliersFromMedFilt(neuriteLengthInfoMovie,10,k);
if ~isempty(outlierIdx)
% for each outlier put flag 
for i = 1:length(outlierIdx)
    frame = outlierIdx(i); 
    analInfo(frame).flagOutlier = 1;   
end 
 





end


