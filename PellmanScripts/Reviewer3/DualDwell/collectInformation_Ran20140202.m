function [ output_args ] = collectInformation( input_args )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dwellIntAll = cellfun(@(x) load([upDirectory(x,2) filesep 'DualAnalysisFinal_20150201' filesep 'dwellInt.mat']),groupList(:,2),'uniformoutput',0);

 numSig = arrayfun(@(x) vertcat(dwellIntAll{x}.dwellInt(:).enoughSig),1:numel(dwellIntAll),'uniformoutput',0)';
 idxST =  arrayfun(@(x) vertcat(dwellIntAll{x}.dwellInt(:).endType),1:numel(dwellIntAll),'uniformoutput',0)'; 
 cons = arrayfun(@(x) vertcat(dwellIntAll{x}.dwellInt(:).consistencyTest),1:numel(dwellIntAll),'uniformoutput',0)';  

 
 
 % collect for all groups
 numSigAll = vertcat(numSig{:}); 
 idxSTAll = vertcat(idxST{:}); 
 consAll = vertcat(cons{:});
 
 % just get the subtracks ending in a terminal or a shrinkage event. 
 numSigFinal = numSigAll((idxSTAll==1| idxSTAll ==3)) ; 
 
 % count the number of significant above the noise estimation from this population. 
 sig = sum(numSigFinal);
 % truncate the consistency test to include only terminal/rescue substracks
 % and signal that is significantly above the noise. 

 
 percentSig = sig/length(numSigFinal); 
 
 
 
 ranges1 = cell(numel(dwellIntAll),1); 
 ranges2 = cell(numel(dwellIntAll),1); 
 slope1 = cell(numel(dwellIntAll),1); 
 for i = 1:numel(dwellIntAll) 
    dwellInt = dwellIntAll{i}.dwellInt; 
    ranges1{i} = arrayfun(@(x) dwellInt(x).channel1.range1,1:length(dwellInt))';
    ranges2{i} = arrayfun(@(x) dwellInt(x).channel2.range2,1:length(dwellInt))';  
    slope1{i} = arrayfun(@(x) dwellInt(x).channel1.meanSlope,1:length(dwellInt))'; 
 end

 
ranges1Final= vertcat(ranges1{:}); 
ranges2Final = vertcat(ranges2{:}); 
slopesFinal = vertcat(slope1{:});
consFinal = consAll((idxSTAll==1| idxSTAll==3)&numSigAll==1& ~isnan(slopesFinal));

ranges1Final = ranges1Final(idxST==1 | idxST ==3); 
ranges2Final = ranges2Final(idxST==1 | idxST==3); 
boxplot([ranges1Final ranges2Final],'notch','on','colorGroup',[1,2],'colors',['g','r'],'outlierSize',1); 
ylabel({'Delta Fluorescence Intensity (AU)' ; 'In Dwell Window'}); 

end

