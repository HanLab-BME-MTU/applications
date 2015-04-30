function [ deltas,neuriteLengths] = GCAplotMultNetOutgrowth( projListC,colorCurrent,emphasize )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


neuriteLengthStruct = cellfun(@(x) load([x filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION\GlobalFunctional\neurite_outgrowth_measurements' ...
    filesep 'neuriteLengthOutput.mat']),projListC);
neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0); 
deltas = cellfun(@(x) (x(end)-x(1)),neuriteLengths);  
 

hold on 
cellfun(@(x)  plotNeuriteOutgrowthInTime(x,colorCurrent,1,5,0,[],[],1),neuriteLengths);

if ~isempty(emphasize)    
     plotNeuriteOutgrowthInTime(neuriteLengths{emphasize},colorCurrent,1,5,0,[],[],1,5); 
end 

end

