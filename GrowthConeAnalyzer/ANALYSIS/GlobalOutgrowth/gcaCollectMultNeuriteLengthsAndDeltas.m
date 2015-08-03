function [ neuriteLengths,deltas] = gcaCollectMultNeuriteLengthsAndDeltas(projListC)

 idxNotRun = cellfun(@(x) exist( [x filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION\GlobalFunctional\neurite_outgrowth_measurements' ...
    filesep 'neuriteLengthOutput.mat'],'file') == 0, projListC); 
projListToRun = projListC(idxNotRun);
if ~isempty(projListToRun)
GCAfindNeuriteLengthWrapper(projListToRun); 
end 
neuriteLengthStruct = cellfun(@(x) load([x filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION' ...
    '\GlobalFunctional\neurite_outgrowth_measurements' ...
    filesep 'neuriteLengthOutput.mat']),projListC);
% find those not run 

neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0); 
deltas = cellfun(@(x) (x(end)-x(1)),neuriteLengths);  


end

