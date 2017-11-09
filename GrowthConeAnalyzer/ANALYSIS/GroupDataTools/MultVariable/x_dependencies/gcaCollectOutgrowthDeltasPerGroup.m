function [ deltasPerGroup ] = gcaCollectOutgrowthDeltasPerGroup(toPlotGroup)
% small function to grap the neurite outgrowth metrics 

for iGroup = 1:numel(toPlotGroup.info.names);
    projListC = toPlotGroup.info.projList{iGroup};
    projListC = projListC(:,1); 
    neuriteLengthStruct = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'MEASUREMENT_EXTRACTION' ...
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteLengthOutput.mat']),projListC);
    
    
    neuriteLengths  = arrayfun(@(x) x.neuriteLength, neuriteLengthStruct,'uniformoutput',0);
    deltasPerGroup{iGroup} = cellfun(@(x) (x(end)-x(1)),neuriteLengths);
end

end

