function [ output_args ] = GCAAnalysis_PlotCorrelationsFilopodiaProperties( toPlot,saveDir )
%GCAAnalysis_PlotCorrelationFilopodiaProperties 
%     Quick check for correlation between filopodia actin intensity and
%     length/orientation in the control population. 
% 
%% Correlate Parameters 
projList = toPlot.info.projList{1}(:,1);
for iProj = 1:numel(projList)    
    % load the filopodia length to the veil 
    load([projList{iProj} filesep 'ANALYSIS' filesep ... 
        'PARAMETER_EXTRACTION_20150305' filesep 'Descriptor' filesep ...
        'Filopodia' filesep 'ConnectToVeil_LengthInt' filesep 'filoLength' filesep ...
        'param_filoLengthToVeil.mat']); 
     filoLength{iProj} = vertcat(paramC{:});
    
     % load the filopodia actin content from tip of the filopodia to the
     % veil
    load([projList{iProj} filesep 'ANALYSIS' filesep ... 
        'PARAMETER_EXTRACTION_20150305' filesep 'Descriptor' filesep ... 
        'Filopodia' filesep 'ConnectToVeil_LengthInt' filesep 'filoAvgIntensity'...
        filesep 'param_filoIntensityToVeil.mat']); 
    filoInt{iProj} = vertcat(paramC{:});      
    
    % load the filopodia actin content for the full bundle 
    
    
    
    % load the filopodia actin content for the embedded structure
    
    
    
    % load the filopodia length to the veil thickness (would need to
    % extract filopodia just from this region to show explicitly) but can
    % just get a feel for now. 
    
    
end 
    
    
save([saveDir filesep 'filoLength.mat']); 
save([saveDir filesep 'filoInt.mat']); 





end

