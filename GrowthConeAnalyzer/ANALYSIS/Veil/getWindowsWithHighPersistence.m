function [output] = getWindowsWithHighPersistence(analysisResults,cutoffPersistence)
%
% INPUT: 
% analysisResults: marcos analysisResultsStructure
% cutoffPersistence: edge movements with persistence times below this persistence will be filtered 
%
% OUTPUT: analysisResults.windows filtered for high persistence has the
% block data (ie the frame of the time series of the protrusion/retraction)
% and the window ID which can be translated back to coordinate. 
%

% remove the fields not associated with the window data 


for iWind = 1:length(analysisResults.protrusionAnalysis.windows)
%idxPer = cellfun(@(x) length(x)>cutoffPersistence, analysisResults.protrusionAnalysis.windows.persTime); 



%analysisResults.protrusionAnalysis.windows(iWind) = analysisResults.protrusionMeasurements.windows(idxPers); 
 idxP = analysisResults.protrusionAnalysis.windows(iWind).persTime>cutoffPersistence;
 idxR = analysisResults.retractionAnalysis.windows(iWind).persTime>cutoffPersistence;
 fields = fieldnames(analysisResults.protrusionAnalysis.windows(iWind));
 fields(strcmpi(fields,'limit')) = []; % get rid of the limit field not the
 % length of the windows. 
 for ifields = 1:length(fields) 
     x = analysisResults.protrusionAnalysis.windows(iWind).(fields{ifields})(idxP); 
     output.protrusionAnalysis.windows(iWind).(fields{ifields}) = x; 
     y = analysisResults.retractionAnalysis.windows(iWind).(fields{ifields})(idxR); 
     output.retractionAnalysis.windows(iWind).(fields{ifields}) = y; 
 end 
 




end

