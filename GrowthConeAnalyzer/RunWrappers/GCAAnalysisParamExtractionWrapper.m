function [ output_args ] = GCAAnalysisParamExtractionWrapper(projList)
%

for iProj = 1:size(projList,1) 
    load([projList{iProj} filesep 'GrowthConeAnalyzer'   filesep 'movieData.mat']); 
    
    %GCAAnalysisExtract_SpatialAC(MD); 
    
    GCAAnalysisExtractFilopodiaMeasurementsMovieForSensitivityAnal(MD); 
    
end 

end

