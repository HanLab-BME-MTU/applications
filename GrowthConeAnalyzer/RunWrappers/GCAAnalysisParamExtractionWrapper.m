function [ output_args ] = GCAAnalysisParamExtractionWrapper(projList)
%

for iProj = 1:numel(projList) 
    load([projList{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']); 
    %GCAGlobalThreshold(MD); % default otsu
    %GCAVisualsMakeOverlaysFilopodiaOldInputMovie(MD); 
    GCAAnalysisExtract_SpatialAC(MD); 
    GCAAddFilopodiaActinContentMetricMovie(MD); 
    GCAAnalysisExtractFilopodiaParamsMovie(MD);
    %GCAAnalysisExtract_MeanGrowthConeFluorescence(MD,true);  
end 

end

