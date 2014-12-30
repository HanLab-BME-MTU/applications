function [ output_args ] = GCAVeilStemRefineScaleActiveContourMovie(movieData)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
imDir = movieData.getChannelPaths{1}; 

refineDir =  [ movieData.outputDirectory_ filesep 'neurite_veilStem_refinements']; 

if isdir(refineDir)
mkdir(refineDir)
end 

 load([movieData.outputDirectory_ filesep 'neurite_body_masks'...
    filesep 'Neurite_Body_Masks_Channel_1' filesep 'analInfoTestSave.mat']); 

framesList = 1:length(analInfo);
GCAveilStemRefineScaleActiveContour(imDir,refineDir,analInfo,framesList); 




