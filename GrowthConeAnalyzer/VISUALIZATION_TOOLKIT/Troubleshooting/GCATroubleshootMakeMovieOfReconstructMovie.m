function [ output_args ] = GCATroubleshootMakeMovieOfReconstructMovie(movieData,frames)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

imDir  = movieData.getChannelPaths{1}; 


for iFrame = 1:length(frames) 
    
saveDir = [movieData.outputDirectory_ filesep 'Reconstruct_Frame' num2str(frames(iFrame),'%03d')]; 

if ~isdir(saveDir)
    mkdir(saveDir)
end 

load([movieData.outputDirectory_ filesep ...
    'SegmentationPackage' filesep 'StepsToReconstruct' ...
    filesep 'IV_veilStem_length' filesep 'Channel_1' filesep 'veilStem.mat']); 
%  load([movieData.outputDirectory_ filesep  'SegmentationPackage' filesep ... 
%     'StepsToReconstruct' filesep 'VI_filopodiaBranch_reconstruction' filesep 'Channel_1' filesep 'filoBranch.mat']); 

load([movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ...
    'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_1' filesep 'filoBranch.mat']);

pixSize_um = movieData.pixelSize_/1000; 

%load([movieData.outputDirectory_ filesep 'filopodia_reconstruct' filesep 'Filopodia_Reconstruct_Channel_1' filesep 'analInfoTestSave.mat']); 
 
%load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' filesep 'analInfoTestSave.mat'])
GCATroubleShootMakeMovieOfReconstructNew(filoBranch,veilStem,frames(iFrame),pixSize_um,saveDir,imDir); 




end 



end

