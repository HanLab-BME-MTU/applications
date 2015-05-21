function [ output_args ] = GCATroubleshootMakeMovieOfReconstructMovie(movieData,frames)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

imDir  = movieData.getChannelPaths{1}; 


for iFrame = 1:length(frames) 
    
saveDir = [movieData.outputDirectory_ filesep 'Reconstruct_Frame' num2str(frames(iFrame),'%03d')]; 

if ~isdir(saveDir)
    mkdir(saveDir)
end 

%load([movieData.outputDirectory_ filesep 'filopodia_reconstruct' filesep 'Filopodia_Reconstruct_Channel_1' filesep 'analInfoTestSave.mat']); 

load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' filesep 'analInfoTestSave.mat'])
GCATroubleShootMakeMovieOfReconstruct(analInfo,frames(iFrame),.216,saveDir,imDir); 




end 



end

