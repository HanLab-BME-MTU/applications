function [ output_args ] = GCACOutgrowthMovieWrapper(projList)
%projList = is a cell array of directories to run 
% 
% 
 for iProj = 1:numel(projList); 
     load([projList{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']); 
% 
     imDir = MD.getChannelPaths{1}; 
     anDir = MD.outputDirectory_; 
     
     load([projList{iProj} filesep 'ANALYSIS' filesep 'neurite_body_masks' filesep ... 
         'Neurite_Body_Masks_Channel_1' filesep 'analInfoTestSave.mat']); 
     
     
GCAfindNeuriteLength(analInfo,anDir,imDir,'k'); 

close all

end

