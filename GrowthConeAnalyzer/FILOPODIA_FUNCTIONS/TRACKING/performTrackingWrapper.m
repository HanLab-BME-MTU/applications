function [ output_args ] = performTrackingWrapper( projList  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
costMatParams.searchRadius = 5; 
costMatParams.predict = 0 ;


for iProj = 1:numel(projList)
  load([projList{iProj} filesep 'analInfoTestSave.mat']); 
   
  if isfield(analInfo,'imgPointer') ; 
      
  imDir = getFilenameBody(analInfo(1).imgPointer); 
 
  mainDir = getFilenameBody(imDir); 
  saveDir = [mainDir filesep 'tracking_noPredict'];
  imDir = [mainDir filesep 'imagesPretty']; 
  if ~isdir(saveDir)
      mkdir(saveDir)
  end 
  
  
      performTracking(analInfo,costMatParams,imDir, saveDir)

  clear imDir saveDir
  end 
end

