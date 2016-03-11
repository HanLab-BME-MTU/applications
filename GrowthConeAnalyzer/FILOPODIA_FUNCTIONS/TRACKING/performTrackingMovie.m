function [ output_args ] = performTrackingMovie(movieData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% params 
costMatParams.searchRadius = 5; 
costMatParams.predict = 0 ;
%% 20140915 for now just load analInfo automatically from reconstruct channel 1 
% and images from donor channel under the current output set-up 

% load from 

  load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep ... 
      'Filopodia_Fits_Channel_1' filesep 'analInfoWithWindInfo.mat']);
  
   imDir = movieData.channels_(1).channelPath_;
  
      
 
 
  
  saveDir = [movieData.outputDirectory_ filesep 'tracking_noPredict'];
  if ~isdir(saveDir)
      mkdir(saveDir) 
  end 
  
  
  
      performTracking(analInfo,costMatParams,imDir, saveDir)

  
  end 


