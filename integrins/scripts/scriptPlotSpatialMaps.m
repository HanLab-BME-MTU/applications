% clear all
% close all

movieName = 'AlphaVPax 121018 Cs2C3 Mode 4';

% masks(:,:,1) = double(imread('refined_mask_121018_Cs2C3_CHO_mEos2Av_PaxStack_0001.tif'));
% masks(:,:,2) = double(imread('refined_mask_121018_Cs2C3_CHO_mEos2Av_PaxStack_0037.tif'));
% 
% cd ../../../../analysisAlphaV/furtherAnalysis/
% 
% load tracksDiffusionLength5InMask.mat
% load diffusionModeClassification.mat
%         
% mkdir spatialMap
% cd spatialMap

plotPropertySpatialMap2D(tracksFinal,diffAnalysisRes,diffModeAnalysisRes,...
    6,1,[5 99],[],movieName,masks);

