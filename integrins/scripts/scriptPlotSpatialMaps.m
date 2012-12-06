clear all
close all

movieName = 'AlphaVFixed 091012 CHO08lofix';

masks(:,:,1) = double(imread('refined_mask_MAX_imagesAlphaV_0001.tif'));
masks(:,:,2) = double(imread('refined_mask_MAX_imagesAlphaV_0002.tif'));

cd ../../../../analysisAlphaV/furtherAnalysis/

load tracksDiffusionLength5InMask.mat
load diffusionModeClassification.mat
        
mkdir spatialMap
cd spatialMap

plotPropertySpatialMap2D(tracksFinal,diffAnalysisRes,diffModeAnalysisRes,...
    6:7,1,[5 99],[],movieName,masks);

