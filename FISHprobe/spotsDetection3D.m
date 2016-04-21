function [nucleiStruc, dataProperties] = spotsDetection3D()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[imageData, dataProperties] = read3dFISH();
mask = segmentNucleiLocalOtsu(imageData.dapi, dataProperties);
imseriesmaskshow(imageData.dapi, mask);
nucleiStruc = findNuclei(imageData, mask, dataProperties);
[nucleiStruc, dataProperties] = singleNucleusSpotDetection(nucleiStruc, dataProperties, imageData);

end

