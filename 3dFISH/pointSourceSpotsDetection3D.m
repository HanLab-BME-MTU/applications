function [pointSourceNucleiStruc, dataProperties, mask, imageData] = pointSourceSpotsDetection3D(varargin)
% pointSourceSpotsDetection3D summarizes the workflow of 3D spots detection
% using Francois' algorithm
%
% OUTPUT:
%   nucleiStruc:
%   dataProperties:
%   mask:
%   imageData:

% 07/2016 Ning Zhang

close all

[imageData, dataProperties] = read3dFISH();
mask = segmentNucleiLocalOtsu(imageData.dapi, dataProperties);
pointSourceNucleiStruc = pointSourceFindNuclei(imageData, mask, dataProperties);
pause(5);
pointSourceNucleiStruc = pointSourceSpotsDetection(pointSourceNucleiStruc, dataProperties);

pointSourceNucleiStruc = spotsPair(pointSourceNucleiStruc, dataProperties);
spotsPlot3(pointSourceNucleiStruc, imageData, dataProperties);
spotsPlot3Paired(pointSourceNucleiStruc, imageData, dataProperties);

end