function [nucleiStruc, dataProperties, mask, imageData] = spotsDetection3D(varargin)
% spotsDetection3D summarizes the workflow of 3D spots detection for 3D
%
% OUTPUT:
%   nucleiStruc:
%   dataProperties:
%   mask:
%   imageData:

% 07/2016 Ning Zhang

close all

p = inputParser;
p.addParameter('detectionMethod', 'mnp', @isstr);
p.addParameter('mannualAdjMode', 0, @isnumeric);
p.parse(varargin{:});
detectionMethod = p.Results.detectionMethod;
mannualAdjMode = p.Results.mannualAdjMode;

[imageData, dataProperties] = read3dFISH();
mask = segmentNucleiLocalOtsu(imageData.dapi, dataProperties);
nucleiStruc = findNuclei(imageData, mask, dataProperties);
pause(5);
nucleiStruc = singleNucleusSpotDetection(nucleiStruc, dataProperties, ...
    'detectionMethod', detectionMethod, 'mannualAdjMode', mannualAdjMode);

nucleiStruc = spotsPair(nucleiStruc, dataProperties);
spotsPlot3(nucleiStruc, imageData, dataProperties);
spotsPlot3Paired(nucleiStruc, imageData, dataProperties);
end

