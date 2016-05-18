function [nucleiStruc, dataProperties, mask, imageData] = spotsDetection3D(varargin)
% spotsDetection3D summarizes the workflow of 3D spots detection for 3D
%
% OUTPUT:
%   nucleiStruc:
%   dataProperties:
%   mask:
%   imageData:

% 05/2016 Ning Zhang

close all

p = inputParser;
p.addParameter('detectionMethod', 'mnp', @isstr);
p.parse(varargin{:});
detectionMethod = p.Results.detectionMethod;

[imageData, dataProperties] = read3dFISH();
mask = segmentNucleiLocalOtsu(imageData.dapi, dataProperties);
nucleiStruc = findNuclei(imageData, mask, dataProperties);
nucleiStruc = singleNucleusSpotDetection(nucleiStruc, dataProperties, ...
    imageData, 'detectionMethod', detectionMethod);

end

