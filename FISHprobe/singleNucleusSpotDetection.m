function [nucleiStruc, dataProperties] = singleNucleusSpotDetection(nucleiStruc, dataProperties, imageData, varargin)
%singleNucleusSpotDetection detects spots from each multi-channel 3D stack
%of single nucleus
%   Detailed explanation goes here

% 04/2016 Ning

p = inputParser;
p.addRequired('nucleiStruc', @(x) true);
p.addRequired('dataProperties', @(x) true);
p.addRequired('imageData', @(x) true);
p.addOptional('detectionMethod', 'mnp', @isstr);
p.parse(nucleiStruc, dataProperties, imageData, varargin{:});

nucleiStruc = p.Results.nucleiStruc;
dataProperties = p.Results.dataProperties;
imageData = p.Results.imageData;
detectionMethod = p.Results.detectionMethod;

% Define dataProperties parameters
% Generate dataProperties.FILTERPRM
% calcFilterParms generates psf size, which is used to define patchsize and
% involved in gaussian filter
% Code borrowed and modified from 
% /home2/nzhang/matlab/applications/FISHprobe/Spots/detect/defaultDataProperties.m

for chaNum = 1:numel(dataProperties.channel)
    chaName = dataProperties.channel(chaNum).name;
    switch chaName
        case 'dapi'
            continue;
            
        case 'green' 
            greenWVL = dataProperties.channel(chaNum).emissionWavelength;
            [FT_XY, FT_Z] = calcFilterParms(greenWVL, dataProperties.NA, ...
                            dataProperties.refractiveIndex, 'gauss', ...
                            dataProperties.sigmaCorrection, ...
                            [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);
            % FT_XY vs psfsigma value from bioformat reader??
            patchXYZ = roundOddOrEven(4*[FT_XY FT_XY FT_Z], 'odd', 'inf');
            dataProperties.channel(chaNum).FILTERPRM = [FT_XY, FT_XY, FT_Z, patchXYZ];
            dataProperties.channel(chaNum).FT_SIGMA = [FT_XY, FT_XY, FT_Z];
            
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).green;
                nucStack = filtermovie(nucStack, dataProperties.channel(chaNum).FILTERPRM, 0);
                switch detectionMethod
                    case 'mnp'
                        spots = spotFindSingleNucMNP(nucStack, dataProperties.channel(chaNum));
                    case 'intensity'
                        spots = spotFindSingleNucIntensity(nucStack, dataProperties.channel(chaNum));
                end
                nucleiStruc(nucNum).greenSpot = spots.sp;
                
            end
            
        case 'red'
            redWVL = dataProperties.channel(chaNum).emissionWavelength;
            [FT_XY, FT_Z] = calcFilterParms(redWVL, dataProperties.NA, ...
                dataProperties.refractiveIndex, 'gauss', ...
                dataProperties.sigmaCorrection, ...
                [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);
            
            patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z], 'odd', 'inf');
            dataProperties.channel(chaNum).FILTERPRM = [FT_XY, FT_XY, FT_Z, patchXYZ];
            dataProperties.channel(chaNum).FT_SIGMA = [FT_XY, FT_XY, FT_Z];
            
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).red;
                nucStack = filtermovie(nucStack, dataProperties.channel(chaNum).FILTERPRM, 0);
                switch detectionMethod
                    case 'mnp'
                        spots = spotFindSingleNucMNP(nucStack, dataProperties.channel(chaNum));
                    case 'intensity'
                        spots = spotFindSingleNucIntensity(nucStack, dataProperties.channel(chaNum));
                end
                nucleiStruc(nucNum).redSpot = spots.sp;
                
            end
            
        otherwise
            error('Unknown channels detected')
    end
   
end

    spotsPlot3(nucleiStruc, imageData, dataProperties);
    
end