function [mask] = maskDetectedStructure_StandAlone(movieFeatures, nFrames, xMax, yMax, radius)
%maskDetectedStructure creates a mask that covers sub cellular structures detected by SubResolutionProcess
%
%SYNOPSIS function [] = maskDetectedStructure(MD, varargin)
%
%INPUT
%   movieFeatures : Contains information about subresolution structures
%                   deteced by SubResolutionProcess.m
%       .xCoord : x coordinates of detected strucutres in image coordinate
%                 system.
%       .ycoord : y coordinates of detected structures
%   nFrames       : number of frames in the movie
%   xMax          : size of the image in image coordinates
%   yMax          : size of the image
%   radius        : radius of imdilate stucture element
%
%OUPUT
%   mask : logical array with value true Mask array is in plot coordinate system [y, x]
%
%Tae H Kim, June 2015

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
%standard process input
ip.addRequired('movieFeatures', @isstruct);
ip.addRequired('nFrame', @isnumeric);
ip.addRequired('xMax', @isnumeric);
ip.addRequired('yMax', @isnumeric);
ip.parse(movieFeatures, nFrames, xMax, yMax);
%% mask creation
%memory allocation
mask = false(yMax, xMax, nFrames);
%set center of identified structure / feature to true
for iFrames = 1:nFrames
    [nStruct, ~] = size(movieFeatures(iFrames).xCoord);
    for iStruct = 1:nStruct
        x=round(movieFeatures(iFrames).xCoord(iStruct,1));
        y=round(movieFeatures(iFrames).yCoord(iStruct,1));
        mask(y,x,iFrames) = true;
    end
end
%imdilate expand around the centers of the identified strucutres to include
%the entirity of the structures
se = strel('disk', radius, 0);
for iFrames = 1:nFrames
    mask(:,:,iFrames) = imdilate(mask(:,:,iFrames), se);
end
end

