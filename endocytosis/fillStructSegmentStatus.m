function [data] = fillStructSegmentStatus(data, fieldname)
% fill experiment structure with vectors designating the segmentation
% status of each trajectory (inside vs. outside pattern) 
% based on the image segmentation in the specified folders
% 
% SYNOPSIS [data] = fillStructSegmentStatus(data);
%
% INPUT     data: experiment structure, which has to contain the fields
%                .source
%                .movieLength
%                .segmentDataFileName 
%                .segmentDataFilePath 
%           fieldname (OPTIONAL)
%                   the segmentation image - or the path to its location - 
%                   can also be directly contained in another field 
%
% OUTPUT    creates new fields in the data structure called
%                .segmentStatus   
%           which is a vector containing the segmentation status for every
%           objects in the trackInfo; 
%           status = 1 means inside the segmented pattern
%           status = 0 means outside the segmented pattern
%           
%
% Dinah Loerke, modified July 30, 2008
% Francois Aguet, modified Feb 2010

for i = 1:length(data)
    
    fprintf('movie #%02d',i);
    
    if (nargin < 2) % use default field names
        SegmentMask = imread([data(i).segmentDataFilePath data(i).segmentDataFileName]);    
    else
        segmentation = data(i).(fieldname);
        if ischar(segmentation)
            SegmentMask = imread(segmentation);
        else
            SegmentMask = segmentation;
        end
    end
        
    if ~islogical(SegmentMask)
        SegmentMask = logical(SegmentMask/max(SegmentMask(:)));
    end    
        
    load([data(i).source 'LifetimeInfo' filesep 'lftInfo.mat']);

    % calculate segmentation status (1=Inside, 0=oustide segmented region)
    [segmentStatusVector, segmentEUdistVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
    
    % lifetime histograms
    data(i).segmentStatus = segmentStatusVector';
    data(i).segmentEUdist = segmentEUdistVector';
end