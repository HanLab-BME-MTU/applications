function [data] = fillStructSegmentStatus(data, fieldname);
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

segmDefault = 0;
if nargin>1
    segmDefault = 1;
end

% number of entries in data structure
lens = length(data);

od = cd;

for i=1:lens
    
    fprintf('movie #%02d',i);
    
    % current path
    path = data(i).source;
    
    % number of frames for this exp
    lenf = data(i).movieLength;
    
    if segmDefault==0
        % load segmentation image from specified location
        SegmFileName = data(i).segmentDataFileName;
        SegmFilePath = data(i).segmentDataFilePath;

        cd(SegmFilePath);
        SegmentMask = imread(SegmFileName);    
        
    else
        segmentation = getfield(data(i), fieldname);
        if isstr(segmentation)
            SegmentMask = imread(segmentation);
        else
            SegmentMask = segmentation;
        end
    end
        
    if ~islogical(SegmentMask)
        SegmentMask = logical(SegmentMask/max(SegmentMask(:)));
    end    
        
    % load lifetime info data file
    lftpath = [path,'/LifetimeInfo'];
    cd(lftpath);
    lftname = 'lftInfo.mat';
    loadfile = load(lftname);
    cd(od);
    
    lftInfo = loadfile.lftInfo;

    % calculate segmentation status (1=Inside, 0=oustide segmented region)
    [segmentStatusVector,segmentEUdistVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
    
    
    % lifetime histograms
    data(i).segmentStatus = segmentStatusVector';
    data(i).segmentEUdist = segmentEUdistVector';
    
        
    fprintf('\b\b\b\b\b\b\b\b\b');       

end % of for

fprintf('\n'); 

end % of function
