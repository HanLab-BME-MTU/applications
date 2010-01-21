function [data] = locateSegmentation(data, outvar)
% locate segmentation data in the specified paths, and store the location 
% as fields in the data structure
%
% INPUT     data    :    experiment structure, which has to contain a field
%                       .source
%                        which is the path to the data location
%           outvar (optional) : if outvar is set to 1 (default is zero),
%                        then the located segmentation data is considered
%                        the segmentation for the OUTSIDE (e.g. cell
%                        outline), as opposed to the INSIDE (e.g. adhesive
%                        pattern outline). If this is the case, then the
%                        results are written to a different field name (see
%                        below)
% OUTPUT    two new fields are added to (or overwritten in) data, which are
%                   data.segmentDataFilePath
%                   data.segmentDataFileName
%           IF the value of input variable outvar==1, then the fields are
%           called
%                   data.segmentDataFilePathOUT
%                   data.segmentDataFileNameOUT
%           to avoid overwriting the pattern segmentation path
%
% NOTE: If the folder containing the segmentation data has a standardized
% name (e.g. if it's always called 'SegmentImages' or something along these
% lines), the function can be adapted to be more user-friendly by
% automatically setting the path to the standard-name subfolder, instead of
% just the source path - this would mean less user clicking
%
% Dinah Loerke, last modified April 20, 2008
% Dinah Loerke, updated July 23, 2008
% Francois Aguet, 01/21/2010

if (nargin == 1)
    outvar = 0;
end

for i = 1:length(data)
    
    [imageName, imagePath] = uigetfile('.tif', 'Select first segmentation image:', data(i).source);
    
    if outvar == 0;
        data(i).segmentDataFileName = imageName;
        data(i).segmentDataFilePath = imagePath;
    else
        data(i).segmentDataFileNameOUT = imageName;
        data(i).segmentDataFilePathOUT = imagePath;
    end
end