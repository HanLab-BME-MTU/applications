function [data] = locateSegmentation(data, outvar);
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

od = cd;

ovariable = 0;
if nargin>1
    if outvar==1
        ovariable = 1;
    end
end

% loop over all entries in the structure to enter the image data necessary
% for the detection input
for i=1:length(data)
    
    path = data(i).source;
    % change to source directory
    cd(path);
    
    [oriImageName, oriImagePath] = uigetfile('.tif',['Select first segmentation image']); 
    
    if ovariable == 0;
        data(i).segmentDataFileName = oriImageName;
        data(i).segmentDataFilePath = oriImagePath;
    else
        data(i).segmentDataFileNameOUT = oriImageName;
        data(i).segmentDataFilePathOUT = oriImagePath;
    end
        
    cd(od);
end


end % of function