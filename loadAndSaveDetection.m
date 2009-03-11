function [exp] = loadAndSaveDetection(exp);
% LOADANDSAVEDETECTION loads the results from Henry's spot detection
% software, and converts them into the movieInfo format suitable for
% tracking, for all the fields in an experiment structure
%
% SYNOPSIS [exp] = loadAndSaveDetection(exp);
%
% INPUT     exp    :    experiment structure, which has to contain a field
%                       .source
%                       which is the path to the data location; at this
%                       location, the function reads the detection data
%                       from a folder called maxdata283, and writes the new
%                       movieInfo format file into a folder called
%                       DetectionStructures
% OUTPUT    exp           
% REMARKS 
%
% Dinah Loerke, last modified Jan 2008


% loop over all entries in the structure
for i=1:length(exp)
    
    % check if detection structure already exists
    path = exp(i).source;
    
    % change to source directory
    od = cd;
    cd(path);
    
    % value of loadvar determines whether the detection data will be loaded
    % and converted; if the detection.mat file already exists, it will stay
    % at value 0 (no loading or conversion), else it's set to 1 (loading
    % and conversion)
    loadvar = 0;
    % if detection structures folder doesn't already exist, create it
    if ~ ( exist('DetectionStructures')==7 )
        mkdir('DetectionStructures');
        loadvar = 1;
    else
        cd('DetectionStructures');
        if ~(exist('detection.mat')==2)
            loadvar = 1;
        end
    end
    
    cd(od);
    
    % if data should be loaded and converted, do it here
    if loadvar == 1
        maxdatapath = [path,'/maxdata283'];
        disp(['loading detection data for movie ',num2str(i)]);
        [detection] = convertDetectDataForTracking(maxdatapath);
        detectiondatapath = [path,'/DetectionStructures'];
        cd(detectiondatapath);
        save('detection.mat','detection');
        cd(od);
    end
end

end % of function