function [exp] = determineMovieLength(exp)
% DETERMINEMOVIELENGTH determines the length of the movie from the
% available detection result file, for all the fields in an 
% experiment structure
%
% SYNOPSIS [exp] = determineMovieLength(exp)
%
% INPUT     exp    :    experiment structure, which has to contain a field
%                       .source
%                       which is the path to the data location; at this
%                       location, the function requires a folder named
%                       DetectionStructures, which in turn has to contain a
%                       variable called detection.mat
% OUTPUT    exp           
% REMARKS   the folder DetectionStructures and the variable detection.mat
% have been created by the function loadAndSaveDetection
%
% Dinah Loerke, last modified Jan 2008



for i=1:length(exp)
    % checl if movieLength value exists 
    loadML = 0;
    if isfield(exp,'movieLength')
        ml = exp(i).movieLength;
        if isempty(ml) | (ml==0)
            loadML = 1;
        end
    else
        loadML = 1;
    end
    
    if loadML==1
        
        od = cd;
        path = exp(i).source;
        detectionpath = [path,filesep,'DetectionStructures'];
        cd(detectionpath);
        
        readDetectionFilename = 'detection.mat';
        
        % load detection fiel, determine its length and save
        loadfile = load(readDetectionFilename);
        if isfield(loadfile,'detection')
            detection = loadfile.detection;
        else
            detection = [];
        end
            
        exp(i).movieLength = length(detection);  
        
        cd(od);
            
    end % of if
end % of for loop

end % of function