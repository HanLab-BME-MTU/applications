function [exp] = runDetection(exp, overwrite)
% runDetection automatically runs Henry's detection software - with iclean
% 1 on all the movies specified in the experiment structure
% SYNOPSIS [exp] = loadAndSaveDetection(exp);
%
% INPUT     exp       : experiment structure, containing a '.source' field (path to data location)
%           overwrite : optional, 1 to overwrite previous detection results (default 0).
%                       
% OUTPUT    none 

% Dinah Loerke, last modified Feb 2008
% Francois Aguet, last modified Mar 2010

if nargin<2
    overwrite = 0;
end

% loop over all entries in the structure to enter the image data necessary for the detection input
nExp = length(exp);
imData(1:nExp) = struct('origImageName', [], 'origImagePath', []);

runStatus = zeros(1,nExp);

for i = 1:nExp
    [origImageName, origImagePath] = uigetfile({'*.tif;*.tiff'}, ['Select first original image for movie #' num2str(i)], exp(i).source);
    
    if (origImageName ~= 0)
        imData(i).origImageName = origImageName;
        imData(i).origImagePath = origImagePath;
        exp(i).firstImageName = origImageName;

        if ~(exist([exp(i).source filesep 'maxdata283'], 'dir') == 7)
            runStatus(i) = 1;
        else
            if (overwrite)
                fprintf('Overwriting detecting results for movie %d.\n', i);
                runStatus(i) = 1;
            end;
        end;
    end;
end;

parfor i = 1:nExp
    if runStatus(i)
        main283AUTO(1, imData(i).origImageName, imData(i).origImagePath);
    end;
end;