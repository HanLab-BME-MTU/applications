function [exp] = runDetection(exp)
% runDetection automatically runs Henry's detection software - with iclean
% 1 on all the movies specified in the experiment structure
% SYNOPSIS [exp] = loadAndSaveDetection(exp);
%
% INPUT     exp    :    experiment structure, which has to contain a field
%                       .source
%                       which is the path to the data location
% OUTPUT               
% REMARKS   The function only performs the detection for a given movie if 
%           there are no detection data (folder called maxdata283) already 
%           present in the specified directory; if you want a partial or
%           faulty detection to be replaced, you need to DELETE the folders
%           first!!!!
%
% Dinah Loerke, last modified Feb 2008
% Francois Aguet, last modified Dec 2009

% loop over all entries in the structure to enter the image data necessary for the detection input
nExp = length(exp);
imData(1:nExp) = struct('origImageName', [], 'origImagePath', []);

runStatus = zeros(1,nExp);

for i = 1:nExp
    [origImageName, origImagePath] = uigetfile('.tif', ['Select first original image for movie #' num2str(i)], exp(i).source);
    
    if (origImageName ~= 0)
        imData(i).origImageName = origImageName;
        imData(i).origImagePath = origImagePath;
        exp(i).firstImageName = origImageName;

        if ~(exist([exp(i).source filesep 'maxdata283'], 'dir') == 7)
            runStatus(i) = 1;
        else
            overwrite = input('Detection results already exist. Overwrite? (y/n): ', 's');
            if strcmp(overwrite, 'y')
                runStatus(i) = 1;
            end;
        end;
    end;
end;

for i = 1:nExp
    if runStatus(i)
        main283AUTO(1, imData(i).origImageName, imData(i).origImagePath);
    end;
end;