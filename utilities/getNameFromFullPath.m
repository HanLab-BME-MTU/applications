function dataName = getNameFromFullPath(fullName)
%GETNAMEFROMFULLPATH is a utility to extract the name of a chromdyn experiment from the full path of the datafile
% 
% SYNOPSIS dataName = getNameFromFullPath(fullName)
%
% INPUT     fullName: full path name (string)
%
% OUTPUT    dataName: name of the experiment (string)
%
% c: jonas 06/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if ~isstr(fullName)
    error('GETNAMEFROMFULLPATH needs string as input')
end

%=============
% find name
%=============

% find last filesep
filesepList = findstr(fullName,filesep);

if isempty(filesepList)
    error('no valid filesep found in full path')
end

lastFilesep = filesepList(end);

% find data
dataStart = findstr(fullName,'data');

% read name (assume e.g,: \ndc10-1_37C_G1_001_corr-data-25-Mar-2004-14-32-00.mat)
dataName = fullName(lastFilesep+1 : dataStart-2);