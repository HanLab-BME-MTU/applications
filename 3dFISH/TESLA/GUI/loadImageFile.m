function dataFilePath = loadImageFile(varargin)
%LOADIMAGEFILE Get the path of the selected file and remember it for the future
%file loading
%   Detailed explanation goes here

% load history
% Search historyFile in the current folder
historyFile = fullfile(fileparts(mfilename('fullpath')),'historyFile.mat');
if exist(historyFile, 'file')
    history = load(historyFile);
    pathName = history.pathName;
    [fullFileName, pathName] = uigetfile({'*.tif'; '*.jpg'}, 'Select the image for analysis', pathName);
else
    [fullFileName, pathName] = uigetfile({'*.tif'; '*.jpg'}, 'Select the image for analysis');
end

if fullFileName ~= 0
    dataFilePath = fullfile(pathName, fullFileName);
else
    dataFilePath = '';
end

% save history
if ~isempty(pathName) && all(pathName ~= 0)
    save(fullfile(fileparts(mfilename('fullpath')),'historyFile.mat'), 'pathName');
end
end

