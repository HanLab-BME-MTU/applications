function [data] = loadConditionDataMultiChannel(condDir, chNames, markers)
% loadConditionData loads the relevant information for all the data
% available for a specific data condition; this requires a specific
% dircetory structure and nomenclature (see below)
%
% SYNOPSIS [data] = loadConditionData()
%
% INPUT
%
% OUTPUT   data: structure with the fields
%                   .source = pathname of the data/movie
%                   .date = date when movie was taken
%                   .framerate = framerate of the movie (2s or 0.4s)
%
%
%
% Francois Aguet, October 2010

if nargin<1
    condDir = [uigetdir(pwd, 'Select the ''condition'' folder:') filesep];
end
if ~strcmp(condDir(end), filesep)
    condDir = [condDir filesep];
end

fprintf('Root directory: %s\n', condDir);


% list of experiments for this condition, each containing one or more 'cell' directories
expDir = dirList(condDir);

% if expDir are 'cell' directories
valid = cell2mat(regexpi(arrayfun(@(x) x.name, expDir, 'UniformOutput', false), 'cell', 'once'));
if ~isempty(valid)
    cellPath = arrayfun(@(x) [condDir x.name filesep], expDir, 'UniformOutput', false);
else
    cellPath = arrayfun(@(x) arrayfun(@(y) [condDir x.name filesep y.name filesep], dirList([condDir x.name]), 'UniformOutput', false), expDir, 'UniformOutput', false);
    cellPath = vertcat(cellPath{:});
end

% check whether directory names contain 'cell'
valid = cellfun(@(x) regexpi(getDirFromPath(x), 'cell'), cellPath);
cellPath = cellPath(valid==1);
nCells = length(cellPath);

% no 'cell' folders are found
if nCells == 0
    error('No data found in directory.');
end

data(1:nCells) = struct('source', [], 'channels', [], 'date', [], 'framerate', [],...
    'imagesize', [], 'movieLength', [], 'markers', []);


% Load/determine channel names
if nargin<2
    nCh = input('Enter the number of channels: ');
    chNames = cell(1,nCh);
    chPath = [uigetdir(cellPath{1}, 'Select first (master) channel:') filesep];
    chNames{1} = chPath(length(cellPath{1})+1:end-1);
    for c = 2:nCh
        chPath = [uigetdir(cellPath{1}, ['Select channel #' num2str(c) ':']) filesep];
        chNames{c} = chPath(length(cellPath{1})+1:end-1);
    end
else
    nCh = length(chNames);
end
for c = 1:nCh
    fprintf('Channel %d name: "%s"\n', c, chNames{c});
end


channels = cell(1,nCh);
for k = 1:nCells
        
    % detect date
    data(k).date = cell2mat(regexp(cellPath{k}, '\d{6}+', 'match'));
    if isempty(date)
        data(k).date = '000000';
    end
    
    % detect frame rate
    fr = regexp(cellPath{k}, '_(\d+)?(.)?\d+s', 'match');
    if ~isempty(fr)
        data(k).framerate = str2double(fr{1}(2:end-1));
    else
        fr = regexp(cellPath{k}, '_\d+ms', 'match');
        if ~isempty(fr)
            data(k).framerate = str2double(fr{1}(2:end-2))/1000;
        end
    end
    
    % assign full channel paths
    for c = 1:nCh
        if ~isempty(chNames{c})
            channels{c} = [cellPath{k} chNames{c} filesep];
        else
            channels{c} = cellPath{k};
        end
        if ~(exist(channels{c}, 'dir')==7)
            channels{c} = [uigetdir(cellPath{k}, ['Select channel #' num2str(c) ':']) filesep];
        end
    end
    data(k).channels = channels;
    data(k).source = channels{1}; % master channel


    % load master channel frames
    tifFiles = dir([data(k).channels{1} '*.tif']);
    data(k).imagesize = size(imread([data(k).channels{1} tifFiles(1).name]));
    data(k).movieLength = length(tifFiles);
    
    if nargin==3
        data(k).markers = markers;
    end
    
    data(k).NA = 1.49;
    data(k).M = 100;
    data(k).pixelSize = 6.7e-6;
    
    fprintf('Loaded: %s\n', cellPath{k});
end