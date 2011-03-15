function [data] = loadConditionData(condDir, chNames, markers, parameters, movieSelector)
% loadConditionData loads the relevant information for all the data
% available for a specific data condition; this requires a specific
% directory structure and nomenclature (see below)
%
% SYNOPSIS [data] = loadConditionData()
%
% INPUT          condDir : root directory where movies are located
%                chNames : cell array of channel names
%                markers : cell array of fluorescent markers
%           {parameters} : optional, vector of microscope parameters: [NA M pixelSize]
%        {movieSelector} : optional, selector string for movie folders
%
% OUTPUT   data: structure with the fields
%                   .source      : path of the data/movie, location of master channel frames
%                   .channels    : cell array of paths for all channels
%                   .date        : date of the acquisition
%                   .framerate   : frame rate of the movie, in seconds
%                   .imagesize   : dimensions of the movie
%                   .movieLength : length of the movie, in frames
%                   .markers     : cell array of fluorescent marker names
%                   .NA          : numerical aperture of the objective
%                   .M           : magnification of the objective
%                   .pixelSize   : pixel size of the CCD, in meters
%
%
% Francois Aguet, October 2010

if nargin<1
    condDir = [uigetdir(pwd, 'Select the ''condition'' folder:') filesep];
end
if ~strcmp(condDir(end), filesep)
    condDir = [condDir filesep];
end

if nargin<3
    markers = [];
end
if nargin<4 || isempty(parameters)
    parameters = [1.49 100 6.7e-6];
end
if nargin<5 || isempty(movieSelector)
    movieSelector = 'cell';
end

fprintf('Root directory: %s\n', condDir);


% list of experiments for this condition, each containing one or more 'cell' directories
expDir = dirList(condDir);

% if condDir is a 'cell' directory
if ~isempty(regexpi(getDirFromPath(condDir), movieSelector, 'once'))
    cellPath{1} = condDir;
% if expDir are 'cell' directories    
elseif ~isempty(cell2mat(regexpi(arrayfun(@(x) x.name, expDir, 'UniformOutput', false), movieSelector, 'once')))
    cellPath = arrayfun(@(x) [condDir x.name filesep], expDir, 'UniformOutput', false);
else
    cellPath = arrayfun(@(x) arrayfun(@(y) [condDir x.name filesep y.name filesep], dirList([condDir x.name]), 'UniformOutput', false), expDir, 'UniformOutput', false);
    cellPath = vertcat(cellPath{:});
end

% check whether directory names contain 'cell'
valid = cellfun(@(x) ~isempty(regexpi(getDirFromPath(x), movieSelector, 'once')), cellPath);
cellPath = cellPath(valid==1);
nCells = length(cellPath);

% no 'cell' folders are found
if nCells == 0
    error('No data found in directory.');
end

data(1:nCells) = struct('source', [], 'channels', [], 'date', [], 'framerate', [],...
    'imagesize', [], 'movieLength', [], 'markers', [],...
    'framePaths', [], 'maskPaths', []);


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
    for c = 1:nCh
        markers{c} = input(['Enter the fluorescent marker for channel ' num2str(c) ': '], 's');
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
    
    % detect frame
    fr = regexp(cellPath{k}, '_(\d+)?(.)?\d+s', 'match');
    if ~isempty(fr)
        data(k).framerate = str2double(fr{1}(2:end-1));
    else
        fr = regexp(cellPath{k}, '_\d+ms', 'match');
        if ~isempty(fr)
            data(k).framerate = str2double(fr{1}(2:end-2))/1000;
        else
            data(k).framerate = 2; % default: 2s
        end
    end
    
    % assign full channel paths
    framePaths = cell(1,nCh);
    for c = 1:nCh
        if ~isempty(chNames{c})
            channels{c} = [cellPath{k} chNames{c} filesep];
        else
            channels{c} = cellPath{k};
        end
        if ~(exist(channels{c}, 'dir')==7)
            channels{c} = [uigetdir(cellPath{k}, ['Select channel #' num2str(c) ':']) filesep];
        end
        framePaths{c} = arrayfun(@(x) [channels{c} x.name], dir([channels{c} '*.tif*']), 'UniformOutput', false);
    end
    data(k).channels = channels;
    data(k).source = channels{1}; % master channel default
    data(k).framePaths = framePaths;

    
    % load master channel frames
    data(k).imagesize = size(imread(framePaths{1}{1}));
    data(k).movieLength = length(framePaths{1});
    data(k).markers = markers;    
    data(k).NA = parameters(1);
    data(k).M = parameters(2);
    data(k).pixelSize = parameters(3);
    
    fprintf('Loaded: %s\n', cellPath{k});
end