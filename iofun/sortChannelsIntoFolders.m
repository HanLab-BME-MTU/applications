% Francois Aguet, 101012

function sortChannelsIntoFolders(condDir, channelID, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('condDir', @ischar);
ip.addRequired('channelID');
ip.addParamValue('MovieSelector', 'cell', @ischar);
ip.parse(condDir, channelID, varargin{:});
movieSelector = ip.Results.MovieSelector;
channelID = ip.Results.channelID;
nCh = numel(channelID);

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

if isempty(cellPath)
    error(['No valid movies found in: ' condDir]);
end

% check whether directory names contain 'cell'
val = @(x) ~isempty(x) && x==1;
valid = cellfun(@(x) val(regexpi(getDirFromPath(x), movieSelector, 'once')), cellPath);
cellPath = cellPath(valid==1);
nCells = length(cellPath);

% no 'cell' folders are found
if nCells == 0
    error('No data found in directory.');
end

% loop through cellPath, put channel frames into subfolders
for i = 1:numel(cellPath)
    for c = 1:nCh
        dest = [cellPath{i} channelID{c}];
        movefile([cellPath{i} '*' channelID{c} '*.tif'], dest)
    end    
end
