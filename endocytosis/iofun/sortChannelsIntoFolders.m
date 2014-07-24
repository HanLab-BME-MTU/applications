% Francois Aguet, 101012

function sortChannelsIntoFolders(rootDir, channelID, channelDir)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('condDir', @ischar);
ip.addRequired('channelID');
ip.addRequired('channelDir');
ip.parse(rootDir, channelID, channelDir);

% get tree of directories from rootDir
dirlist = recursiveDir(rootDir);
nc = numel(channelID);

% parse each dir for TIF files that contain channel IDs
for i = 1:numel(dirlist)
    files = [dir([dirlist{i} '*.tif*']) dir([dirlist{i} '*.TIF*'])];
    if ~isempty(files)
        files = arrayfun(@(x) x.name, files, 'unif', 0);
        cfiles = cell(1,nc);
        for c = 1:nc
            idx = ~cellfun(@isempty, regexpi(files, channelID{c}));
            cfiles{c} = files(idx);
        end
        % if files found for all channels
        if ~any(cellfun(@isempty, cfiles))
            for c = 1:nc
                % create channel directory
                dest = [dirlist{i} channelDir{c} filesep];
                fprintf('Moving files to %s\n', dest);
                [~,~] = mkdir(dest);
                cfiles{c} = cellfun(@(x) [dirlist{i} x], cfiles{c}, 'unif', 0);
                movefile([dirlist{i} '*' channelID{c} '*'], dest);
            end
        end
    end
end



function p = recursiveDir(d)

p = {};% path to be returned

% Generate path based on given root directory
files = dir(d);
if isempty(files)
  return
end

isdir = logical(cat(1,files.isdir));
dirs = files(isdir);

% Remove invisible
idx = arrayfun(@(i) ~strcmp(i.name(1), '.'), dirs);
dirs = dirs(idx);

for i = 1:length(dirs)
    tmp = [d dirs(i).name filesep];    
    p = [p tmp recursiveDir(tmp)]; %#ok<AGROW>
end



% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('condDir', @ischar);
% ip.addRequired('channelID');
% ip.addParamValue('MovieSelector', 'cell', @ischar);
% ip.parse(condDir, channelID, varargin{:});
% movieSelector = ip.Results.MovieSelector;
% channelID = ip.Results.channelID;
% nCh = numel(channelID);
% 
% % list of experiments for this condition, each containing one or more 'cell' directories
% expDir = dirList(condDir);
% 
% % if condDir is a 'cell' directory
% if ~isempty(regexpi(getDirFromPath(condDir), movieSelector, 'once'))
%     cellPath{1} = condDir;
%     % if expDir are 'cell' directories
% elseif ~isempty(cell2mat(regexpi(arrayfun(@(x) x.name, expDir, 'UniformOutput', false), movieSelector, 'once')))
%     cellPath = arrayfun(@(x) [condDir x.name filesep], expDir, 'UniformOutput', false);
% else
%     cellPath = arrayfun(@(x) arrayfun(@(y) [condDir x.name filesep y.name filesep], dirList([condDir x.name]), 'UniformOutput', false), expDir, 'UniformOutput', false);
%     cellPath = vertcat(cellPath{:});
% end
% 
% if isempty(cellPath)
%     error(['No valid movies found in: ' condDir]);
% end
% 
% % check whether directory names contain 'cell'
% val = @(x) ~isempty(x) && x==1;
% valid = cellfun(@(x) val(regexpi(getDirFromPath(x), movieSelector, 'once')), cellPath);
% cellPath = cellPath(valid==1);
% nCells = length(cellPath);
% 
% % no 'cell' folders are found
% if nCells == 0
%     error('No data found in directory.');
% end
% 
% % loop through cellPath, put channel frames into subfolders
% for i = 1:numel(cellPath)
%     for c = 1:nCh
%         dest = [cellPath{i} channelID{c}];
%         movefile([cellPath{i} '*' channelID{c} '*.tif'], dest)
%     end    
% end
