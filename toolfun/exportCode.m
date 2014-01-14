%exportCode(masterList, varargin) is a generic function for exporting 
% Matlab functions with their dependent functions

% Francois Aguet, 10/2013

function [ignoreList] = exportCode(masterList, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('masterList');
ip.addOptional('destPath', [], @ischar);
ip.addParamValue('IncludeSources', false, @islogical);
ip.addParamValue('AddList', []); % path to sources, header files etc.
ip.addParamValue('IgnoreList', []);
ip.addParamValue('ExternList', []);
ip.addParamValue('IncludeGPL', true, @islogical);
ip.addParamValue('Verbose', true, @islogical);
ip.parse(masterList, varargin{:});
destPath = ip.Results.destPath;

% this is the master list of all functions to include
if ~iscell(masterList)
    masterList = {masterList};
end

addList = ip.Results.AddList;
if ~isempty(addList)
    if ~iscell(addList)
        addList = {addList};
    end
    addListSource = cellfun(@which, addList, 'unif', 0);
end

if isempty(destPath)
    [~,destPath] = fileparts(masterList{1});
end
if ~strcmpi(destPath(end), filesep)
    destPath = [destPath filesep];
end

nf = numel(masterList);
fctList = cell(1,nf);
toolboxList = cell(1,nf);
for i = 1:nf;
    [fctList{i}, toolboxList{i}] = getFunDependencies(masterList{i});
end
fctList = unique(vertcat(fctList{:}));
toolboxList = unique(vertcat(toolboxList{:}));
[fnames, fpaths, mexNames, mexPaths, sourceNames, sourcePaths, ignoreList] = ...
    parseFilesForExport(fctList, ip.Results.IgnoreList, ip.Results.ExternList);

[~,~] = mkdir(destPath);

% copy core functions
for i = 1:numel(fnames)
    copyfile([fpaths{i} filesep fnames{i}], [destPath fnames{i}]);
end

% copy MEX functions
if numel(mexNames)>0
    mdest = [destPath 'mex' filesep];
    [~,~] = mkdir(mdest);
    for i = 1:numel(mexNames)
        copyfile([mexPaths{i} filesep mexNames{i}], [mdest mexNames{i}]);
    end
end

% copy sources into MEX directory (this does not include external headers!!)
if ip.Results.IncludeSources
    for i = 1:numel(sourceNames)
        copyfile([sourcePaths{i} filesep sourceNames{i}], [mdest sourceNames{i}]);
    end
    for i = 1:numel(addList)
        copyfile(addListSource{i}, [mdest addList{i}]);
    end
end

if ip.Results.IncludeGPL
    gplPath = which('GPL-License.txt');
    copyfile(gplPath, [destPath 'GPL-License.txt']);
end

% set permissions
cmd = ['chmod -R 755 ' destPath];
system(cmd);

if ~isempty(toolboxList) && ip.Results.Verbose
    disp('The package uses the following toolboxes:')
    disp(toolboxList);
end
