function idlistList = loadIdlistList(startDirectory, condition,selectAll)
%LOADIDLISTLIST loads a list of idlists
%
% SYNOPSIS: idlistList = loadIdlistList(startDirectory)
%
% INPUT startDirectory (opt): directory that is given as first choice to
%                             the user
%       condition (opt): string with a condition that the idlist has to
%                        meet to be acceptable. For example
%                        'length(idlist(1).stats.labelcolor == 4'
%       selectAll(opt) : true if all idlists found should be selected {0}
%
% OUTPUT idlistList: structure array with fields idlist, dataProperties,
%                    name, dirName
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 22-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 ||isempty(startDirectory)
    startDirectory = pwd;
end
if nargin < 2 || isempty(condition)
    condition = 'true';
end
if nargin < 3 || isempty(selectAll)
    selectAll = false;
end

% select top-directory
if isunix
    oldDiag = getappdata(0,'UseNativeSystemDialogs');
    setappdata(0,'UseNativeSystemDialogs',false)
    topDir = uigetdir(startDirectory);
    setappdata(0,'UseNativeSystemDialogs',oldDiag)
else
    topDir = uigetdir(startDirectory);
end


% find all -data- files
fileList = searchFiles('-data-','log',topDir,1);

if ~selectAll
    % launch listSelectGUI
    selectIdx = listSelectGUI(fileList(:,1),[],'move');
    % shorten fileList
    fileList = fileList(selectIdx,:);
end

% loop through selected files
nFiles = length(fileList);
idlistList(1:nFiles) = struct('idlist',[],'dataProperties',[],'name',[],'dirName',[]);
idlistCt = 1;
for iFile = 1:nFiles
    dataFileName = fullfile(fileList{iFile,2},fileList{iFile,1});
    % load lastResult
    load(dataFileName,'lastResult');
    % load dataProperties
    load(dataFileName,'dataProperties');
    % check if idlist
    if strfind(lastResult,'idlist')
        % if yes: load it
        tmp=load(dataFileName,lastResult);
        idlist = tmp.(lastResult);
        % check condition
        if eval(condition)
            % store idlist
            idlistList(idlistCt).idlist = idlist;
            % store projectName
            dirName = fileList{iFile,2};
            fileSepIdx = strfind(dirName,filesep);
            idlistList(idlistCt).name = dirName(fileSepIdx(end)+1:end);

            % store dataProperties
            idlistList(idlistCt).dataProperties = dataProperties;

            % store dirName
            idlistList(idlistCt).dirName = dirName;

            % update counter
            idlistCt = idlistCt + 1;
        end
    end
end
% remove empty entries
idlistList(idlistCt:end) = [];