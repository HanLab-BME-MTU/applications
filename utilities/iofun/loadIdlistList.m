function idlistList = loadIdlistList(startDirectory, condition,selectAll,movie)
%LOADIDLISTLIST loads a list of idlists
%
% SYNOPSIS: idlistList = loadIdlistList(startDirectory)
%
% INPUT startDirectory (opt): directory that is given as first choice to
%                             the user
%       condition (opt): string with a condition that the idlist has to
%                        meet to be acceptable. For example
%                        'length(idlist(1).stats.labelcolor == 4'
%                        For alternative ways to state the condition, and
%                        for how to ask the user for input via GUI, see
%                        help checkIdlist
%       selectAll(opt) : true if all idlists found should be selected {0}
%       movie (opt) : 0/[] - nothing
%                     1 - loadStruct for corr/raw
%                     2 - loadStruct for filtered
%
% OUTPUT idlistList: structure array with fields idlist, dataProperties,
%                    name, dirName, dataFileName (includes dirName), slist,
%                    loadStruct
%                    idlist(1).stats.idname is the idlist-type
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
if nargin < 4 || isempty(movie)
    movie = 0;
end
% transform movie to input for loadProjectData
movie = 3-movie;

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

if ~selectAll && size(fileList,1) > 1
    % launch listSelectGUI
    selectIdx = listSelectGUI(fileList(:,1),[],'move');
    % shorten fileList
    fileList = fileList(selectIdx,:);
end

% loop through selected files
nFiles = size(fileList,1);
idlistList(1:nFiles) = struct('idlist',[],'dataProperties',[],'name',[],...
    'dirName',[],'slist',[],'loadStruct',[],'dataFileName',[]);
idlistCt = 1;
for iFile = 1:nFiles
    dataFileName = fullfile(fileList{iFile,2},fileList{iFile,1});
    try
        % load files
        [idlist,dataProperties,projectProperties,slist,loadStruct] = ...
            loadProjectData(fileList{iFile,1},fileList{iFile,2},'last',0,[],movie);

        % check condition
        [goodIdlist,errorMessage,goodTimes] = checkIdlist(idlist,condition);

        % check condition
        if goodIdlist
            % remove bad frames
            idlist(~goodTimes).linklist = [];
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

            % store slist
            idlistList(idlistCt).slist = slist;

            % store loadStruct
            idlistList(idlistCt).loadStruct = loadStruct;

            % store datafile-name
            idlistList(idlistCt).dataFileName = dataFileName;

            % update counter
            idlistCt = idlistCt + 1;
        end
    catch
    end
end
% remove empty entries
idlistList(idlistCt:end) = [];

% remove checkIdlist from memory to clear persistent variables
clear checkIdlist