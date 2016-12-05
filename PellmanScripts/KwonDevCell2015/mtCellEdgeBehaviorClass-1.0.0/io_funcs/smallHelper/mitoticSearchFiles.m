function [listOfFiles,tokenList] = searchFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode,fullName)
%searchFiles is an utility to search for files containing a specific string
%
%SYNOPSIS listOfFiles = searchFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode,fullName)
%
%INPUT    includeString: string contained in the filenames you are looking
%                               for (can be regular expression; is case
%                               insensitive)
%         excludeString (opt): string not contained in the filenames you
%                               are looking for (can be regular expression;
%                               is case insensitive)
%         directory (opt): directory to search. if empty, current directory is searched (default)
%                               if 'ask', program asks for directory
%         includeSubDirectories (opt): whether to search subdirectories or not (0/{1})
%         selectionMode (opt): which file(s) to select if there are several files matching the
%                               includeString within one directory
%                              {'all'}-all; 'new'-newest; 'old'-oldest; 'GUI'-open GUI to select one file
%         fullName (opt): if 1, listOfFiles will be a nFiles-by-1 cell
%                               array with full path. Default: 0
%
%OUTPUT  listOfFiles: cell array with {[fileName] [directoryName]}
%        tokenList  : if includeString is a regular expression with tokens
%                       (see regexp for help), tokens are returned in a
%                       tokenList (cell array of strings)
%                           tmp1 = ...
%                           regexp(listOfFiles(:,1),includeString,'tokens')
%                           tmp2 = cat(1,tmp1{:});
%                           tokenList = cat(1,tmp2{:});
%                       if no tokens are indicated in the input, tokenList will
%                       be empty.
%                       Assuming that all the tokens are numbers,
%                       converting tokenList into a matrix of doubles can
%                       be done as follows:
%                           tokenListAsDoubles = str2double(tokenList);
%
%c: 7-03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---test input---
if nargin<1||isempty(includeString)
    error('not enough input arguments or empty includeString')
end

%includeString
if ~ischar(includeString)
    error('includeString has to be a string')
end

%excludeString
if nargin>1 && ~isempty(excludeString)
    if ~ischar(excludeString)
        error('excludeString has to be a string')
    end
else
    excludeString = [];
end

%directory (ask if necessary)
if nargin>2
    if ~exist('directory','var') || isempty(directory)
        directory = pwd;
    elseif strcmp(directory,'ask')

        % on linux, uigetdir requires the selection of a file, which is
        % rater inconvenient, if there are no files in the directory. Use
        % workaround

        if isunix
            % remember system dialog setting
            sysDialogState = getappdata(0,'UseNativeSystemDialogs');
            setappdata(0,'UseNativeSystemDialogs',false)
        end

        directory = uigetdir(pwd,'select a directory to search');

        if isunix
            % reset system dialog setting
            setappdata(0,'UseNativeSystemDialogs',sysDialogState)
        end

        if directory == 0
            error('searchFiles aborted by user')
        end
    elseif ~isdir(directory)
        error([directory,' is not a valid directory!'])
    end
else
    directory = pwd;
end

%includeSubDirectories
if nargin<4||isempty(includeSubDirectories)
    includeSubDirectories = 1;
end

%selectionMode
if nargin<5||isempty(selectionMode)
    selectionMode = 'all';
else
    if ~any([strcmp(selectionMode,'new'),...
            strcmp(selectionMode,'old'),...
            strcmp(selectionMode,'GUI'),...
            strcmp(selectionMode,'all')])
        error('wrong selectionMode')
    end
end

% fullfile
if nargin<6||isempty(fullName)
    fullName = false;
end

%---end test input---


%---collect files---

%look for directories and wanted files in current directory, then check all subdirs

%init variables

dirs2checkInitLength = 100;
dirs2checkLength     = 100;
dirs2check           = cell(100,1);

dirs2check{1} = directory;
dirs2checkCt  = 1;
topDir2check  = 1;

listOfFilesInitLength = 100;
listOfFilesLength     = 100;
listOfFiles           = cell(100,2-fullName);
listOfFilesCt         = 0;


while ~isempty(dirs2check) && ~isempty(dirs2check{topDir2check})
    %init/empty var
    pathCell = {};

    %read currentDir
    currentDir = dirs2check{topDir2check};

    %list all files of current directory
    currentDirList = dir(currentDir);

    %look for subdirectories and add to dirs2check
    isDirList = cat(1,currentDirList.isdir);
    if length(isDirList)>2
        subDirIdx = find(isDirList(3:end))+2;
        for i=1:length(subDirIdx)

            % we will add a directory - count how many this will make
            dirs2checkCt = dirs2checkCt + 1;

            % make sure that the dirList is long enough - make sure we do
            % not get problems with topDir2check
            if dirs2checkCt + 1 > dirs2checkLength
                tmpDirs2check = dirs2check;
                newDirs2checkLength = dirs2checkLength + dirs2checkInitLength;
                dirs2check = cell(newDirs2checkLength,1);
                dirs2check(1:dirs2checkLength) = tmpDirs2check;
                dirs2checkLength = newDirs2checkLength;
            end

            dirs2check{dirs2checkCt} = [currentDir,filesep,currentDirList(subDirIdx(i)).name];
        end
    end %if length(isDirList)>2

    %look for files in current directory and store them
    newFiles = chooseFile(includeString,currentDir,selectionMode,excludeString);
    if ~isempty(newFiles)
        if ~iscell(newFiles)
            newFiles = cellstr(newFiles);
        end
        [pathCell{1:size(newFiles,1)}] = deal(currentDir);

        % we will add a file - count how many this will make
        numNewFiles = length(newFiles);
        listOfFilesCtStart = listOfFilesCt + 1;
        listOfFilesCt = listOfFilesCt + numNewFiles;

        % make sure that the dirList is long enough
        while listOfFilesCt > listOfFilesLength
            tmpListOfFiles = listOfFiles;
            newListOfFilesLength = listOfFilesLength + listOfFilesInitLength;
            listOfFiles = cell(newListOfFilesLength,2-fullName);
            listOfFiles(1:listOfFilesLength,:) = tmpListOfFiles;
            listOfFilesLength = newListOfFilesLength;
        end

        if fullName
            for f=1:numNewFiles
                listOfFiles{listOfFilesCtStart-1+f} = fullfile(currentDir,newFiles{f});
            end
        else
        listOfFiles(listOfFilesCtStart:listOfFilesCt,:) = [newFiles, pathCell'];
        end

    end %if ~isempty(newFiles)

    % move on to next directory
    topDir2check = topDir2check + 1;

    %check wheter we want to look at subDirs
    if ~includeSubDirectories
        dirs2check = {};
    end


end %while ~isempty(dirs2check)

%---end collect files---

% remove placeholders
listOfFiles(listOfFilesCt+1:end,:)=[];



%----------- find tokens --------------
if nargout > 1
    % we need to make a tokenList
    if any(findstr(includeString,'(')) && any(findstr(includeString,')'))
        % only call regexp if there are tokens at all

        tmp = regexp(listOfFiles(:,1),includeString,'tokens');
        tmp = cat(1,tmp{:});
        if ~isempty(tmp)
            tokenList = cat(1,tmp{:});
        else
            tokenList = [];
        end
    else
        tokenList = [];
    end
end