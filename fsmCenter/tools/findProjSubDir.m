function dirList = findProjSubDir(parDir,dirSpec)
%findDir This function finds all the project subdirectories for specified 
%        package common name.
%
% SYNOPSIS : dirList = findProjSubDir(parDir,dirSpec)
%    Search in 'parDir' all directories that meet the name specification in
%    'dirSpec' and return a cell array of directories. If 'parDir' does not
%    exit or is empty, an empty array is returned.
%    Alternatively, parDir can already be a structure returned by the DIR command. 
%    Since MATLAB's DIR is slow, this can speed up functions which call FINDPROJSUBDIR
%    many times.
%    In this case, parDir is a structure with fields:
%       name
%       date
%       bytes
%       isdir
%
% INPUT :
%    parDir  : A string that specifies the parent directory. Use '.' for the 
%              current directory.
%    dirSpec : A string for package common name. For example, 'tack' means all
%              directories whose first four characters are 'tack'.

dirList = {};

if isempty(parDir)
   return;
end

% Default
readDir=1;

if isstruct(parDir)
    % Check fields
    fields=fieldnames(parDir);
    if ~ ( strcmp(char(fields(1)),'name') & ...
            strcmp(char(fields(2)),'date') & ...
            strcmp(char(fields(3)),'bytes') & ...
            strcmp(char(fields(4)),'isdir') )
        error('Invalid structure parDir (should be a structure as returned by the MATLAB command DIR.)');
    end
    readDir=0;
else
    % Check that parDir is a string pointing to an existing directory
    if ~strcmp(class(parDir),'char')
        error('parDir is not a valid string.');
    else
        if ~isdir(parDir)
            error('parDir is not pointing to a valid directory.');
        end
    end
end

dirSpecLen = length(dirSpec);

if readDir==1
    parDir = dir(parDir);
end
allDirList = {parDir(find([parDir.isdir])).name};

for k = 1:length(allDirList)
   if length(allDirList{k}) >= dirSpecLen
      if strcmp(allDirList{k}(1:dirSpecLen),dirSpec)
         dirList = {dirList{:} allDirList{k}};
      end
   end
end

