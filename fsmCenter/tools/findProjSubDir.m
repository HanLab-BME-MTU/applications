function dirList = findProjSubDir(parDir,dirSpec)
%findDir This function finds all the project subdirectories for specified 
%        package common name.
%
% SYNOPSIS : dirList = findDir(parDir,dirSpec)
%    Search in 'parDir' all directories that meet the name specification in
%    'dirSpec' and return a cell array of directories.
%
% INPUT :
%    parDir  : A string that specifies the parent directory. If it is empty,
%              the current directory is searched.
%    dirSpec : A string for package common name. For example, 'tack' means all
%              directories whose first four characters are 'tack'.

if isempty(parDir)
   parDir = pwd;
end

if ~isdir(parDir)
   error([parDir 'does not exist.']);
end

dirSpecLen = length(dirSpec);

wholeList = dir(parDir);
allDirList = {wholeList(find([wholeList.isdir])).name};

dirList = {}:
for k = 1:length(allDirList)
   if length(allDirList{k}) >= dirSpecLen
      if strcmp(allDirList{k}(1:dirSpecLen),dirSpec)
         dirList = {dirList{:} allDirList{k}};
      end
   end
end

