function dirList = findProjSubDir(parDir,dirSpec)
%findDir This function finds all the project subdirectories for specified 
%        package common name.
%
% SYNOPSIS : dirList = findDir(parDir,dirSpec)
%    Search in 'parDir' all directories that meet the name specification in
%    'dirSpec' and return a cell array of directories. If 'parDir' does not
%    exit or is empty, an empty array is returned.
%
% INPUT :
%    parDir  : A string that specifies the parent directory. Use '.' for the 
%              current directory.
%    dirSpec : A string for package common name. For example, 'tack' means all
%              directories whose first four characters are 'tack'.

dirList = {};

if isempty(parDir) | ~isdir(parDir)
   return;
end

dirSpecLen = length(dirSpec);

wholeList = dir(parDir);
allDirList = {wholeList(find([wholeList.isdir])).name};

for k = 1:length(allDirList)
   if length(allDirList{k}) >= dirSpecLen
      if strcmp(allDirList{k}(1:dirSpecLen),dirSpec)
         dirList = {dirList{:} allDirList{k}};
      end
   end
end

