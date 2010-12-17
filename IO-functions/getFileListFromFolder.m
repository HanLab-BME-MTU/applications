function [sortedFileList]=getFileListFromFolder(path_folder,pattern)
% get the first bead image:
posExt={'.mat','.tif','.jpg','.TIF'};
entriesDir = dir(path_folder);

% get all relevant filenames
iEntry = 1;
fileList = {};
frameNoList=[];
for i = 1:length(entriesDir)
   if nargin >1 && ~isempty(pattern)
       check=~isempty(strfind(entriesDir(i).name,pattern));
   else
       check=true;
   end   
   if(~entriesDir(i).isdir) && check
      fileList(iEntry) = {strcat(path_folder,filesep,entriesDir(i).name)};
      [~,~,fno,ext]=getFilenameBody(entriesDir(i).name);
      if sum(strcmp(ext,posExt)==1)
        frameNoList(iEntry)=str2double(fno);
        % Here one should check that all images are from the same kind of data.
        iEntry = iEntry + 1;
      end
   end
end

if ~isempty(frameNoList)
    %The outputFileList might be unsorted, this is fixed in the following:
    [~,newIndx] = sort(frameNoList);
    sortedFileList=fileList(newIndx);
else
    sortedFileList=[];
end

