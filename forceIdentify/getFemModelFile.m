function [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,varargin)
%getFemModelFile: Get all field boundary files or the most appropriate field boundary file 
%                 for the given image index.
%
% SYNOPSIS: 
%    [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir)
%    [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm)
%
% AUTHOR: Lin Ji.
% DATE  : June 20, 2006.

%Search for field boundary files.
dirList  = dir(femModelDir);
fileList = {dirList(find([dirList.isdir]~=1)).name};

if isempty(fileList)
   femModelFile     = {};
   femModelImgIndex = [];
   return;
end

femModelFileI    = strmatch('femModel',fileList);
femModelFileList = fileList(femModelFileI);
if isempty(femModelFileList)
   femModelFile     = {};
   femModelImgIndex = {};
   return;
end

%Get the corresponding image index associated with the field boundaries.
modelFileImgIndex = [];
modelFileList     = {};
indPat = 'femModel(\d+)\.mat';
for kk = 1:length(femModelFileList)
   indexTokens = regexp(femModelFileList{kk},indPat,'tokens');

   if ~isempty(indexTokens)
      modelFileImgIndex = [modelFileImgIndex str2num(indexTokens{1}{1})];
      modelFileList     = [modelFileList femModelFileList{kk}];
   end
end
if ~isempty(modelFileList)
   %Sort the image index.
   [modelFileImgIndex,sortI] = sort(modelFileImgIndex);
   modelFileList = modelFileList(sortI);
end

if nargin == 1
   for kk = 1:length(modelFileList)
      femModelFile{kk} = [femModelDir filesep modelFileList{kk}];
   end
   femModelImgIndex = modelFileImgIndex;
   return;
end

imgIndex     = varargin{1};
imgIndexForm = varargin{2};

if isempty(modelFileList)
   femModelFile     = '';
   femModelImgIndex = [];
else
   %Look for the most appropriate boundary for the current 'imgIndex'.
   kk = 1;
   while kk <= length(modelFileImgIndex) && modelFileImgIndex(kk) < imgIndex
      kk = kk+1;
   end
   if kk == 1 || modelFileImgIndex(kk) <= imgIndex
      femModelFile = [femModelDir filesep 'femModel' ...
         sprintf(imgIndexForm,modelFileImgIndex(kk)) '.mat'];
      femModelImgIndex = modelFileImgIndex(kk);
   elseif abs(modelFileImgIndex(kk-1)-imgIndex) <= abs(modelFileImgIndex(kk)-imgIndex)
      femModelFile = [femModelDir filesep 'femModel' ...
         sprintf(imgIndexForm,modelFileImgIndex(kk-1)) '.mat'];
      femModelImgIndex = modelFileImgIndex(kk-1);
   end
end

