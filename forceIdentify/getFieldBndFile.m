function [fieldBndFile,fieldBndFileImgIndex] = getFieldBndFile(fieldBndDir,varargin)
%getFieldBndFile: Get all field boundary files or the most appropriate filed boundary file 
%                 for the given image index.
%
% SYNOPSIS: 
%    [fieldBndFile,fieldBndFileImgIndex] = getFieldBndFile(fieldBndDir)
%    [fieldBndFile,fieldBndFileImgIndex] = getFieldBndFile(fieldBndDir,imgIndex,imgIndexForm)
%
% AUTHOR: Lin Ji.
% DATE  : June 20, 2006.

%Search for field boundary files.
dirList  = dir(fieldBndDir);
fileList = {dirList(find([dirList.isdir]~=1)).name};

if isempty(fileList)
   fieldBndFile = {};
   fieldBndFileImgIndex = {};
   return;
end

fieldBndFileI    = strmatch('fieldBnd',fileList);
fieldBndFileList = fileList(fieldBndFileI);
if isempty(fieldBndFileList)
   fieldBndFile         = {};
   fieldBndFileImgIndex = {};
   return;
end

%Get the corresponding image index associated with the field boundaries.
indPat = '\D(\d+)\.mat';
bndFileImgIndex = [];
bndFileList     = {}; %This list does not include 'fieldBnd.mat'.
for kk = 1:length(fieldBndFileList)
   indexTokens = regexp(fieldBndFileList{kk},indPat,'tokens');

   if ~isempty(indexTokens)
      bndFileImgIndex = [bndFileImgIndex str2num(indexTokens{1}{1})];
      bndFileList     = [bndFileList fieldBndFileList{kk}];
   end
end

if ~isempty(bndFileList)
   %Sort the image index.
   [bndFileImgIndex,sortI] = sort(bndFileImgIndex);
   bndFileList = bndFileList(sortI);
end

if nargin == 1
   for kk = 1:length(bndFileList)
      fieldBndFile{kk} = [fieldBndDir filesep bndFileList{kk}];
   end
   fieldBndFileImgIndex = bndFileImgIndex;
   return;
end

imgIndex     = varargin{1};
imgIndexForm = varargin{2};

if isempty(bndFileList)
   fieldBndFile         = '';
   fieldBndFileImgIndex = [];
   return;
else
   %Look for the most appropriate boundary for the current 'imgIndex'.
   kk = 1;
   while kk <= length(bndFileImgIndex) && bndFileImgIndex(kk) < imgIndex
      kk = kk+1;
   end
   if kk == 1 || bndFileImgIndex(kk) <= imgIndex
      fieldBndFile = [fieldBndDir filesep 'fieldBnd' ...
         sprintf(imgIndexForm,bndFileImgIndex(kk)) '.mat'];
      fieldBndFileImgIndex = bndFileImgIndex(kk);
   elseif abs(bndFileImgIndex(kk-1)-imgIndex) <= abs(bndFileImgIndex(kk)-imgIndex)
      fieldBndFile = [fieldBndDir filesep 'fieldBnd' ...
         sprintf(imgIndexForm,bndFileImgIndex(kk-1)) '.mat'];
      fieldBndFileImgIndex = bndFileImgIndex(kk-1);
   end
end

