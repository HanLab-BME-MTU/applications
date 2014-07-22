function cameraBackground = getCameraBackground()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image');
if ~ischar(filename) || ~ischar(pathname)
     return;
end
imageFileList = getFileStackNames([pathname filesep filename]);
[m nFiles]=size(imageFileList);

[filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Mask');
if ~ischar(filename) || ~ischar(pathname)
     return;
end
maskFileList = getFileStackNames([pathname filesep filename]);
[m nMasks]=size(maskFileList);


for j=1:nFiles
    j
    img=imread(imageFileList{j});
    maskCell=imread(maskFileList{j});
    outsideMask = ~maskCell;
    outsideSample=img(outsideMask(:));
    cameraBackground(j)=mean(outsideSample);
end



