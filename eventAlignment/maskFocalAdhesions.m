function maskFA = maskFocalAdhesions(minSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image');
   
if ~ischar(filename) || ~ischar(pathname)
     return;
end
   
inputFileList = getFileStackNames([pathname filesep filename]);

[m nFiles]=size(inputFileList);

for j=1:nFiles
    j
    img{j}=imread(inputFileList{j});
    maskFA{j}=blobSegmentThreshold(img{j},minSize,0);
    img{j}=double(img{j});
    img{j}=img{j}.*maskFA{j};
end
writeImageFiles(img)


