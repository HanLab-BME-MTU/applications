function [y1Norm laMean]= normalizeLamellipodiaIntensity(y1,samplesAvg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n m p]=size(samplesAvg);

nWin=n;
nLayer=m;
nTime=p;

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
    cameraBackground=mean(outsideSample);
    
    %for l=1:nWin
    %    minVal(l)=min(samplesAvg(l,3:4,j));
    %end
    %laMean=nanmean(minVal);
    
    laMean(j)=nanmean(samplesAvg(:,3,j));
    %laMean=nanmean(nanmean(samplesAvg(:,6:8,j)));
    
    %laMean=nanmean(nanmean(samplesAvg(:,4:6,j)));
    %laMean=getLamellarBackground(img,maskCell);
    
    for k=1:nWin
        y1Norm(k,j)=(y1(k,j)-cameraBackground)/(laMean(j)-cameraBackground);
    end
end



