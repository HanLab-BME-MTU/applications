function T = calcPreRegByXCorr(inputFileList,path_PreRegT)
% Calculates the shift T between images and a reference frame. But the
% transformation is not performed!
%
% INPUT:
% inputFileList     File list of images that are to be registered. Images
%                   at odd positions in the stack are considered to be the
%                   reference frames.
%
% OUTPUT:           Linear Transformation that can be used with
%                   perfRegInPixStep to register the images.


sortedFileList=getFileListFromFolder(inputFileList);

% read in the reference frame:
refImg   = double(imread(sortedFileList{1}));

% generate the template:
[rows, cols] = size(refImg);
xmin = floor(cols/4);
ymin = floor(rows/4);
width  = floor(cols/2)-1;
height = floor(rows/2)-1;
rect=[xmin ymin width height] ;
template = imcrop(refImg,rect);

% figure 
% imagesc(refImg);
% colormap('gray');
% 
% figure
% imagesc(template);
% colormap('gray');

% find the position of the template in the reference frame. This has to be
% done only once. In principle this could be calculated from above but in
% this way we get around all problems caused by round off effects and so
% on:
selfCorr = normxcorr2(template,refImg);
%figure, surf(selfCorr), shading flat

[~ , imax] = max(abs(selfCorr(:)));
[rowPosInRef, colPosInRef] = ind2sub(size(selfCorr),imax(1));



for i=2:2:length(sortedFileList)
    % Read in the shifted image:
    currImg= double(imread(sortedFileList{i}));
    
    % Find the maximum of the cross correlation between the template and
    % the current image:
    xCorr  = normxcorr2(template,currImg);
    %figure, surf(xCorr), shading flat
    [~ , imax] = max(abs(xCorr(:)));
    [rowPosInCurr, colPosInCurr] = ind2sub(size(xCorr),imax(1));
    
    % The shift is thus:
    rowShift = rowPosInCurr-rowPosInRef;
    colShift = colPosInCurr-colPosInRef;
    
    % and the Transformation is:
    T(round(i/2),:)=-[rowShift colShift];
end
% save it into the current folder
save(path_PreRegT, 'T');