function [ roiMask ] = getSegmentedRois(projList)
% Currently there is no automatic cell-edge segmentation in the plusTipTracker. In order to 
% perform analysis relative to cell edge one therefore needs to call a 
% separate function to segment the image
% this function will obtain a mask corresponding to the first image 
% of the respective movie for all movies in a give projList. 
% (Note uses same segmentation function as the windowing package
% 'firstMinAfterFirstMaxSeg')
% Used currently (03/11) for Mijung's micropattern experiments.
% MB 03/11 



for iProj = 1:length(projList)
    
    
[listOfImages] = searchFiles('.tif',[],projList(iProj).imDir,0);

fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));

[imL imW] = size(img);

[roiMask] = firstMinAfterFirstMaxSeg(img,0,1,0);

roiMask = roiMask(:,:,1);
%centerRoiXCoord = zeros(2,1);

[path body no ext] = getFilenameBody(projList(iProj).imDir);

[y1,x1]= ind2sub([imL,imW],find(roiMask,1));
roiYX = bwtraceboundary(roiMask(:,:,1),[y1,x1],'N');

save([path filesep 'roi_1' filesep 'roiYXSeg.mat'],'roiYX');

imwrite(roiMask,[path filesep 'roi_1' filesep 'roiMaskSeg.tif']);

%for i = 1:2
%stats = regionprops(bwlabel(roiMask(:,:,i)),'centroid');
%centerRoiXCoord(i,1) = stats.Centroid(1,1);
%end 

%[value, idxMax] = max(centerRoiXCoord);
%[value, idxMin] = min(centerRoiXCoord);


%imwrite(roiMask(:,:,idxMin),[path filesep 'roi_3' filesep 'roiMaskSeg.tif']);
%roiMaskLeft = roiMask(:,:,idxMin);
%[y1,x1]= ind2sub([imL,imW],find(roiMaskLeft,1));
%roiYX = bwtraceboundary(roiMaskLeft,[y1,x1],'N');

%save([path filesep 'roi_3' filesep 'roiYXSeg.mat'],'roiYX');


%imwrite(roiMask(:,:,idxMax),[path filesep 'roi_4' filesep 'roiMaskSeg.tif']);

%[y1,x1]= ind2sub([imL,imW],find(roiMask(:,:,idxMax),1));
%roiYX = bwtraceboundary(roiMask(:,:,idxMax),[y1,x1],'N');

%save([path filesep 'roi_4' filesep 'roiYXSeg.mat'],'roiYX');

end 
end 













