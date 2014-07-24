function [ roiMask ] = getSegmentedRois(projList,useFirstImage)
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

    if useFirstImage == 1 
    numberImages = 1 ;
    else 
        numberImages = length(listOfImages); 
    end 
    
for iImg = 1:numberImages
    
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(iImg,1))];
img = double(imread(fileNameIm));

[imL imW] = size(img);

[roiMask] = firstMinAfterFirstMaxSeg(img,0,1,0);

roiMask = roiMask(:,:,1);
%centerRoiXCoord = zeros(2,1);



 
if iImg == 1
roiMaskAvg = roiMask ;
else 
    roiMaskAvg = roiMaskAvg + roiMask; 
end 

end 

[path body no ext] = getFilenameBody(projList(iProj).imDir);

save([path filesep 'roi_1' filesep 'roiSegMaskAvg.mat'], 'roiMaskAvg')

figure; imagesc(roiMaskAvg); 

saveas(gcf,[path filesep 'roi_1' filesep 'avgMask.tif']); 

roiMaskAvg(roiMaskAvg > 1) = 1;



[y1,x1]= ind2sub([imL,imW],find(roiMaskAvg,1));
roiYX = bwtraceboundary(roiMaskAvg(:,:,1),[y1,x1],'N');



save([path filesep 'roi_1' filesep 'roiYXSeg.mat'],'roiYX');


imwrite(roiMaskAvg,[path filesep 'roi_1' filesep 'roiMaskSeg.tif']);

close(gcf)
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













