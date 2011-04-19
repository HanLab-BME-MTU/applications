function [ roiMask ] = makeDonutRoi_new(projRoiAll,newProj,bitDepth,pixInt,pixArea)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

roiMaskAll = double(imread([projRoiAll filesep 'roiMask.tif']));
roiYX = load([projRoiAll filesep 'roiYX']);
roiYX = roiYX.roiYX;
imDir = [pwd filesep 'images'];
[listOfImages] = searchFiles('.tif',[],imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm))./((2^bitDepth)-1);

% Make a mask from the image where normalized values above 0.7 
% and thus considered saturated are set to 1 and all non-saturated pixels
% are set to zero. 
% This saturated region will be subtracted from the cell region mask to 
% create a cellular mask which likewise masks out the saturated pixels
% (in our case the pixels corresponding to the centrosome)
img(img>pixInt) = 1;
img(img~=1) = 0;
%roiMaskDonutHole = img;
bw1 = bwmorph(img,'open');
bw2 = bwareaopen(bw1,pixArea);
bw3 = bwmorph(bw2,'thicken');
bw4 = bwmorph(bw3,'bridge');
roiMaskDonutHole = imfill(bw4,'holes');


roiMask = roiMaskAll - roiMaskDonutHole; 
roiMask(roiMask<0) = 0;
%- roiMaskDonutHole;
mkdir(newProj);
imwrite(roiMask,[newProj filesep 'roiMask.tif']);
imwrite(roiMaskDonutHole,[newProj filesep 'roiMaskDonutHole.tif']);
save([newProj filesep 'roiYX'],'roiYX');
end

