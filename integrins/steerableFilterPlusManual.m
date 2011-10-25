function mask = steerableFilterPlusManual(image)

%% Output

mask = [];

%% Steerable filter to enhance edges

%get image size
[imSizeX] = size(image,1);
imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));

%use steerable filter to get potential edges
image = double(image);
[~,~,nms] = steerableDetector(image, 3, 1.5);

%% Seeds

%display image
h = figure;
imshow(image,[]);

userEntry = 'y';
xEdge = [];
yEdge = [];
while strcmp(userEntry,'y')
    
    %close figure
    close(h);
    
    %call subfunction that gets edges from manual seeds and NMS image
    [x0,y0] = edgesFromManualSeeds(image,nms);
    
    %add coordinates to existing coordinates and retain only unique elements
    xEdge = [xEdge; x0];
    yEdge = [yEdge; y0];
    xyEdge = [xEdge yEdge];
    xyEdge = unique(xyEdge,'rows');
    xEdge = xyEdge(:,1);
    yEdge = xyEdge(:,2);
    
    if isempty(xEdge)
        disp('No edge points selected ... exiting ...')
        return
    end
    
    %make edge mask
    edgeMask1 = zeros(size(image));
    linIndx = (yEdge-1)*imSizeX + xEdge;
    edgeMask1(linIndx) = 1;
    
    %display edge mask
    imageComp1(:,:,1) = edgeMask1;
    imageComp1(:,:,2) = edgeMask1;
    imageComp1(:,:,3) = imageNorm;
    h = figure;
    imshow(imageComp1,[]);
    
    %ask user if they want to select more seeds
    userEntry = input('Choose more seeds? y/n ','s');
    
end

%% Edge gap closing

userEntry = 'y';
while strcmp(userEntry,'y')
    
    %ask user for dilation radius
    userEntry = input('Choose square edge size to close all gaps in edge. ','s');
    dilationRadius = str2double(userEntry);
    SE = strel('square',dilationRadius);
    edgeMask2 = imclose(edgeMask1,SE);
    
    %display edge mask
    close(h);
    h = figure;
    imageComp2(:,:,1) = edgeMask2;
    imageComp2(:,:,2) = edgeMask2;
    imageComp2(:,:,3) = imageNorm;
    imshow(imageComp2,[]);
    
    %ask user if they want a different dilation radius
    userEntry = input('Choose a different square edge size? y/n ','s');
    
end
close(h);

%% Mask generation

%first fill holes
mask = logical(edgeMask2);
mask = imfill(mask,'holes');

%then ask user to specify points in cell interior for floodfill
figure
imshow(mask);
disp('Please click on points in cell interior to fill')
[yIn,xIn] = getpts;
mask = imfill(mask,[xIn yIn]);
close(h);

%show final mask
SE = strel('square',3);
maskEdge = mask - imerode(mask,SE);
imageComp3(:,:,1) = double(maskEdge);
imageComp3(:,:,2) = double(maskEdge);
imageComp3(:,:,3) = imageNorm;
figure
imshow(imageComp3,[]);

%% ~~~ the end ~~~


%% Sub-function

function [x0,y0,errFlag] = edgesFromManualSeeds(image,nms)

errFlag = 0;

%get image size
[imSizeX,imSizeY] = size(image);

%make composite image for display
nmsNorm = (nms-min(nms(:))) / (max(nms(:))-min(nms(:)));
imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));
imageComp(:,:,1) = nmsNorm;
imageComp(:,:,2) = nmsNorm;
imageComp(:,:,3) = imageNorm;

%allow user to click on image to indicate edges of interest
h = figure;
imshow(imageComp,[]);
userEntry = input('select points? y/n ','s');
x = [];
y = [];
while strcmp(userEntry,'y')
    [yT,xT] = getpts;
    x = [x; xT];
    y = [y; yT];
    userEntry = input('select points again? y/n ','s');
end
close(h);

%starting seed
x0 = round(x);
y0 = round(y);

%keep only points where nms is not zero
linIndx = (y0-1)*imSizeX + x0;
nmsVal = nms(linIndx);
x0 = x0(nmsVal~=0);
y0 = y0(nmsVal~=0);

if isempty(x0);
    disp('No edge points selected ... exiting ...')
    errFlag = 1;
    return
end

%get edges from seeds
lengthIncrease = 1;
while lengthIncrease
    
    %seed after expansion
    x1 = [repmat(x0-1,3,1); repmat(x0,3,1); repmat(x0+1,3,1)];
    y1 = repmat([y0-1; y0; y0+1],3,1);
    indxKeep = find( x1>=1 & x1<=imSizeX & y1>=1 & y1<=imSizeY );
    x1 = x1(indxKeep);
    y1 = y1(indxKeep);
    linIndx = (y1-1)*imSizeX + x1;
    nmsVal = nms(linIndx);
    x1 = x1(nmsVal~=0);
    y1 = y1(nmsVal~=0);
    
    %retain only unique entries
    seed1 = [x1 y1];
    seed1 = unique(seed1,'rows');
    x1 = seed1(:,1);
    y1 = seed1(:,2);
    
    %calculate increase in seed length
    lengthIncrease = length(x1) - length(x0);
    x0 = x1;
    y0 = y1;
    
end

