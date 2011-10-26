function [mask,filterSigma,closureRadius] = steerableFilterPlusManual(image)
%STEERABLEFILTERPLUSMANUAL is an interactive tool that allows the user to find cell edges based on intensity gradients
%
%SYNOPSIS [mask,filterSigma,closureRadius] = steerableFilterPlusManual(image)
%
%INPUT  image        : Cell image to be segmented.
%
%OUTPUT mask         : Mask (1 inside cell, 0 outside).
%       filterSigma  : Gaussian standard deviation used in line filter.
%       closureRadius: Radius used to close gaps in edge.
%
%Khuloud Jaqaman, October 2011

%% Output

mask = [];

%% Input

%get image size and normalize it for later display
image = double(image);
[imSizeX] = size(image,1);
imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));

%% Steerable filter to enhance edges

userEntry = 'y';
while strcmp(userEntry,'y')
    
    %get filter sigma from user
    userEntry = input('Filter sigma: ','s');
    filterSigma = str2double(userEntry);
    
    %use steerable filter to get potential edges
    [~,~,nms] = steerableDetector(image, 3, filterSigma);
    
    %make composite image for display
    nmsNorm = (nms-min(nms(:))) / (max(nms(:))-min(nms(:)));
    imageComp(:,:,2) = nmsNorm;
    imageComp(:,:,1) = 0;
    imageComp(:,:,3) = imageNorm;
    
    %display composite image
    h = figure;
    imshow(imageComp,[]);
    
    %ask user if they want to use a different sigma
    userEntry = input('Choose a different filter sigma? y/n ','s');
    close(h);
    
end
    
%% Seeds

%display image
h = figure;
imshow(image,[]);

userEntry = 'y';
xEdge = [];
yEdge = [];
edgeMask1 = zeros(size(image));
while strcmp(userEntry,'y')
    
    %close figure
    close(h);
    
    %call subfunction that gets edges from manual seeds and NMS image
    [x0,y0] = edgesFromManualSeeds(imageNorm,nms,edgeMask1);
    
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
    imageComp(:,:,1) = edgeMask1;
    imageComp(:,:,2) = 0;
    h = figure;
    imshow(imageComp,[]);
    
    %ask user if they want to select more seeds
    userEntry = input('Select more edge points? y/n ','s');
    
end

%% Edge gap closing

userEntry = 'y';
while strcmp(userEntry,'y')
    
    %ask user for closure radius
    userEntry = input('Closure radius to close gaps in edge: ','s');
    closureRadius = str2double(userEntry);
    SE = strel('disk',closureRadius);
    edgeMask2 = imclose(edgeMask1,SE);
    edgeMask2 = bwmorph(edgeMask2,'bridge');
    
    %display edge mask
    close(h);
    h = figure;
    imageComp(:,:,1) = edgeMask2;
    imageComp(:,:,2) = 0;
    imshow(imageComp,[]);
    
    %ask user if they want a different dilation radius
    userEntry = input('Choose another closure radius? y/n ','s');
    
end
close(h);

%% Mask generation

%first fill holes
mask = logical(edgeMask2);
mask = imfill(mask,'holes');

%then ask user to specify points in cell interior for floodfill
figure
imshow(mask);
disp('Click on points in cell interior to fill the whole cell.')
close(h);
[yIn,xIn] = getpts;
mask = imfill(mask,round([xIn yIn]));

%some final polishing of final mask ...
mask = imfill(mask,'holes');
SE = strel('disk',5);
mask = imerode(mask,SE);
mask = imdilate(mask,SE);

%show final mask
SE = strel('square',3);
maskEdge = mask - imerode(mask,SE);
imageComp(:,:,1) = double(maskEdge);
imageComp(:,:,2) = 0;
figure
imshow(imageComp,[]);

%% ~~~ the end ~~~


%% Sub-function

function [x0,y0,errFlag] = edgesFromManualSeeds(image,nms,edgeMask1)

errFlag = 0;

%get image size
[imSizeX,imSizeY] = size(image);

%make composite image for display
nmsNorm = (nms-min(nms(:))) / (max(nms(:))-min(nms(:)));
nmsNorm(find(edgeMask1==1)) = 0; %#ok<FNDSB>
imageComp(:,:,1) = edgeMask1;
imageComp(:,:,2) = nmsNorm;
imageComp(:,:,3) = image;

%allow user to click on image to indicate edges of interest
h = figure;
imshow(imageComp,[]);
disp('Select edge points.')
[yT,xT] = getpts;
close(h);

%starting seed
x0 = round(xT);
y0 = round(yT);

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

