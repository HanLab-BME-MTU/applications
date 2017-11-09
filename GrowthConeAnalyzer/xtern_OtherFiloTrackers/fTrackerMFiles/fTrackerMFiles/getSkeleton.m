function [skeleton, centroid, singleCone,closedCone,imageBackground,dilatedCone,OtherImage]=getSkeleton(imagen, threshold, pixelsDilate, ...
    show, skelFill, diskSize, closeSize, byThresholding,EdgeStrongness)

% Outputs the skeleton of an image after some magic touches

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


% Kernel for imcolse
disco = strel('disk', closeSize);

if byThresholding==1   
    %% subtract background
    imageBackground = imopen(imagen,strel('disk',diskSize));
    newImage = imsubtract(imagen,imageBackground);
    newImage = newImage+max(max(imageBackground))-1; 

    %% First threshold
    binaryImage=zeros(size(imagen));
    binaryImage(newImage>threshold*mean(mean(newImage))) = 1;
    binaryImage=bwmorph(binaryImage, 'clean');               %% Clean single pixels
    singleCone=makeSingleObjectImage(binaryImage);           %% Remove stuff outside the biggest object   
    dilatedCone=bwmorph(singleCone, 'dilate', pixelsDilate); %% Dilate
    
    filteredCone = imfilter(dilatedCone, ones(pixelsDilate)/sum(sum(ones(pixelsDilate)))); %% Convolve with a square
    filteredCone(filteredCone<0.7)=0;    
    
    closedCone=imclose(filteredCone,disco); %% Enlarge with imclose
    if skelFill==1
        filledCone = imfill(closedCone,'holes');
    else
        filledCone = closedCone;
    end
    OtherImage{1} = binaryImage;
    OtherImage{2} = filteredCone;
else
    imageBackground = 0*imagen+1;
    newImage = imagen;
    
    imageBackground = imopen(imagen,strel('disk',diskSize));
    newImage = imsubtract(imagen,imageBackground);
    
    [edgesImage,thresEDGE] = edgeCHANGED(newImage,'zerocross',EdgeStrongness,'log');

%     filteredCone=imfilter(edgesImage, ones(pixelsDilate)/sum(sum(ones(pixelsDilate))));     %% Convolve with a square
%     filteredCone(filteredCone<0.7)=0;

    closedCone  = imclose(edgesImage,disco);                      %% Enlarge with imclose
    dilatedCone = bwmorph(closedCone, 'dilate', pixelsDilate);    %% Dilate
    singleCone  = makeSingleObjectImage(dilatedCone);
    if skelFill==1
        filledCone = imfill(singleCone,'holes');
    else
        filledCone = singleCone;
    end
    OtherImage{1}= edgesImage;
    OtherImage{2}= filledCone;
end
    centroid = floor(getCenterOfMass(newImage, closedCone));  %% WHY WITH CLOSEDCONE and NOT with singleCONE ???
    skeleton = bwmorph(filledCone, 'thin', Inf);
if show==1
    h=figure('Name', 'Background Substracted');h=imshow(imageBackground);
    h=figure('Name', 'Image Analized');h=imshow(newImage);
    if byThresholding==1
        h = figure('Name', 'Binary Image');h=imshow(binaryImage);
        h = figure('Name', 'Convolved');h=imshow(filteredCone);
    else
        h=figure('Name', 'Edges Found');h=imshow(edgesImage);
        h=figure('Name', 'Image Dilated');h=imshow(dilatedCone);
    end 
    h = figure('Name', 'Dilated Image');h=imshow(dilatedCone);
    h = figure('Name', 'Filled In');h=imshow(filledCone);
    h = figure('Name', 'Single Cone');h=imshow(singleCone);
    h = figure('Name', 'Closed Cone');h=imshow(closedCone);
end