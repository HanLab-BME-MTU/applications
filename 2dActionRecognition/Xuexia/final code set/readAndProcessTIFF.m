function [matStretchedStack,NUM_STACKS] = readAndProcessTIFF(TIFF_FILENAME)
%Reads and contrast stretches the image TIFF_FILENAME 

%Measure properties of of the image
InfoImage = imfinfo(TIFF_FILENAME);

N_IMAGE = InfoImage(1).Height;
M_IMAGE = InfoImage(1).Width;
NUM_STACKS = length(InfoImage);

%Intializes the storage matrix
matStretchedStack = zeros(N_IMAGE,M_IMAGE,NUM_STACKS,'uint16');

%Loops through the image stack to read them
for i = 1:NUM_STACKS
    %Opens the image frame i
    temp_matImageFrame = imread(TIFF_FILENAME,'Index',i);
    %Contrast stretches the image frame i and stores it
    matStretchedStack(:,:,i) = imadjust(temp_matImageFrame,stretchlim(temp_matImageFrame),[]);
end

%for some reason math fails w/ uint16 so convert to double
matStretchedStack = double(matStretchedStack);
end