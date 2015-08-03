function [cellFragments] = splitVideo(VIDEO_FILENAME,TIME_WINDOW)
%Splits the tiff video "VIDEO_FILENAME", into fragments of TIME_WINDOW
%number of frames

%reads the video
matStretchedStack = readAndProcessTIFF(VIDEO_FILENAME);

%finds the length of the video
NUM_STACKS = size(matStretchedStack,3);

%initializes the storage cell for video fragments
cellFragments = cell(floor(NUM_STACKS/TIME_WINDOW),1);

%loops through video in TIME_WINDOW steps
prevIdx = 1;%initializes the start position for the fragment
for k = TIME_WINDOW:TIME_WINDOW:NUM_STACKS
    temp_matFragment = matStretchedStack(:,:,prevIdx:k);
    cellFragments{k/TIME_WINDOW} = temp_matFragment;
    
    prevIdx = k;%updates the start position for the next fragment
end
end