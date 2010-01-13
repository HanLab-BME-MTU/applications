function make3DMaskMovieImaris(movieData,iChan)

% 
% make3DMaskMovieImaris(movieData,iChan)
% 
% This function uses Imaris to make a movie overlaying the segmentation on
% the fluorescence for the input movie and channel.
% (REQUIRES IMARIS!)
% 
% Input:
% 
%   movieData - A movieData structure describing a 3D movie as created with
%   setup3DMovieData.m
% 
%   iChan - The index of the channel to make the movie from.
%   Optional. Default is 1.
% 
% 
% Output:
% 
% This function doesn't output shit. It just loads the images and masks
% into imaris and overlays them. I was unable to find a way to force imaris
% to generate a movie from matlab, or to take screenshots of it (the
% saveSnapShot method doesn't seem to work....)
%
% To make the movie, just pick a point of view and then click the little
% red circle in the bottom left corner of the screen.
%
% Hunter Elliott
% 1/2010
%

%% ------- Input ------ %%

%Init / validate the movieData
movieData = setup3DMovieData(movieData);

if nargin < 2 || isempty(iChan)
    iChan = 1;
end

%Make sure it is a 3D movie
if ~ismovie3D(movieData)
    error('This function can only be used with 3D movies! Please create the movieData with setup3dMovieData!')
end

%Make sure the masks have been created
if ~checkMovieMasks(movieData,1);
    error('Problem with masks in input movie! Check mask files and movieData!')
end


%% ----- Init ----- %%


%Get the mask and image file names
maskFiles = getMovieMaskFileNames(movieData,iChan);
imageFiles = getMovieImageFileNames(movieData,iChan);

%Open the imaris application
% disp('Opening Imaris, please wait...');
% imarisApp = actxserver('Imaris.Application');
% imarisApp.mVisible = 1; %Make the interface visible
imarisApp = imarisStartNew(true); %Use this so it doesn't close when execution is finished

%Create the scene in imaris
imarisScene = imarisApp.mFactory.CreateDataContainer;

%Add lighting and frame objects to scene
imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene

imSize = movieData.imSize;

%Initialize the volume data
volData = imarisApp.mFactory.CreateDataSet;
volData.Create('eTypeUint16',...
                imSize(1),... %Image sizes
                imSize(2),...
                imSize(3),...
                2,... %Channel number
                movieData.nImages(iChan)); %Number of timepoints                        
            
%Add the volume data in the scene
imarisApp.mDataSet = volData;   
%Set the color of the volume data (chanID, R, G , B, Alpha)
imarisApp.mDataSet.SetChannelColor(0,1,0,0,.75); %Sets fluorescence channel to red, partially opaque
imarisApp.mDataSet.SetChannelColor(1,0,1,0,.33); %Sets mask channel to green, mostly transparent

%Set range for the channels 
volData.SetChannelRange(1,0,2); %Masks will have values from 0-1, but setting this to two makes the fluorescence easier to see
volData.SetChannelRange(0,0,2^15) %Yeah i know, this assumes 16bit but i dont' care right now

%% ---- Load all the images into imaris ----- %%

nImages = movieData.nImages(iChan);

wtBar = waitbar(0,'Please wait, loading all masks and images...');

for iImage = 1:nImages;

    %Load image & mask for this timepoint
    currIm = double(stkRead(imageFiles{1}{iImage}));    
    currMask = tif3Dread(maskFiles{1}{iImage});
    
    %Normalize the image
    currIm = currIm - min(currIm(:));
    currIm = currIm * round((2^16 ./ max(currIm(:))));
    currIm = uint16(currIm);
    
    %Add the image & mask data to the scene
    volData.SetDataVolume(currIm,0,iImage-1); %Image is chan 0
    volData.SetDataVolume(uint16(currMask),1,iImage-1); %mask is chan 1
        
    waitbar(iImage/nImages,wtBar);
    
end

%Check if the user closed the waitbar, and if not close it
if ishandle(wtBar)
    close(wtBar);
end
