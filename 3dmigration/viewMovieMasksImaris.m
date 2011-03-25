function imarisApp = viewMovieMasksImaris(movieData3D,iChannel)
%VIEWMOVIEMASKSIMARIS views the masks overlain on the images for the input 3D movie in imaris
% 
% viewMovieMasksImaris(movieData3D)
% viewMovieMasksImaris(movieData3D,iChannel)
% imarisApp = viewMovieMasksImaris(...);
%
% This function opens Imaris and loads both the masks and images for a
% specified channel of the input movieData for 3D viewing within imaris.
% This enables visual inspection/validation of the masks. Requires that
% imaris is installed locally and that the specified channel has valid
% masks.
% 
% Input:
% 
%   movieData3D - The MovieData3D object describing the movie to view.
% 
%   iChannel - The index of the channel to view in imaris. This
%   corresponds to the channel object's location within the channel array
%   in the MovieData. This channel must have valid masks.
% 
% Output:
%   
%   imarisApp - The handle to the imaris application. If deleted, imaris
%   will close.
%
%   Additionally, Imaris will be opened and then the images and masks
%   loaded and displayed. The images will be displayed in red, while the
%   masks will be displayed as white.
%   
%Hunter Elliott
%10/2010
%

%% -------- Input -------- %%

if nargin < 1 || isempty(movieData3D) || ~isa(movieData3D,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(iChannel)
    iChannel = 1;
elseif ~isequal(round(abs(iChannel)),iChannel) || numel(iChannel) > 1
    error('The iChannel argument must be a single, positive integer.')
end

%Make sure the movie has been segmented
iSegProc = movieData3D.getProcessIndex('SegmentationProcess3D',1,1);
if isempty(iSegProc)
    error('The input MovieData must have a valid SegmentationProcess3D! Please segment the movie prior to viewing masks!');
end
%Make sure the selected channel has valid masks
if ~(movieData3D.processes_{iSegProc}.checkChannelOutput(iChannel))
    error('The specified channel does not have valid masks! Check the masks and channel argument!')
end

%% -------- Init ------ %%

%Get the image file names etc
imageNames = movieData3D.getImageFileNames(iChannel);
imagePath = movieData3D.channels_(iChannel).channelPath_;
nImages = movieData3D.nFrames_;

%Get mask file names
maskNames = movieData3D.processes_{iSegProc}.getOutMaskFileNames(iChannel);
maskPath = movieData3D.processes_{iSegProc}.outMaskPaths_{iChannel};

%Start imaris and get app handle
imarisApp = imarisStartNew(true);

%Create a blank scene
imarisScene = imarisApp.mFactory.CreateDataContainer;

%Add lighting and frame objects to scene
imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene

%Initialize the volume data
volData = imarisApp.mFactory.CreateDataSet;
volData.Create('eTypeUint16',...
                movieData3D.imSize_(1),... %Image sizes
                movieData3D.imSize_(2),...
                movieData3D.nSlices_,...
                2,... %Number of Channels
                nImages); %Number of timepoints 
            
            
%Add this volume data to the scene
imarisApp.mDataSet = volData;            
%Set the pixel sizes
imarisApp.mDataSet.mExtendMinX = 0;
imarisApp.mDataSet.mExtendMinY = 0;
imarisApp.mDataSet.mExtendMinZ = 0;
imarisApp.mDataSet.mExtendMaxX = imarisApp.mDataSet.mSizeX * movieData3D.pixelSize_;
imarisApp.mDataSet.mExtendMaxY = imarisApp.mDataSet.mSizeY * movieData3D.pixelSize_;
imarisApp.mDataSet.mExtendMaxZ = imarisApp.mDataSet.mSizeZ * movieData3D.zSpacing_;
volData.mUnit = 'nm'; %Set units to nanometers

%Set fluorescence channel to red, set display range to max, opaque
imarisApp.mDataSet.SetChannelColor(0,1,0,0,0); 
volData.SetChannelRange(0,0,2^16-1);
%Set mask channel to white, display range 0 to 1, opaque
imarisApp.mDataSet.SetChannelColor(1,1,1,1,0); 
volData.SetChannelRange(1,0,1);

%String for setting frame times. 
oneFrame= movieData3D.timeInterval_/(24*60*60); %Fraction of day corresponding to one frame sec, for datestr.m use
yearString = '2000-01-01 ';%The year/date doesn't seem to matter, so I just use this generic one
msString = '.000'; %The miliseconds portion of the date string


%% ------- Load and Display all Images ----- %%

wtBar = waitbar(0,'Please wait, loading all images and masks...');

for iImage = 1:nImages

    %Load the image and mask
    currIm = stackRead([imagePath filesep imageNames{1}{iImage}]);
    currMask = tif3Dread([maskPath filesep maskNames{1}{iImage}]);
    
    %Add it to the imaris scene
    volData.SetDataVolume(currIm,0,iImage-1); %Imaris indexes start at 0
    volData.SetDataVolume(uint16(currMask),1,iImage-1); %Imaris indexes start at 0
    
    %Get the time string for this frame
    if iImage == 1
        %Datestr.m doesn't return the last portion if its all zeros...
        secString = '00:00:00';        
    else    
        %Use datestr to convert to hr:min:sec
        secString = datestr(1+(iImage-1)*oneFrame);
        secString = secString(13:end);
    end
    tString = [yearString secString msString];
    volData.SetTimePoint(iImage-1,tString);
    
    
    waitbar(iImage/nImages,wtBar);
    
end


%% ----- Finalization ----- %%

%Adjust the camera so all the data is in view
imarisApp.mSurpassCamera.Fit;

if ishandle(wtBar)
    close(wtBar);
end


    
