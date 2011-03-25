function imarisApp = viewMovieImaris(movieData3D,iChannel)
%VIEWMOVIEIMARIS opens the input 3D movie for viewing in Imaris
% 
% viewMovieImaris(movieData3D)
% viewMovieImaris(movieData3D,iChannel)
% imarisApp = viewMovieImaris(...);
%
% This function opens Imaris and loads a channel of the input movieData for
% 3D viewing within imaris. Requires that imaris is installed locally.
% 
% Input:
% 
%   movieData3D - The MovieData3D object describing the movie to view.
% 
%   iChannel - The index of the channel to view in imaris. This
%   corresponds to the channel object's location within the channel array
%   in the MovieData.
% 
% Output:
%   
%   imarisApp - The handle to the imaris application. If deleted, imaris
%   will close.
%
%   Additionally, Imaris will be opened and then the images loaded and displayed.
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

%% -------- Init ------ %%

%Get the image file names etc
imageNames = movieData3D.getImageFileNames(iChannel);
imagePath = movieData3D.channels_(iChannel).channelPath_;
nImages = movieData3D.nFrames_;

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
                1,... %Number of Channels
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

%Set fluorescence channel to red, set display range to max
imarisApp.mDataSet.SetChannelColor(0,1,0,0,1); 
volData.SetChannelRange(0,0,2^16-1);

%String for setting frame times. 
oneFrame= movieData3D.timeInterval_/(24*60*60); %Fraction of day corresponding to one frame sec, for datestr.m use
yearString = '2000-01-01 ';%The year/date doesn't seem to matter, so I just use this generic one
msString = '.000'; %The miliseconds portion of the date string

%% ------- Load and Display all Images ----- %%

wtBar = waitbar(0,'Please wait, loading all images...');

for iImage = 1:nImages

    %Load the image
    currIm = stackRead([imagePath filesep imageNames{1}{iImage}]);
    
    %Add it to the imaris scene
    volData.SetDataVolume(currIm,0,iImage-1); %Imaris indexes start at 0
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


    
