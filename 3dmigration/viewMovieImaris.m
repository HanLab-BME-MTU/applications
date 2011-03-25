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

nChan = 1;

%Get the image file names etc
imageNames = movieData3D.getImageFileNames(iChannel);
imagePaths{1} = movieData3D.channels_(iChannel).channelPath_;
nImages = movieData3D.nFrames_;
chanCols = [1 0 0];
chanRange = [0 2e16-2];

%Start imaris and get app handle
imarisApp = imarisStartNew(nargout==0);

%Create a blank scene
imarisScene = imarisApp.mFactory.CreateDataContainer;

%Add lighting and frame objects to scene
imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene

%Check if we have masks or skeletons as these will be additional channels
chanNames = {'Fluorescence','Masks','Skeleton'};
               
%Check for masks 
iSegProc = movieData3D.getProcessIndex('SegmentationProcess3D',1,1);
if ~isempty(iSegProc) && movieData3D.processes_{iSegProc}.checkChannelOutput(iChannel)
    disp('Masks found - displaying as additional channel.')
    nChan = nChan + 1;
    imagePaths{nChan} = movieData3D.processes_{iSegProc}.outMaskPaths_{iChannel};
    imageNames{nChan} = movieData3D.processes_{iSegProc}.getOutMaskFileNames(iChannel);
    imageNames{nChan} = imageNames{nChan}{1};%De-cell this element
    chanCols = vertcat(chanCols,[1 1 1]);
    chanRange = vertcat(chanRange,[0 2]);%We make the range go to slightly above 1 so the masks are transparent
end

%Check for skeletons
iSkelProc = movieData3D.getProcessIndex('SkeletonizationProcess',1,1);
if ~isempty(iSkelProc) && movieData3D.processes_{iSkelProc}.checkChannelOutput(iChannel)
    disp('Masks found - displaying as additional channel.')
    nChan = nChan + 1;
    imagePaths{nChan} = movieData3D.processes_{iSkelProc}.outImagePaths_{iChannel};
    imageNames{nChan} = movieData3D.processes_{iSegProc}.getOutImageFileNames(iChannel);
    imageNames{nChan} = imageNames{nChan}{1};%De-cell this element
    chanCols = vertcat(chanCols,[0 0 1]);
    chanRange = vertcat(chanRange,[0 1]);
end


%Initialize the volume data
volData = imarisApp.mFactory.CreateDataSet;
volData.Create('eTypeUint16',...
                movieData3D.imSize_(1),... %Image sizes
                movieData3D.imSize_(2),...
                movieData3D.nSlices_,...
                nChan,... %Number of Channels
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

%String for setting frame times. 
oneFrame= movieData3D.timeInterval_/(24*60*60); %Fraction of day corresponding to one frame sec, for datestr.m use
yearString = '2000-01-01 ';%The year/date doesn't seem to matter, so I just use this generic one
msString = '.000'; %The miliseconds portion of the date string

%% ------- Load and Display all Images ----- %%

wtBar = waitbar(0,'Please wait, loading all images...');

for iImage = 1:nImages

    for iChan = 1:nChan
    
        if iImage == 1
            %Set channel color and range
            imarisApp.mDataSet.SetChannelColor(iChan-1,...
                                        chanCols(iChan,1),...
                                        chanCols(iChan,2),...
                                        chanCols(iChan,3),0);                                 
            imarisApp.mDataSet.SetChannelRange(iChan-1,...
                                            chanRange(iChan,1),...
                                            chanRange(iChan,2));
                                        
            imarisApp.mDataSet.SetChannelName(iChan-1,chanNames{iChan});
            
        end
        
        %Load the image
        if iChan == 1
            currIm = stackRead([imagePaths{iChan} filesep imageNames{iChan}{iImage}]);
        else
            %Stackread doesn't support the bitpacking compression of
            %binary tifs
            currIm = uint16(tif3Dread([imagePaths{iChan} filesep imageNames{iChan}{iImage}]));
        end
    
        %Add it to the imaris scene
        imarisApp.mDataSet.SetDataVolume(currIm,iChan-1,iImage-1); %Imaris indexes start at 0
        
    end
        
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
    imarisApp.mDataSet.SetTimePoint(iImage-1,tString);

    waitbar(iImage/nImages,wtBar);
    
end

%Adjust the camera so all the data is in view
imarisApp.mSurpassCamera.Fit;

if ishandle(wtBar)
    close(wtBar);
end






    
