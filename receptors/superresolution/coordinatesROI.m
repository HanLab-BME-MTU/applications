function movieInfoCropped = coordinatesROI(movieInfo,doPlot)
%CROPDETECTIONRESULTS keeps only coordinates within a region of interest
%
%SYNPOSIS movieInfoCropped = coordinatesROI(movieInfo,doPlot)
%
%INPUT  movieInfo   : Coordinates in the format of the output of
%                     detectSubResFeatures2D_StandAlone.
%       doPlot      : 1 to make a plot in the end, 0 otherwise.
%                     Optional. Default: 0.
%
%OUTPUT movieInfoCropped : Coordinates in region of interest, in same
%                          format as input.
%
%Khuloud Jaqaman, August 2011

%% Input
if nargin < 1
    disp('Wrong number of input arguments')
    return
end

if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

%% Output
movieInfoCropped = movieInfo;

%% Crop

%get number of frames
numFrames = length(movieInfo);

%extract all coordinates
xCoord = vertcat(movieInfo.xCoord);
yCoord = vertcat(movieInfo.yCoord);

%get image for plotting and cropping
[fName,dirName] = uigetfile('*.tif','Choose image to help with cropping');
currentImage = double(imread(fullfile(dirName,fName)));

%display image with coordinates overlaid
currentImage = (currentImage-min(currentImage(:)))/(max(currentImage(:))-min(currentImage(:)));
h = figure;
imshow(currentImage);
hold on
plot(xCoord(:,1),yCoord(:,1),'g.')

%crop and get range of interest
[~,rect] = imcrop(h);
close all
rect = round(rect);
imageRange = [rect(2) rect(2)+rect(4); rect(1) rect(1)+rect(3)];

%go over all frames and keep only coordinates withing image range of
%interest
for iFrame = 1 : numFrames
    
    xCoord = movieInfoCropped(iFrame).xCoord;
    yCoord = movieInfoCropped(iFrame).yCoord;
    amp    = movieInfoCropped(iFrame).amp;
    
    indxKeep = find( xCoord(:,1) >= (imageRange(2,1)-0.5) & ...
        xCoord(:,1) <= (imageRange(2,2)+0.5) & ...
        yCoord(:,1) >= (imageRange(1,1)-0.5) & ...
        yCoord(:,1) <= (imageRange(1,2)+0.5) );
    
    movieInfoCropped(iFrame).xCoord = xCoord(indxKeep,:);
    movieInfoCropped(iFrame).yCoord = yCoord(indxKeep,:);
    movieInfoCropped(iFrame).amp    = amp(indxKeep,:);
    
end

%% Plot

%plot cropped coordinates
if doPlot
    xCoord = vertcat(movieInfoCropped.xCoord);
    yCoord = vertcat(movieInfoCropped.yCoord);
    imshow(currentImage);
    hold on
    plot(xCoord(:,1),yCoord(:,1),'g.')
end

%% ~~~ the end ~~~
