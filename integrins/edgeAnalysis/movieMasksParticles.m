function movieMasksParticles(MD,movieInfo,numFramesSPT,saveMovie,movieName,...
    movieType,plotFullScreen)
%MOVIEMASKSPARTICLES makes a movie of cell masks and detected particles
%
%SYNPOSIS movieMasksParticles(MD,movieInfo,numFramesSPT,saveMovie,movieName,...
%    movieType,plotFullScreen)
%
%INPUT  MD          : movieData with cell mask.
%       movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       numFramesSPT: Number of particle frames between mask frames.
%       saveMovie   : 1 to save movie, 0 otherwise.
%                     Optional. Default: 0
%       movieName   : filename for saving movie.
%                     Optional. Default: masksParticles (if saveMovie = 1).
%       movieType   : 'mov' to make a Quicktime movie using MakeQTMovie,
%                     'avi' to make AVI movie using Matlab's movie2avi,
%                     'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                     using ImageMagick and ffmpeg. These options works
%                     only under linux or mac.
%                     Optional. Default: 'mov'.
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%
%OUTPUT the movie.
%
%Khuloud Jaqaman, September 2012

%% input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--movieMasksParticles: Incorrect number of input arguments!');
    return
end

%check whether to save movie
if nargin < 4 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 5 || isempty(movieName))
    movieName = 'masksParticles';
end

%decide on movie type
if nargin < 6 || isempty(movieType)
    movieType = 'mov';
end

%check whether to use full screen for plotting
if nargin < 7 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

%% preamble

%get images
imageDir = MD.channels_.channelPath_;
imageFileListing = dir(imageDir);
imageFileListing = keepOnlyTiffFiles(imageFileListing);

%get masks
maskDir = [MD.movieDataPath_ filesep 'refined_masks' filesep 'refined_masks_for_channel_1'];
maskFileListing = dir(maskDir);
maskFileListing = keepOnlyTiffFiles(maskFileListing);

%get number of frames
numFramesMovie = length(imageFileListing);

%read first image to get image size
currentImage = imread(fullfile(imageDir,imageFileListing(1).name));
[isx,isy] = size(currentImage);
imageRange = [1 isx; 1 isy];

%determine where to save movie
dir2saveMovie = MD.movieDataPath_;

%% make movie

%initialize movie if it is to be saved
if saveMovie
    movieVar = struct('cdata',[],'colormap',[]);
    movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%go over all frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h = figure();
    set(h,'Position',scrsz);
else
    figure
end
for iFrame = 1 : numFramesMovie - 1
    
    %read image + masks
    imageStack = imread(fullfile(imageDir,imageFileListing(iFrame).name));
    maskStack1 = imread(fullfile(maskDir,maskFileListing(iFrame).name));
    maskStack2 = imread(fullfile(maskDir,maskFileListing(iFrame+1).name));
    
    %collect particle positions
    xCoordRange = vertcat(movieInfo((iFrame-1)*numFramesSPT+1:iFrame*numFramesSPT).xCoord);
    yCoordRange = vertcat(movieInfo((iFrame-1)*numFramesSPT+1:iFrame*numFramesSPT).yCoord);
    
    clf;
    
    %plot cell image + mask
    axes('Position',[0 0 0.495 1]);
    imshow(imageStack,[]);
    hold on;
    maskBounds = bwboundaries(maskStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    textDeltaCoord = min(diff(imageRange,[],2))/20;
    text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
        textDeltaCoord,num2str(iFrame),'Color','white');
    %     text(imageRange(1,1)+textDeltaCoord-1,imageRange(2,1)+...
    %         textDeltaCoord+2,[num2str((iFrame-1)*10) ' s'],'Color','white','FontSize',30);
    
    %plot particles + mask
    axes('Position',[0.505 0 0.495 1]);
    imshow(ones(isx,isy));
    hold on    
    
    %plot particles
    if ~isempty(xCoordRange)
        plot(xCoordRange(:,1),yCoordRange(:,1),'.');
    end
    
    %overlay masks
    maskBounds = bwboundaries(maskStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    maskBounds = bwboundaries(maskStack2);
    cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',2)),maskBounds);
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFramesMovie,movieVar,iFrame);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
end

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%% Sub-function

function fileListing = keepOnlyTiffFiles(fileListing)

numFiles = length(fileListing);
goodFile = zeros(numFiles,1);
for iFile = 1 : numFiles
    fileName = fileListing(iFile).name;
    if (length(fileName) > 4) && ((strcmp(fileName(end-2:end),'tif')) || ...
            (strcmp(fileName(end-3:end),'tiff')))
        goodFile(iFile) = 1;
    end
end
fileListing = fileListing(goodFile==1);

%% ~~~ end ~~~

