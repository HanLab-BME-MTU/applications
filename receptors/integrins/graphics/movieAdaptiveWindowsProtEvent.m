function movieAdaptiveWindowsProtEvent(MD,eventWindows,winPositions,...
    frameOnset,frameStart,frameEnd,saveMovie,movieName,movieType,...
    plotFullScreen,axisLim)

%% input

%check whether correct number of input arguments was used
if nargin < 6
    disp('--movieMasksParticles: Incorrect number of input arguments!');
    return
end

%check whether to save movie
if nargin < 7 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 8 || isempty(movieName))
    movieName = 'adaptiveWindows';
end

%decide on movie type
if nargin < 9 || isempty(movieType)
    movieType = 'mov';
end

%check whether to use full screen for plotting
if nargin < 10 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

%check image range to plot
if nargin < 11 || isempty(axisLim)
    axisLim = [];
end

%% Collect windows

%get number of slices
numSlice = length(eventWindows);

%static windows
winsOnsetFrame = winPositions(frameOnset,:);
for iSlice = 1 : numSlice
    indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
    winsOnset{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,1}(:,1));
    winsBand2{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,2}(:,1));
    winsBand3{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,3}(:,1));
    tmp = [];
    for iBigBand = 4 : 10
        if ~isempty(eventWindows(iSlice).onset{1,iBigBand})
            tmp = [tmp; eventWindows(iSlice).onset{1,iBigBand}(:,1)];
        end
    end
    winsLamella{iSlice} = winsOnsetFrame{indxSlice}(tmp);
end

%dynamic windows before
for iFrame = frameStart : frameOnset-1
    winsCurrentFrame = winPositions(iFrame,:);
    frameDiff = frameOnset - iFrame;
    for iSlice = 1 : numSlice
        indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
        winsBefDynamic{iSlice,frameDiff} = winsCurrentFrame{indxSlice}(eventWindows(iSlice).befDynamic{frameDiff,1}(:,1));
    end
end

%dynamic windows after
for iFrame = frameOnset + 1 : frameEnd
    winsCurrentFrame = winPositions(iFrame,:);
    frameDiff = iFrame - frameOnset;
    for iSlice = 1 : numSlice
        indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
        winsAftDynamic{iSlice,frameDiff} = winsCurrentFrame{indxSlice}(eventWindows(iSlice).aftDynamic{frameDiff,frameDiff}(:,1));
    end    
end
 
%% Prepare for movie

%get images
imageDir = MD.channels_.channelPath_;
imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end

%get masks
maskDir = MD.processes_{2}.outFilePaths_{1};
maskFileListing = dir([maskDir filesep '*.tif']);
if isempty(maskFileListing)
    maskFileListing = dir([maskDir filesep '*.tiff']);
end

%get number of frames
numFramesMovie = length(imageFileListing);

%read first image to get image size
currentImage = imread(fullfile(imageDir,imageFileListing(1).name));
[isx,isy] = size(currentImage);
imageRange = [1 isx; 1 isy];

%determine where to save movie
dir2saveMovie = MD.movieDataPath_;

%% Make movie

%initialize movie if it is to be saved
if saveMovie
    movieVar = struct('cdata',[],'colormap',[]);
    movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%initialize figure
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h = figure();
    set(h,'Position',scrsz);
else
    figure
end

%before protrusion onset
for iFrame = frameStart : frameEnd
    
    clf;
    
    %get delta time from onset
    frameDiff = iFrame - frameOnset;
    
    %read image + mask
    imageStack = imread(fullfile(imageDir,imageFileListing(iFrame).name));
    maskStack1 = imread(fullfile(maskDir,maskFileListing(iFrame).name));
    
    %plot cell image + mask
    imshow(imageStack,[prctile(imageStack(:),5) prctile(imageStack(:),95)]);
    hold on;
    maskBounds = bwboundaries(imdilate(maskStack1,strel('disk',1)));
    cellfun(@(x)(plot(x(:,2),x(:,1),'Color',[1 0.7 0],'LineWidth',1)),maskBounds);
    axisLim2 = [axisLim(1:2); axisLim(3:4)];
    textDeltaCoord = min(diff(axisLim2,[],2))/20;
    text(axisLim2(1,1)+textDeltaCoord,axisLim2(2,1)+...
        textDeltaCoord,[num2str(frameDiff) ' s'],'Color','white');
    
    if frameDiff < 0 %before protrusion onset
        
        %switch sign
        frameDiff = -frameDiff;
        
        %right behind edge
        plotWindows([winsBefDynamic{:,frameDiff}],{'c','FaceAlpha',0.5,'EdgeColor','c'});
        
    elseif frameDiff == 0 % at onset
        
        %right behind edge
        plotWindows(winsOnset,{'c','FaceAlpha',0.5,'EdgeColor','c'});
        
    else %during protrusion
        
        %new protrusion area
        plotWindows([winsAftDynamic{:,1:frameDiff}],{'r','FaceAlpha',0.5,'EdgeColor','r'});
        
        %right behind edge
        plotWindows([winsAftDynamic{:,frameDiff}],{'c','FaceAlpha',0.5,'EdgeColor','c'});
        
    end
    
    %static windows ...
    %onset
    plotWindows(winsOnset,{'k','FaceAlpha',0.5,'EdgeColor','k'});
    %band 2
    plotWindows(winsBand2,{'b','FaceAlpha',0.5,'EdgeColor','b'});
    %band 3
    plotWindows(winsBand3,{'g','FaceAlpha',0.5,'EdgeColor','g'});
    %lamella
    plotWindows(winsLamella,{'m','FaceAlpha',0.5,'EdgeColor','m'});
    
    %zoom in on area of interest
    axis(axisLim);
    
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

%% ~~~ the end ~~~
