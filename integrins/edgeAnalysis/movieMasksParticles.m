function movieMasksParticles(movieInfo,numFramesSPT,firstMaskFile,...
    firstImageFile,saveMovie,movieName,dir2saveMovie,movieType,plotFullScreen)
%MOVIEMASKSPARTICLES makes a movie of cell masks and detected particles
%
%SYNPOSIS makeMovieMasksParticles(movieInfo,firstMaskFile,saveMovie,movieName,...
%    movieType,plotFullScreen)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       numFramesSPT: Number of particle frames between mask frames.
%       firstMaskFile: Name, including full path, of the first mask file.
%                     Optional. Default: [].
%       firstImageFile: Name, including full path, of the first image file.
%                     Optional. Default: [].
%       saveMovie   : 1 to save movie, 0 otherwise.
%                     Optional. Default: 0
%       movieName   : filename for saving movie.
%                     Optional. Default: masksParticles (if saveMovie = 1).
%       dir2saveMovie: Directory where to save output movie.
%                     If not input, movie will be saved in directory where
%                     masks are located.
%                     Optional. Default: [].
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

%% input - basic

%check whether correct number of input arguments was used
if nargin < 2
    disp('--movieMasksParticles: Incorrect number of input arguments!');
    return
end

%ask user for masks
if nargin < 3 || isempty(firstMaskFile)
    [fName,dirName] = uigetfile('*.tif','Specify first mask file');
else
    if iscell(firstMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstMaskFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstMaskFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFramesMovie = length(outFileList);
    
    %read first mask to get mask size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    
else %else, exit
    
    disp('--movieMasksParticles: Bad file selection');
    return
    
end

imageRange = [1 isx; 1 isy];

%ask user for images
if nargin < 4 || isempty(firstImageFile)
    [fName,dirName] = uigetfile('*.tif','Specify first image file');
else
    if iscell(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileListImages = getFileStackNames([dirName,fName]);
    
else %else, exit
    
    disp('--movieMasksParticles: Bad file selection');
    return
    
end

%% input - additional parameters

%check whether to save movie
if nargin < 5 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 6 || isempty(movieName))
    movieName = 'masksParticles';
end

%check where to save resulting movie
if saveMovie && (nargin < 7|| isempty(dir2saveMovie))
    dir2saveMovie = dirName;
end

%decide on movie type
if nargin < 8 || isempty(movieType)
    movieType = 'mov';
end

%check whether to use full screen for plotting
if nargin < 9 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

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
    imageStack = imread(outFileListImages{iFrame});
    imageStack1 = imread(outFileList{iFrame});
    imageStack2 = imread(outFileList{iFrame+1});
    
    %collect particle positions
    xCoordRange = vertcat(movieInfo((iFrame-1)*numFramesSPT+1:iFrame*numFramesSPT).xCoord);
    yCoordRange = vertcat(movieInfo((iFrame-1)*numFramesSPT+1:iFrame*numFramesSPT).yCoord);
    
    clf;
    
    %plot cell image + mask
    axes('Position',[0 0 0.495 1]);
    imshow(imageStack,[]);
    hold on;
    maskBounds = bwboundaries(imageStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    textDeltaCoord = min(diff(imageRange,[],2))/20;
    text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
        textDeltaCoord,num2str(iFrame),'Color','white');
    
    %plot particles + mask
    axes('Position',[0.505 0 0.495 1]);
    imshow(ones(isx,isy));
    hold on    
    
    %plot particles
    if ~isempty(xCoordRange)
        plot(xCoordRange(:,1),yCoordRange(:,1),'.');
    end
    
    %overlay masks
    maskBounds = bwboundaries(imageStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    maskBounds = bwboundaries(imageStack2);
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

%% ~~~ end ~~~

