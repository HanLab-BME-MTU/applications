function movieFAMasks(firstImageFile,firstMaskFile,saveMovie,movieName,...
    movieType,dir2saveMovie,plotFullScreen,filterNoise,filterBackground)
%movieFAMasks makes a movie of focal adhesion masks and displays filtered image that was segmented
%
%SYNPOSIS 
%
%INPUT  firstImageFile: Name and location of first image file in stack.
%       firstMaskFile : Name and location of first mask file in stack.
%       saveMovie     : 1 to save movie, 0 otherwise.
%                       Optional. Default: 0
%       movieName     : filename for saving movie.
%                       Optional. Default: masksParticles (if saveMovie = 1).
%       movieType     : 'mov' to make a Quicktime movie using MakeQTMovie,
%                       'avi' to make AVI movie using Matlab's movie2avi,
%                       'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                       using ImageMagick and ffmpeg. These options works
%                       only under linux or mac.
%                       Optional. Default: 'mov'.
%       dir2saveMovie : Directory where to save output movie.
%                       If not input, movie will be saved in directory where
%                       masks are located.
%                       Optional. Default: [].
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%       filterNoise   : Either 0 to not filter noise or filter sigma > 0 to
%                       filter noise.
%                       Optional. Default: 1.
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 10.
%                       The idea is to use here the same filtering
%                       parameters applied to the image when it was
%                       segmented, in order to display the filtered image.
%       
%OUTPUT the movie.
%
%Khuloud Jaqaman, January 2013

%% input

%first image file
if nargin < 1 || isempty(firstImageFile)
    [fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image');
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
    outFileList = getFileStackNames([dirName,fName]);
    numFrames = length(outFileList);
    
    %read images
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    imageStack = NaN(isx,isy,numFrames);
    imageStack(:,:,1) = currentImage;
    for iFrame = 2 : numFrames
        imageStack(:,:,iFrame) = imread(outFileList{iFrame});
    end
    
else %else, exit
    disp('--movieFAMasks: Bad file selection');
    return
end

%first mask file
if nargin < 2 || isempty(firstMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first mask in the stack - specify very first mask');
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
    
    %read images
    maskStack = NaN(isx,isy,numFrames);
    for iFrame = 1 : numFrames
        maskStack(:,:,iFrame) = imread(outFileList{iFrame});
    end
    
else %else, exit
    disp('--movieFAMasks: Bad file selection');
    return
end

%check whether to save movie
if nargin < 3 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 4 || isempty(movieName))
    movieName = 'masksFAs';
end

%decide on movie type
if nargin < 5 || isempty(movieType)
    movieType = 'mov';
end

%check where to save resulting movie
if saveMovie && (nargin < 6 || isempty(dir2saveMovie))
    dir2saveMovie = dirName;
end

%check whether to use full screen for plotting
if nargin < 7 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

%noise filtering
if nargin < 8 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 9 || isempty(filterBackground)
    filterBackground = 10;
end

%% make movie

%initialize movie if it is to be saved
if saveMovie
    movieVar = struct('cdata',[],'colormap',[]);
    movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
        movieName,numFrames,movieVar,[]);
end

%go over all frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h = figure();
    set(h,'Position',scrsz);
else
    figure
end
for iFrame = 1 : numFrames
    
    clf;
    
    %read image + mask
    image = imageStack(:,:,iFrame);
    mask = maskStack(:,:,iFrame);
    
    %remove noise by filtering image with a narrow Gaussian
    if filterNoise > 0
        imageFiltered = filterGauss2D(image,filterNoise);
    else
        imageFiltered = image;
    end
    
    %estimate background by filtering image with a wide Gaussian
    if filterBackground > 0
        imageBackground = filterGauss2D(image,filterBackground);
    else
        imageBackground = zeros(size(image));
    end
    
    %calculate noise-filtered and background-subtracted image
    imageFilteredMinusBackground = imageFiltered - imageBackground;
    
    %plot filtered image
    axes('Position',[0 0 0.495 1]);
    imshow(imageFilteredMinusBackground,[]);
    %     hold on
    %     text(10,20,[num2str((iFrame-1)*10) ' s'],'Color','white');
    %     plot([80 207],[50 50],'y:','LineWidth',0.5)
    %     plot([80 207],[177 177],'y:','LineWidth',0.5)
    %     plot([80 80],[50 177],'y:','LineWidth',0.5)
    %     plot([207 207],[50 177],'y:','LineWidth',0.5)
    
    %plot mask boundaries on top of original image
    axes('Position',[0.505 0 0.495 1]);
    imshow(image,[]);
    hold on;
    maskBounds = bwboundaries(mask);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',1)),maskBounds);
    %     plot([80 207],[50 50],'y:','LineWidth',0.5)
    %     plot([80 207],[177 177],'y:','LineWidth',0.5)
    %     plot([80 80],[50 177],'y:','LineWidth',0.5)
    %     plot([207 207],[50 177],'y:','LineWidth',0.5)
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFrames,movieVar,iFrame);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
end

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,numFrames,movieVar,[]);
end

%% ~~~ end ~~~

