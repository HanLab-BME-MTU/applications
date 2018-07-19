function movieAlternateCellAndSingleMol(firstCellFile,firstSingleMolFile,...
    numSMFilesBetween,saveMovie,movieName,dir2saveMovie,movieType,intensityScale,...
    plotFullScreen)
%MOVIEALTERNATECELLANDSINGLEMOL makes a movie with alternating cell and single molecule frames
%
%Khuloud Jaqaman, October 2012

%% input

%ask user for cell images
if nargin < 1 || isempty(firstCellFile)
    [fName,dirName] = uigetfile('*.tif','specify first image in the CELL stack - specify very first image, even if not to be plotted');
else
    if iscell(firstCellFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstCellFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstCellFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstCellFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileListCell = getFileStackNames([dirName,fName]);
    numFilesCell = length(outFileListCell);
    
    %determine which frames the files correspond to, and generate the inverse map
    %indicate missing frames with a zero
    frame2fileMapCell = zeros(numFilesCell,1);
    for iFile = 1 : numFilesCell
        [~,~,frameNumStr] = getFilenameBody(outFileListCell{iFile});
        frameNum = str2double(frameNumStr);
        frame2fileMapCell(frameNum) = iFile;
    end
    
    %assign as number of frames the last frame number observed
    numFramesCell = frameNum;
    
    %read first image to get image size
    currentImage = imread(outFileListCell{1});
    [isx,isy] = size(currentImage);
    
else %else, exit
    
    disp('--overlayFeaturesMovie: Bad file selection');
    return
    
end

%ask user for single molecule images
if nargin < 2 || isempty(firstSingleMolFile)
    [fName,dirName] = uigetfile('*.tif','specify first image in the CELL stack - specify very first image, even if not to be plotted');
else
    if iscell(firstSingleMolFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstSingleMolFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstSingleMolFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstSingleMolFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileListSM = getFileStackNames([dirName,fName]);
    numFilesSM = length(outFileListSM);
    
    %determine which frames the files correspond to, and generate the inverse map
    %indicate missing frames with a zero
    frame2fileMapSM = zeros(numFilesSM,1);
    for iFile = 1 : numFilesSM
        [~,~,frameNumStr] = getFilenameBody(outFileListSM{iFile});
        frameNum = str2double(frameNumStr);
        frame2fileMapSM(frameNum) = iFile;
    end
    
    %assign as number of frames the last frame number observed
    numFramesSM = frameNum;
    
else %else, exit
    
    disp('--overlayFeaturesMovie: Bad file selection');
    return
    
end

%get number of single molecule files between cell files
if nargin < 3 || isempty(numSMFilesBetween)
    numSMFilesBetween = 400;
end

%check startend and assign default if necessary
% if nargin < 4 || isempty(startend)
startend = [1 numFramesCell];
% else
%     startend(2) = min(startend(2),numFramesCell); %make sure that last frame does not exceed real last frame
% end

%keep only the frames of interest
outFileListCell = outFileListCell(frame2fileMapCell(startend(1)):frame2fileMapCell(startend(2)));
frame2fileMapCell = frame2fileMapCell(startend(1):startend(2));
indxNotZero = find(frame2fileMapCell~=0);
frame2fileMapCell(indxNotZero) = frame2fileMapCell(indxNotZero) - frame2fileMapCell(indxNotZero(1)) + 1;

%get number of frames in movie to be made
numFramesMovie = diff(startend) + 1;

%get image size
imageRange = [1 isx; 1 isy];

%check whether to save movie
if nargin < 4 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 5 || isempty(movieName))
    movieName = 'alternatingMovie';
end

%check where to save resulting movie
if saveMovie && (nargin < 6 || isempty(dir2saveMovie))
    dir2saveMovie = dirName;
end

%decide on movie type
if nargin < 7 || isempty(movieType)
    movieType = 'mov';
end

%check how to scale image intensity
if nargin < 8 || isempty(intensityScale)
    intensityScale = 1;
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
        movieName,2*numFramesMovie-1,movieVar,[]);
end

%go over all specified frames and find minimum and maximum intensity in all
%of them combined
% switch intensityScale
%     case 0
%         intensityMinMax = [];
%     case 1
%         meanIntensity = zeros(numFramesMovie,1);
%         stdIntensity = meanIntensity;
%         for iFrame = 1 : numFramesMovie
%             if frame2fileMap(iFrame) ~= 0
%                 imageStack = double(imread(outFileList{frame2fileMap(iFrame)}));
%                 meanIntensity(iFrame) = mean(imageStack(:));
%                 stdIntensity(iFrame) = std(imageStack(:));
%             end
%         end
%         meanIntensity = mean(meanIntensity);
%         stdIntensity = mean(stdIntensity);
%         intensityMinMax = [meanIntensity-2*stdIntensity meanIntensity+6*stdIntensity];
%     case 2
%         minIntensity = zeros(numFramesMovie,1);
%         maxIntensity = minIntensity;
%         for iFrame = 1 : numFramesMovie
%             if frame2fileMap(iFrame) ~= 0
%                 imageStack = double(imread(outFileList{frame2fileMap(iFrame)}));
%                 minIntensity(iFrame) = min(imageStack(:));
%                 maxIntensity(iFrame) = max(imageStack(:));
%             end
%         end
%         minIntensity = min(minIntensity);
%         maxIntensity = max(maxIntensity);
%         intensityMinMax = [minIntensity maxIntensity];
% end
intensityMinMax = [];

%go over all specified frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h     = figure();
    set(h,'Position',scrsz);
else
    figure
end
for iFrame = 1 : numFramesMovie
    
    iCellFrame = frame2fileMapCell(iFrame);
    iSMFrameStart = (iCellFrame-1)*numSMFilesBetween+2;
    iSMFrameEnd = iCellFrame*numSMFilesBetween;
    
    if iCellFrame ~= 0 %if frame exists
        
        %read specified images
        imageStackCell = imread(outFileListCell{iCellFrame});
        if iFrame ~= numFramesMovie
            imageStackSM = imread(outFileListSM{iCellFrame});
        end
        
    else %otherwise
        
        %make empty frame
        imageStackCell = zeros(isx,isy);
        imageStackSM = zeros(isx,isy);
        
    end
    
    %plot cell image
    clf;
    
    axes('Position',[0 0 1 1]);
    imshow(imageStackCell,intensityMinMax);
    xlim(imageRange(2,:));
    ylim(imageRange(1,:));
    hold on;
    textDeltaCoord = min(diff(imageRange,[],2))/20;
    text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
        textDeltaCoord,['Frame ' num2str(iSMFrameStart-1) ' / Elapsed time: ' num2str((iFrame+startend(1)-2)*10) 's'],'Color','white');
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,2*numFramesMovie-1,movieVar,iFrame*2-1);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
    if iFrame ~= numFramesMovie
        
        %plot single molecule image
        clf;
        
        axes('Position',[0 0 1 1]);
        imshow(imageStackSM,intensityMinMax);
        xlim(imageRange(2,:));
        ylim(imageRange(1,:));
        hold on;
        textDeltaCoord = min(diff(imageRange,[],2))/20;
        text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
            textDeltaCoord,['Frames ' num2str(iSMFrameStart) '-' num2str(iSMFrameEnd)],'Color','white');
        
        %add frame to movie if movie is saved
        if saveMovie
            movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
                movieName,2*numFramesMovie-1,movieVar,iFrame*2);
        end
        
        %pause for a moment to see frame
        pause(0.1);
        
    end
    
end

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,2*numFramesMovie-1,movieVar,[]);
end

%% ~~~ end ~~~

