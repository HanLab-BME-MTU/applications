function trackDiffPaxInfo = sptDiffvsPaxInt(tracksFinal,diffModeAnRes,...
    paxFrames,firstPaxFile,firstMaskFile,saveDir,doPlot)
%SPTDIFFVSPAXINT correlates single molecule diffusion properties to local paxillin intensity
%
%SYNOPSIS trackDiffPaxInfo = sptDiffvsPaxInt(tracksFinal,diffModeAnRes,...
%    paxFrames,firstPaxFile,firstMaskFile,saveDir,doPlot)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       diffModeAnRes  : Diffusion mode analysis results, output of
%                        trackDiffModeAnalysis.
%       paxFrames      : The frames at which there are paxillin images.
%       firstPaxFile   : Location and name of first paxillin file.
%                        If not input, user will be prompted to choose.
%       firstMaskFile  : Location and name of first FA mask file.
%                        If not input, user will be prompted to choose.
%       saveDir        : Directory where results are to be saved.
%                        If not input, user will be prompted to choose.
%       doPlot         : 1 to plot results, 0 otherwise. Plots will be
%                        saved in saveDir.
%                        Optional. Default: 0.
%
%OUTPUT 
%
%REMARKS 
%
%Khuloud Jaqaman, January 2013

%% Input

if nargin < 3
    disp('--sptDiffVsPaxInt: Incorrect number of input arguments!');
    return
end

%ask user for paxillin images
if nargin < 4 || isempty(firstPaxFile)
    [fName,dirName] = uigetfile('*.tif','specify first paxillin image in the stack - specify very first image');
else
    if iscell(firstPaxFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstPaxFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numPaxFrames = length(outFileList);
    
    %read images
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    paxImageStack = NaN(isx,isy,numPaxFrames);
    paxImageStack(:,:,1) = currentImage;
    for iPaxFrame = 2 : numPaxFrames
        paxImageStack(:,:,iPaxFrame) = imread(outFileList{iPaxFrame});
    end
    
else %else, exit
    disp('--sptDiffVsPaxInt: Bad file selection');
    return
end

%ask user for FA masks
if nargin < 5 || isempty(firstMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first FA mask in the stack - specify very first mask');
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
    faMaskStack = NaN(isx,isy,numPaxFrames);
    for iPaxFrame = 1 : numPaxFrames
        faMaskStack(:,:,iPaxFrame) = imread(outFileList{iPaxFrame});
    end
    
else %else, exit
    disp('--sptDiffVsPaxInt: Bad file selection');
    return
end

%ask user for saving directory
if nargin < 6 || isempty(saveDir)
    saveDir = uigetdir([],'Choose directory to save results');
end

%plot flag
if nargin < 7 || isempty(doPlot)
    doPlot = 0;
end

%generate paxFrameMin and paxFramesMax for parfor loop
paxFramesMin = paxFrames(1:end-1);
paxFramesMax = paxFrames(2:end);

%% Track pre-processing

%get number of tracks
numTracksCompound = size(tracksFinal,1);

%get track segment start and end times
trackSEL = getTrackSEL(tracksFinal,1);

%get number of track segments
numTracks = size(trackSEL,1);

%calculate the "average" time at which a track exists
%group tracks based on what paxillin frame range they fall into
trackTimeMean = mean(trackSEL(:,1:2),2);
trackGroup = repmat(struct('indx',[]),numPaxFrames-1,1);
for iPaxFrame = 1 : numPaxFrames - 1
   
    %get current spt frame number and next spt frame number
    minFrame = paxFramesMin(iPaxFrame);
    maxFrame = paxFramesMax(iPaxFrame);
        
    %find tracks whose "average" time is between minFrame and maxFrame
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);

    %store this information for later use
    trackGroup(iPaxFrame).indx = indxFrameRange;
    
end

%get average track positions
if isstruct(tracksFinal) %if compound tracks
    
    %initilize mean position
    [xCoordMean,yCoordMean] = deal(NaN(numTracks,1));
    
    %get mean position based on track segments
    iSeg = 0;
    for iTrack = 1 : numTracksCompound
        xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
        yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
        numSeg = size(xCoordAll,1);
        xCoordMean(iSeg+1:iSeg+numSeg) = nanmean(xCoordAll,2);
        yCoordMean(iSeg+1:iSeg+numSeg) = nanmean(yCoordAll,2);
        iSeg = iSeg + numSeg;
    end
    
else %if simple tracks
    
    %extract x- and y-coordinates
    xCoordAll = tracksFinal(:,1:8:end);
    yCoordAll = tracksFinal(:,2:8:end);
    
    %calculate the average coordinates
    xCoordMean = nanmean(xCoordAll,2);
    yCoordMean = nanmean(yCoordAll,2);
    
end

%round average track positions to the closest pixel
xCoordMean = round(xCoordMean);
yCoordMean = round(yCoordMean);

%get track diffusion coefficient and mode
trackMode = vertcat(diffModeAnRes.diffMode);
trackDiffCoef = vertcat(diffModeAnRes.diffCoef);

%% Paxillin and FA readout

%filter paxillin images to dampen noise and subtract background
paxImageFilter = NaN(size(paxImageStack));
for iPaxFrame = 1 : numPaxFrames
    image0 = paxImageStack(:,:,iPaxFrame);
    imageF = filterGauss2D(image0,1);
    imageB = filterGauss2D(image0,10);
    paxImageFilter(:,:,iPaxFrame) = imageF - imageB;
end

%go over track groups and read out filtered paxillin intensity per track
%also read whether each track is in an FA or not
trackDiffPaxInfo = NaN(numTracks,6);
trackDiffPaxInfo(:,1:2) = [trackMode trackDiffCoef];
for iPaxFrame = 1 : numPaxFrames - 1
    
    %get tracks in this time window and their positions
    trackIndx = trackGroup(iPaxFrame).indx;
    trackPosX = xCoordMean(trackIndx);
    trackPosY = yCoordMean(trackIndx);
    
    %convert position to a linear index
    %transform from image to matrix coordinates
    linearInd = sub2ind([isx isy],trackPosY,trackPosX);
    
    %read paxillin intensity in frame before and frame after
    paxBef = paxImageFilter(:,:,iPaxFrame);
    paxAft = paxImageFilter(:,:,iPaxFrame+1);
    paxIntBef = paxBef(linearInd);
    paxIntAft = paxAft(linearInd);
    
    %read whether in FA or not in frame before and frame after
    maskBef = faMaskStack(:,:,iPaxFrame);
    maskAft = faMaskStack(:,:,iPaxFrame+1);
    inFaBef = maskBef(linearInd);
    inFaAft = maskAft(linearInd);
    
    %store in big vector of all tracks
    trackDiffPaxInfo(trackIndx,3:end) = [paxIntBef paxIntAft inFaBef inFaAft];
    
end
    
%% Plotting and saving results

if doPlot
    
    %get paxillin intensity at locations of tracks in the different modes
    paxIntMode1 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==1,3);
    paxIntMode2 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==2,3);
    paxIntMode3 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==3,3);
    paxIntMode4 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==4,3);
    
    %get intensity range and calculate histogram bins
    minInt = floor(min(trackDiffPaxInfo(:,3)));
    maxInt = ceil(max(trackDiffPaxInfo(:,3)));
    xPaxInt = linspace(minInt,maxInt,50);
    
    %get paxillin intensity corresponding to these tracks
    nPaxInt1 = hist(paxIntMode1,xPaxInt);
    nPaxInt2 = hist(paxIntMode2,xPaxInt);
    nPaxInt3 = hist(paxIntMode3,xPaxInt);
    nPaxInt4 = hist(paxIntMode4,xPaxInt);
    
    %plot figure and save
    h1 = figure('Name','Paxillin intensity distribution per diffusion mode');
    hold on
    plot(xPaxInt,nPaxInt1/sum(nPaxInt1),'k')
    plot(xPaxInt,nPaxInt2/sum(nPaxInt2),'b')
    plot(xPaxInt,nPaxInt3/sum(nPaxInt3),'c')
    plot(xPaxInt,nPaxInt4/sum(nPaxInt4),'m')
    legend('Mode 1','Mode 2','Mode 3','Mode 4')
    xlabel('Paxillin intensity (a.u)')
    ylabel('Normalized frequency')
    saveas(h1,fullfile(saveDir,'paxIntPerDiffMode'),'fig');
    
    %get track modes inside or outside focal adhesions
    trackModeInFA = trackDiffPaxInfo(trackDiffPaxInfo(:,5)==1|trackDiffPaxInfo(:,6)==1,1);
    trackModeOutFA = trackDiffPaxInfo(trackDiffPaxInfo(:,5)==0&trackDiffPaxInfo(:,6)==0,1);
    nModeIn = hist(trackModeInFA,1:4);
    nModeOut = hist(trackModeOutFA,1:4);
    
    %plot figure and save
    h2 = figure('Name','Diffusion modes inside and outside FAs');
    hold on
    plot(1:4,nModeIn/sum(nModeIn),'b')
    plot(1:4,nModeOut/sum(nModeOut),'r')
    legend('Inside FAs','Outside FAs')
    xlabel('Diffusion mode')
    ylabel('Normalized frequency')
    saveas(h2,fullfile(saveDir,'diffModeInOutFAs'),'fig');
    
    %save results
    save(fullfile(saveDir,'sptDiffVsPaxIntRes'),'trackDiffPaxInfo',...
        'paxIntMode1','paxIntMode2','paxIntMode3','paxIntMode4',...
        'xPaxInt','nPaxInt1','nPaxInt2','nPaxInt3','nPaxInt4',...
        'trackModeInFA','trackModeOutFA','nModeIn','nModeOut');
    
else
    
    %save results
    save(fullfile(saveDir,'sptDiffVsPaxIntRes'),'trackDiffPaxInfo');
    
end