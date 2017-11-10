function trackDiffPaxInfo = sptDiffvsPaxInt(tracksFinal,diffModeAnRes,...
    numNeighbors,paxFrames,firstPaxImFile,firstPaxMaskFile,...
    firstCellMaskFile,saveDir,doPlot,lengthMinMax)
%SPTDIFFVSPAXINT correlates single molecule diffusion properties to local paxillin intensity
%
%SYNOPSIS trackDiffPaxInfo = sptDiffvsPaxInt(tracksFinal,diffModeAnRes,...
%    numNeighbors,paxFrames,firstPaxImFile,firstPaxMaskFile,...
%    firstCellMaskFile,saveDir,doPlot,lengthMinMax)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       diffModeAnRes  : Diffusion mode analysis results, output of
%                        trackDiffModeAnalysis.
%       numNeighbors   : Number of neighbors, output of numNeighborsTrack.
%       paxFrames      : The SPT frames at which there are paxillin images.
%       firstPaxImFile : Location and name of first paxillin image file.
%                        If not input, user will be prompted to choose.
%       firstPaxMaskFile: Location and name of first FA mask file.
%                        If not input, user will be prompted to choose.
%       firstCellMaskFile: Location and name of first cell mask file.
%                        If not input, user will be prompted to choose.
%       saveDir        : Directory where results are to be saved.
%                        If not input, user will be prompted to choose.
%       doPlot         : 1 to plot results, 0 otherwise. Plots will be
%                        saved in saveDir.
%                        Optional. Default: 0.
%       lengthMinMax   : Minimum and maximum length of trajectories to
%                        include in analysis.
%                        Optional. Default: [5 99]
%
%OUTPUT 
%
%REMARKS 
%
%Khuloud Jaqaman, January 2013

%% Input

if nargin < 4
    disp('--sptDiffVsPaxInt: Incorrect number of input arguments!');
    return
end

%ask user for paxillin images
if nargin < 5 || isempty(firstPaxImFile)
    [fName,dirName] = uigetfile('*.tif','specify first paxillin image in the stack - specify very first image');
else
    if iscell(firstPaxImFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxImFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstPaxImFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxImFile);
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
if nargin < 6 || isempty(firstPaxMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first FA mask in the stack - specify very first mask');
else
    if iscell(firstPaxMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxMaskFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstPaxMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstPaxMaskFile);
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

%ask user for cell masks
if nargin < 7 || isempty(firstCellMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first cell mask in the stack - specify very first mask');
else
    if iscell(firstCellMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstCellMaskFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstCellMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstCellMaskFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    
    %read images
    cellMaskStack = NaN(isx,isy,numPaxFrames);
    for iPaxFrame = 1 : numPaxFrames
        cellMaskStack(:,:,iPaxFrame) = imread(outFileList{iPaxFrame});
    end
    
else %else, exit
    disp('--sptDiffVsPaxInt: Bad file selection');
    return
end

%ask user for saving directory
if nargin < 8 || isempty(saveDir)
    saveDir = uigetdir([],'Choose directory to save results');
end

%plot flag
if nargin < 9 || isempty(doPlot)
    doPlot = 0;
end

%trajectory length
if nargin < 10 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end

%generate paxFrameMin and paxFramesMax for for loop
paxFramesMid = floor((paxFrames(1:end-1) + paxFrames(2:end))/2);
paxFramesMin = [1 paxFramesMid];
paxFramesMax = [paxFramesMid paxFrames(end)+1];

%% Track pre-processing

%keep only tracks of appropriate length
criteria.lifeTime.min = lengthMinMax(1);
criteria.lifeTime.max = lengthMinMax(2);
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx);
diffModeAnRes = diffModeAnRes(indx);
numNeighbors = numNeighbors(indx);

%get number of tracks
numTracksCompound = size(tracksFinal,1);

%get track segment start and end times
trackSEL = getTrackSEL(tracksFinal,1);

%get number of track segments
numTracks = size(trackSEL,1);

%calculate the "average" time at which a track exists
%group tracks based on what paxillin frame range they fall into
trackTimeMean = mean(trackSEL(:,1:2),2);
trackGroup = repmat(struct('indx',[]),numPaxFrames,1);
for iPaxFrame = 1 : numPaxFrames
   
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

%get track diffusion coefficient, mode and number of neighbors
trackMode = vertcat(diffModeAnRes.diffMode);
trackDiffCoef = vertcat(diffModeAnRes.diffCoef);
trackNeighbors = vertcat(numNeighbors.value);

%% Paxillin and FA readout

%filter paxillin images to dampen noise and subtract background
paxImageFilter = NaN(size(paxImageStack));
for iPaxFrame = 1 : numPaxFrames
    image0 = paxImageStack(:,:,iPaxFrame);
    imageF = filterGauss2D(image0,1);
    imageB = filterGauss2D(image0,10);
    paxImageFilter(:,:,iPaxFrame) = imageF - imageB;
end

%assemble a matrix of track information
%diffusion mode, diffusion coefficient, number of neighbors, associated
%paxillin intensity, inside or outside an FA
trackDiffPaxInfo = NaN(numTracks,5);
trackDiffPaxInfo(:,1:3) = [trackMode trackDiffCoef trackNeighbors];
for iPaxFrame = 1 : numPaxFrames
    
    %get tracks associated with this paxillin frame and their positions
    trackIndx = trackGroup(iPaxFrame).indx;
    trackPosX = xCoordMean(trackIndx);
    trackPosY = yCoordMean(trackIndx);
    
    %convert position to a linear index
    %transform from image to matrix coordinates
    linearInd = sub2ind([isx isy],trackPosY,trackPosX);
    
    %read associated paxillin intensity
    paxCurrent = paxImageFilter(:,:,iPaxFrame);
    paxInt = paxCurrent(linearInd);
    
    %read whether in FA or not
    maskCurrent = faMaskStack(:,:,iPaxFrame);
    inFaOrNot = maskCurrent(linearInd);
    
    %store in big vector of all tracks
    trackDiffPaxInfo(trackIndx,4:end) = [paxInt inFaOrNot];
    
end

%% Analysis of readout

%Analysis 1: Paxillin intensity distribution per diffusion mode
paxIntMode1 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==1,4);
paxIntMode2 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==2,4);
paxIntMode3 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==3,4);
paxIntMode4 = trackDiffPaxInfo(trackDiffPaxInfo(:,1)==4,4);

%get intensity range and calculate histogram bins
minInt = floor(min(trackDiffPaxInfo(:,4)));
maxInt = ceil(max(trackDiffPaxInfo(:,4)));
xPaxInt = linspace(minInt,maxInt,51);
xHalfInc = (xPaxInt(2) - xPaxInt(1))/2;
xPaxInt = xPaxInt(1:end-1) + xHalfInc;

%get paxillin intensity distributions
nPaxInt1 = hist(paxIntMode1,xPaxInt);
nPaxInt1 = nPaxInt1 / sum(nPaxInt1);
nPaxInt2 = hist(paxIntMode2,xPaxInt);
nPaxInt2 = nPaxInt2 / sum(nPaxInt2);
nPaxInt3 = hist(paxIntMode3,xPaxInt);
nPaxInt3 = nPaxInt3 / sum(nPaxInt3);
nPaxInt4 = hist(paxIntMode4,xPaxInt);
nPaxInt4 = nPaxInt4 / sum(nPaxInt4);

%Analysis 2: Diffusion modes inside and outside FAs
indxInFA = find(trackDiffPaxInfo(:,5)==1);
indxOutFA = find(trackDiffPaxInfo(:,5)==0);
trackModeInFA = trackDiffPaxInfo(indxInFA,1);
trackModeOutFA = trackDiffPaxInfo(indxOutFA,1);
xMode = 1 : 4;
nModeIn = hist(trackModeInFA,xMode);
nModeIn = nModeIn / sum(nModeIn);
nModeOut = hist(trackModeOutFA,xMode);
nModeOut = nModeOut / sum(nModeOut);

%Analysis 3: Paxillin intensity distribution per number of neighbors
indxNN0 = find( trackDiffPaxInfo(:,3)==0 );
fracNN0 = length(indxNN0) / numTracks;
paxIntNN0 = trackDiffPaxInfo(indxNN0,4);
indxNN1 = find( trackDiffPaxInfo(:,3)==1 );
fracNN1 = length(indxNN1) / numTracks;
paxIntNN1 = trackDiffPaxInfo(indxNN1,4);
indxNN23 = find( trackDiffPaxInfo(:,3)==2 | trackDiffPaxInfo(:,3)==3 );
fracNN23 = length(indxNN23) / numTracks;
paxIntNN23 = trackDiffPaxInfo(indxNN23,4);
indxNNMax = find(trackDiffPaxInfo(:,3)>3);
fracNNMax = length(indxNNMax) / numTracks;
paxIntNNMax = trackDiffPaxInfo(indxNNMax,4);

%get paxillin intensity distributions
nPaxIntNN0 = hist(paxIntNN0,xPaxInt);
nPaxIntNN0 = nPaxIntNN0 / sum(nPaxIntNN0);
nPaxIntNN1 = hist(paxIntNN1,xPaxInt);
nPaxIntNN1 = nPaxIntNN1 / sum(nPaxIntNN1);
nPaxIntNN23 = hist(paxIntNN23,xPaxInt);
nPaxIntNN23 = nPaxIntNN23 / sum(nPaxIntNN23);
nPaxIntNNMax = hist(paxIntNNMax,xPaxInt);
nPaxIntNNMax = nPaxIntNNMax / sum(nPaxIntNNMax);

%Analysis 4: Number of neighbors inside and outside FAs
trackNumNeighborsInFA = trackDiffPaxInfo(indxInFA,3);
trackNumNeighborsOutFA = trackDiffPaxInfo(indxOutFA,3);
neighborsMax = ceil(max(trackDiffPaxInfo(:,3)));
xNeighbors = 0 : neighborsMax;
nNeighborsIn = hist(trackNumNeighborsInFA,xNeighbors);
nNeighborsIn = nNeighborsIn / sum(nNeighborsIn);
nNeighborsOut = hist(trackNumNeighborsOutFA,xNeighbors);
nNeighborsOut = nNeighborsOut / sum(nNeighborsOut);

%Analysis 5: Particle density inside and outside FAs
numTracksIn = length(indxInFA);
numTracksOut = length(indxOutFA);
areaInFA = sum(faMaskStack(:));
areaTot = sum(cellMaskStack(:));
areaOutFA = areaTot - areaInFA;
densityInFA = numTracksIn / areaInFA;
densityOutFA = numTracksOut / areaOutFA;
densityTot = numTracks / areaTot;

%save results
save(fullfile(saveDir,'sptDiffVsPaxIntRes'),'trackDiffPaxInfo',...
    'paxIntMode1','paxIntMode2','paxIntMode3','paxIntMode4',...
    'xPaxInt','nPaxInt1','nPaxInt2','nPaxInt3','nPaxInt4',...
    'trackModeInFA','trackModeOutFA','xMode','nModeIn','nModeOut',...
    'paxIntNN0','paxIntNN1','paxIntNN23','paxIntNNMax',...
    'nPaxIntNN0','nPaxIntNN1','nPaxIntNN23','nPaxIntNNMax',...
    'trackNumNeighborsInFA','trackNumNeighborsOutFA',...
    'xNeighbors','nNeighborsIn','nNeighborsOut',...
    'densityTot','densityInFA','densityOutFA');

%% Plotting

if doPlot
    
    %Figure 1: Paxillin intensity distribution per diffusion mode
    h1 = figure('Name','Paxillin intensity distribution per diffusion mode');
    hold on
    plot(xPaxInt,nPaxInt1,'k')
    plot(xPaxInt,nPaxInt2,'b')
    plot(xPaxInt,nPaxInt3,'c')
    plot(xPaxInt,nPaxInt4,'m')
    legend('Mode 1','Mode 2','Mode 3','Mode 4')
    xlabel('Paxillin intensity (a.u)')
    ylabel('Normalized frequency')
    saveas(h1,fullfile(saveDir,'paxIntPerDiffMode'),'fig');
    
    %Figure 2: Diffusion modes inside and outside FAs
    h2 = figure('Name','Diffusion modes inside and outside FAs');
    hold on
    plot(xMode,nModeIn,'b','Marker','.')
    plot(xMode,nModeOut,'r','Marker','.')
    legend('Inside FAs','Outside FAs')
    xlabel('Diffusion mode')
    ylabel('Normalized frequency')
    saveas(h2,fullfile(saveDir,'diffModeInOutFAs'),'fig');
    
    %Figure 3: Paxillin intensity distribution per number of neighbors
    h3 = figure('Name','Paxillin intensity distribution per number of neighbors');
    hold on
    plot(xPaxInt,nPaxIntNN0,'k')
    plot(xPaxInt,nPaxIntNN1,'b')
    plot(xPaxInt,nPaxIntNN23,'c')
    plot(xPaxInt,nPaxIntNNMax,'m')
    legend(['NN 0 (frac tracks ' num2str(fracNN0) ')'],['NN 1(frac tracks ' num2str(fracNN1) ')'],['NN 2&3 (frac tracks ' num2str(fracNN23) ')'],['NN>=4 (frac tracks ' num2str(fracNNMax) ')'])
    xlabel('Paxillin intensity (a.u)')
    ylabel('Normalized frequency')
    saveas(h3,fullfile(saveDir,'paxIntPerNumNeighbors'),'fig');
    
    %Figure 4: Number of neighbors inside and outside FAs
    h4 = figure('Name','Number of neighbors inside and outside FAs');
    hold on
    plot(xNeighbors,nNeighborsIn,'b','Marker','.')
    plot(xNeighbors,nNeighborsOut,'r','Marker','.')
    legend('Inside FAs','Outside FAs')
    xlabel('Number of neighbors')
    ylabel('Normalized frequency')
    saveas(h4,fullfile(saveDir,'numNeighborsInOutFAs'),'fig');
    
    %Figure 5: Particle density inside and outside FAs
    h5 = figure('Name','Particle density inside and outside FAs');
    hold on
    bar([densityTot densityInFA densityOutFA])
    xlabel('Total Inside Outside')
    ylabel('Particle density (/pixel^2)');
    saveas(h5,fullfile(saveDir,'densityInOutFAs'),'fig');

end


