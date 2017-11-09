function moviePropertySpatialMap2D(MD,numFramesSPT,tracksFinal,diffAnalysisRes,...
    diffModeAnRes,numNeighbors,property2plot,position2plot,lengthMinMax,...
    fixedSegLength,saveMovie,movieName,movieType,plotFullScreen)
%MOVIEPROPERTYSPATIALMAP2D creates movie with spatial maps of trajectory properties
%
%SYNOPSIS moviePropertySpatialMap2D(MD,numFramesSPT,tracksFinal,diffAnalysisRes,...
%    diffModeAnRes,numNeighbors,property2plot,position2plot,lengthMinMax,...
%    fixedSegLength,saveMovie,movieName,movieType,plotFullScreen)
%
%INPUT  MD             : movieData with cell mask.
%       numFramesSPT   : Number of particle frames between mask frames.
%       tracksFinal    : Output of trackCloseGapsKalman.
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%       diffModeAnRes  : Output of trackDiffModeAnalysis.
%       numNeighbors   : Output of numNeighborsTrack.
%       property2plot  : Property to plot:
%                        1 - MSS Diffusion classification.
%                        2 - MSS diffusion coefficient.
%                        3 - MSS Confinement radius.
%                        4 - Track lifetime.
%                        5 - Average frame-to-frame displacement.
%                        6 - Diffusion mode.
%                        7 - Diffusion coefficient from mode analysis.
%                        8 - Number of neighbors.
%                        Optional. Default: diffusion mode (6).
%	    position2plot  : Trajectory position to plot:
%                        1 - Center position.
%                        2 - Start position.
%                        3 - End position.
%                        Optional. Default: Center position (1).
%       lengthMinMax   : Minimum and maximum length of trajectories to
%                        include in plots.
%                        Optional. Default: [5 99].
%       fixedSegLength : 1 to divide a property range into segments of
%                        fixed length, 0 to divide a property range into
%                        segments of variable length such that they contain
%                        an equal number of elements.
%                        Optional. Default: 1;
%       saveMovie      : 1 to save movie, 0 otherwise.
%                        Optional. Default: 0
%       movieName      : filename for saving movie.
%                        Optional. Default: movieSpatialMap (if saveMovie = 1).
%       movieType      : 'mov' to make a Quicktime movie using MakeQTMovie,
%                        'avi' to make AVI movie using Matlab's movie2avi,
%                        'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                        using ImageMagick and ffmpeg. These options works
%                        only under linux or mac.
%                        Optional. Default: 'mov'.
%       plotFullScreen : 1 the figure will be sized to cover the whole
%                        screen. In this way the movie will be of highest
%                        possible quality. default is 0.
%
%OUTPUT the movie.
%
%Khuloud Jaqaman, September 2014

%% Input

if nargin < 6
    disp('--moviePropertySpatialMap2D: Incorrect number of input arguments!');
    return
end

if nargin < 7 || isempty(property2plot)
    property2plot = 6;
end

if nargin < 8 || isempty(position2plot)
    position2plot = 1;
end

if nargin < 9 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end

if nargin < 10 || isempty(fixedSegLength)
    fixedSegLength = 1;
end

%check whether to save movie
if nargin < 11 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 12 || isempty(movieName))
    movieName = 'movieSpatialMap';
end

%decide on movie type
if nargin < 13 || isempty(movieType)
    movieType = 'mov';
end

%check whether to use full screen for plotting
if nargin < 14 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

%% Preparation for plotting

%define the default of 10 segments for plotting some of the properties
numSegments = 10;

%assign segment colors
segmentColor = [0 0 0; 0 0 1; 0.2 0.7 0.7; 0 1 1; 0 1 0; ...
    0.6824 0.4667 0; 1 0.7 0; 1 0 0; 1 0 1; 0.7 0 1];

%construct parts of figure titles
plottedProperty = {'MSS classification','MSS diffusion coefficient',...
    'Confinement radius','Lifetime','Frame-to-frame displacement',...
    'Diffusion mode','Mode analysis diffusion coefficient','Number of neighbors'};
plottedPosition = {'center position','start position','end position'};

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

%% Trajectory pre-processing

%keep only trajectories of acceptable length
criteria.lifeTime.min = lengthMinMax(1);
criteria.lifeTime.max = lengthMinMax(2);
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx);
diffAnalysisRes = diffAnalysisRes(indx);
diffModeAnRes = diffModeAnRes(indx);
numNeighbors = numNeighbors(indx);

%save tracksFinal into a new variable name
inputStructure = tracksFinal;

%if tracksFinal is a structure ...
if isstruct(tracksFinal)
    
    %get number of compound tracks
    numCompTrack = length(tracksFinal);
    
    %get the number of segments per compound track
    numTrackSeg = getNumSegmentsPerTrack(tracksFinal);
    
    %determine the location of the 1st segment of each compound track in
    %the "big matrix" of all segments
    segLoc = [1; cumsum(numTrackSeg(1:end-1))+1];
    
    %get the start, end and life time information of each segment
    trajSEL = getTrackSEL(tracksFinal,1);
    
    %get the mean existence time of each segment
    trajMeanTime = mean(trajSEL(:,1:2),2);
    
    %from this get total number of segments (called "traj")
    numTraj = size(trajSEL,1);
    
    %get the start, end and center position of trajectories (segments)
    %also calculate the average frame-to-frame displacement
    trajXCoord = NaN(numTraj,3);
    trajYCoord = NaN(numTraj,3);
    frame2frameDisp = NaN(numTraj,1);
    for iCompTrack = 1 : numCompTrack
        
        %get current compound track's information over its lifetime
        infoCompTrack = tracksFinal(iCompTrack).tracksCoordAmpCG;
        xCoordCompTrack = infoCompTrack(:,1:8:end);
        yCoordCompTrack = infoCompTrack(:,2:8:end);
        
        %get "local" start, end and life times
        segSEL = getTrackSEL(infoCompTrack);
        
        for iSegment = 1 : numTrackSeg(iCompTrack)
            
            %get current segment's positions over its lifetime
            xCoordCurrent = xCoordCompTrack(iSegment,segSEL(iSegment,1):segSEL(iSegment,2));
            yCoordCurrent = yCoordCompTrack(iSegment,segSEL(iSegment,1):segSEL(iSegment,2));
            
            %calculate start, end and center positions
            startPos = [xCoordCurrent(1) yCoordCurrent(1)];
            endPos = [xCoordCurrent(end) yCoordCurrent(end)];
            centerPos = [nanmean(xCoordCurrent) nanmean(yCoordCurrent)];
            
            %assemble the position information, ordered as instructed for the input
            %variable position2plot
            trajXCoord(segLoc(iCompTrack)+iSegment-1,:) = [centerPos(1) startPos(1) endPos(1)];
            trajYCoord(segLoc(iCompTrack)+iSegment-1,:) = [centerPos(2) startPos(2) endPos(2)];
            
            %calculate the average frame-to-frame displacement (root mean
            %square)
            frame2frameDisp(segLoc(iCompTrack)+iSegment-1) = sqrt(nanmean(diff(xCoordCurrent).^2+diff(yCoordCurrent).^2));
            
        end
    end
    
else %if tracksFinal is a matrix ...
    
    %get number of trajectories
    numTraj = size(tracksFinal,1);
    
    %extract the x- and y-coordinates from the big matrix
    xCoord = tracksFinal(:,1:8:end);
    yCoord = tracksFinal(:,2:8:end);
    
    %get the start, end and life time information of trajectories
    trajSEL = getTrackSEL(tracksFinal);
    
    %get the mean existence time of each trajectory
    trajMeanTime = mean(trajSEL(:,1:2),2);
    
    %get the start, end and center position of trajectories
    %also calculate the average frame-to-frame displacement
    trajXCoord = NaN(numTraj,3);
    trajYCoord = NaN(numTraj,3);
    frame2frameDisp = NaN(numTraj,1);
    for iTraj = 1 : numTraj
        
        %get current track's positions over its lifetime
        xCoordCurrent = xCoord(iTraj,trajSEL(iTraj,1):trajSEL(iTraj,2));
        yCoordCurrent = yCoord(iTraj,trajSEL(iTraj,1):trajSEL(iTraj,2));
        
        %calculate start, end and center positions
        startPos = [xCoordCurrent(1) yCoordCurrent(1)];
        endPos = [xCoordCurrent(end) yCoordCurrent(end)];
        centerPos = [nanmean(xCoordCurrent) nanmean(yCoordCurrent)];
        
        %assemble the position information, ordered as instructed for the input
        %variable position2plot
        trajXCoord(iTraj,:) = [centerPos(1) startPos(1) endPos(1)];
        trajYCoord(iTraj,:) = [centerPos(2) startPos(2) endPos(2)];
        
        %calculate the average frame-to-frame displacement (root mean
        %square)
        frame2frameDisp(iTraj) = sqrt(nanmean(diff(xCoordCurrent).^2+diff(yCoordCurrent).^2));
        
    end
    
end

%find x-coordinate limits
minXCoord = min(floor(min(trajXCoord(:))),0);
maxXCoord =  ceil(max(trajXCoord(:)));

%find y-coordinate limits
minYCoord = min(floor(min(trajYCoord(:))),0);
maxYCoord =  ceil(max(trajYCoord(:)));

%% Property extraction and pre-processing

if any(property2plot==1)
    
    %get classifications from MSS diffusion analysis results
    trajClass = vertcat(diffAnalysisRes.classification);
    trajClass = trajClass(:,2);
    
    %calculate the fraction of trajectories in each classification
    fracTrajClass = hist(trajClass,1:3);
    fracTrajClass = fracTrajClass / sum(fracTrajClass);
    
end

if any(property2plot==2)
    
    %get diffusion coefficients from diffusion analysis results
    % diffCoefNorm = catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef');
    diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');
    
    %divide the range of diffusion coefficients into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    % [diffCoefSegment,segmentEdgesDC,fracInSegmentsDC] = divideRangeIntoSegments(...
    %     diffCoefNorm,numSegments,fixedSegLength);
    [diffCoefSegment,segmentEdgesDC,fracInSegmentsDC] = divideRangeIntoSegments(...
        diffCoefGen,numSegments,fixedSegLength);
    
end

if any(property2plot==3)
    
    %get confinement radii from diffusion analysis results
    confRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');
    
    %divide the range of confinement radii into segments, determine which
    %segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [confRadSegment,segmentEdgesCR,fracInSegmentsCR] = divideRangeIntoSegments(...
        confRad,numSegments,fixedSegLength);
    
end

if any(property2plot==4)
    
    %get trajectory lifetimes
    trajLft = trajSEL(:,3);
    
    %divide the range of lifetimes into segments, determine which
    %segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [trajLftSegment,segmentEdgesLft,fracInSegmentsLft] = divideRangeIntoSegments(...
        trajLft,numSegments,fixedSegLength);
    
end

if any(property2plot==5)
    
    %divide the range of frame-to-frame displacements into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [f2fDispSegment,segmentEdgesFFD,fracInSegmentsFFD] = divideRangeIntoSegments(...
        frame2frameDisp,numSegments,fixedSegLength);
    
end

if any(property2plot==6)
    
    %get diffusion modes from diffusion mode analysis results
    trajDiffMode = vertcat(diffModeAnRes.diffMode);
    
    %calculate the fraction of trajectories in each mode
    fracDiffMode = hist(trajDiffMode,1:4);
    fracDiffMode = fracDiffMode / sum(fracDiffMode);
    
end

if any(property2plot==7)
    
    %get diffusion coefficients from diffusion mode analysis results
    diffCoefDisp2 = vertcat(diffModeAnRes.diffCoef);
    
    %divide the range of diffusion coefficients into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [diffCoef2Segment,segmentEdgesDC2,fracInSegmentsDC2] = divideRangeIntoSegments(...
        diffCoefDisp2,numSegments,fixedSegLength);
    
end

if any(property2plot==8)
    
    %get number of neighbors from numNeighbors variable
    smDensity = vertcat(numNeighbors.value);
    
    %divide the range of number of neighbors into segements, determine which segment
    %each trajectory falls into, and calculate the fraction of trajectories
    %in each segement
    [smDensitySegment,segmentEdgesDensity,fracInSegmentsDensity] = divideRangeIntoSegments(...
        smDensity,numSegments,fixedSegLength);
    
end

% %% Some numbers for output
% 
% %get indices of trajectories in the various classifications
% indxConf = find(trajClass==1);
% indxBrown = find(trajClass==2);
% indxDir = find(trajClass==3);
% 
% f2fDispRMS = [nanmean(frame2frameDisp(indxConf)) nanmean(frame2frameDisp(indxBrown)) ...
%     nanmean(frame2frameDisp(indxDir)) nanmean(frame2frameDisp([indxConf;indxBrown;indxDir])) ...
%     nanmean(frame2frameDisp)];

%% Movie

%get property and position index
iProperty = property2plot;
iPos = position2plot;

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
for iFrame = 1 : numFramesMovie-1
    
    clf;
    
    %read image + mask
    imageStack = imread(fullfile(imageDir,imageFileListing(iFrame).name));
    maskStack1 = imread(fullfile(maskDir,maskFileListing(iFrame).name));
    maskStack2 = imread(fullfile(maskDir,maskFileListing(iFrame+1).name));

    %plot cell image + mask
    axes('Position',[0 0 0.495 1]);
    imshow(imageStack,[prctile(imageStack(:),5) prctile(imageStack(:),95)]);
    %     imshow(imageStack,[prctile(imageStack(:),25) prctile(imageStack(:),75)]);
    hold on;
    maskBounds = bwboundaries(maskStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    maskBounds = bwboundaries(maskStack2);
    cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',2)),maskBounds);
    maskBounds = bwboundaries(maskStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
    
    textDeltaCoord = min(diff(imageRange,[],2))/20;
    %     text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
    %         textDeltaCoord,num2str(iFrame),'Color','white');
    text(imageRange(1,1)+textDeltaCoord-1,imageRange(2,1)+...
        textDeltaCoord+2,[num2str((iFrame-1)*10) ' s'],'Color','white'); %,'FontSize',20);
    
    legend1 = legend('Current frame','Next frame');
    set(legend1,'TextColor',[1 1 1],'Box','Off','Fontsize',8);
    
    %make space for plotting particles + mask
    axes('Position',[0.505 0 0.495 1]);
    imshow(ones(isx,isy));
    hold on
    
    if iFrame < numFramesMovie
        
        %read mask of next image
        maskStack2 = imread(fullfile(maskDir,maskFileListing(iFrame+1).name));
        
        %collect particle positions
        indxTime = find( trajMeanTime>=(iFrame-1)*numFramesSPT+1 & trajMeanTime < iFrame*numFramesSPT );
        
        %initialize figure legend text
        legendText = [];
        
        %plot particles color-coded based on property
        switch iProperty
            
            %             case 1 %classification
            %
            %                 %unclassified
            %                 indx = find(isnan(trajClass));
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color',[0.7 0.7 0.7]);
            %                     legendText{end+1} = 'unclassified'; %#ok<AGROW>
            %                 end
            %
            %                 %free diffusion
            %                 indx = find(trajClass == 2);
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color','c');
            %                 end
            %
            %                 %confined
            %                 indx = find(trajClass == 1);
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color','b');
            %                 end
            %
            %                 %directed
            %                 indx = find(trajClass == 3);
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color','m');
            %                 end
            %
            %             case 2 %diffusion coefficient MSS
            %
            %                 %trajectories without a diffusion coefficient i.e.
            %                 %unclassified trajectories
            %                 indx = find(isnan(diffCoefSegment));
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color',[0.7 0.7 0.7]);
            %                     legendText{end+1} = 'unclassified'; %#ok<AGROW>
            %                 end
            %
            %                 %go over the different diff. coef. segments
            %                 for iSegment = 1 : numSegments
            %                     indx = find(diffCoefSegment==iSegment);
            %                     if ~isempty(indx)
            %                         plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                             '.','Color',segmentColor(iSegment,:));
            %                     end
            %                 end
            %
            %             case 3 %confinement radius
            %
            %                 %unclassified trajectories
            %                 indx = find(isnan(confRadSegment) & isnan(trajClass));
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color',[0.7 0.7 0.7]);
            %                     legendText{end+1} = 'unclassified'; %#ok<AGROW>
            %                 end
            %
            %                 %trajectories that are not confined
            %                 indx = find(isnan(confRadSegment) & ~isnan(trajClass));
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '+','Color',[0.7 0.7 0.7]);
            %                     legendText{end+1} = 'not confined'; %#ok<AGROW>
            %                 end
            %
            %                 %go over the different conf. rad. segment
            %                 for iSegment = 1 : numSegments
            %                     indx = find(confRadSegment==iSegment);
            %                     if ~isempty(indx)
            %                         plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                             '.','Color',segmentColor(iSegment,:));
            %                     end
            %                 end
            %
            %             case 4 %lifetime
            %
            %                 %go over the different lifetime segments
            %                 for iSegment = 1 : numSegments
            %                     indx = find(trajLftSegment==iSegment);
            %                     if ~isempty(indx)
            %                         plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                             '.','Color',segmentColor(iSegment,:));
            %                     end
            %                 end
            %
            %             case 5 %frame-to-frame displacement
            %
            %                 %trajectories without a frame-to-frame displacement
            %                 indx = find(isnan(f2fDispSegment));
            %                 if ~isempty(indx)
            %                     plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                         '.','Color',[0.7 0.7 0.7]);
            %                     legendText{end+1} = 'unclassified'; %#ok<AGROW>
            %                 end
            %
            %                 %go over the different frame-to-frame displacement segments
            %                 for iSegment = 1 : numSegments
            %                     indx = find(f2fDispSegment==iSegment);
            %                     if ~isempty(indx)
            %                         plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
            %                             '.','Color',segmentColor(iSegment,:));
            %                     end
            %                 end
                
            case 6 %diffusion mode
                
                %get current particles' information
                trajDiffModeTime = trajDiffMode(indxTime);
                trajXCoordTime = trajXCoord(indxTime,:);
                trajYCoordTime = trajYCoord(indxTime,:);
                
                %first plot modes in correct order for legend
                
                %mode 1
                indx = find(trajDiffModeTime == 1);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','k','LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'Immobile'; %#ok<AGROW>
                end
                
                %mode 2
                indx = find(trajDiffModeTime == 2);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','b','LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'Slow'; %#ok<AGROW>
                end
                
                %mode 3
                indx = find(trajDiffModeTime == 3);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','c','LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'Fast'; %#ok<AGROW>
                end
                
                %mode 4
                indx = find(trajDiffModeTime == 4);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','m','LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'Very fast'; %#ok<AGROW>
                end
                
                %unclassified
                indx = find(isnan(trajDiffModeTime));
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'Unclassified'; %#ok<AGROW>
                end
                
                %then plot modes in best order for display
                
                %mode 3
                indx = find(trajDiffModeTime == 3);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','c','LineWidth',2,'MarkerSize',1.5);
                end
                
                %mode 2
                indx = find(trajDiffModeTime == 2);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','b','LineWidth',2,'MarkerSize',1.5);
                end
                
                %mode 4
                indx = find(trajDiffModeTime == 4);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','m','LineWidth',2,'MarkerSize',1.5);
                end
                
                %mode 1
                indx = find(trajDiffModeTime == 1);
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color','k','LineWidth',2,'MarkerSize',1.5);
                end
                
            case 7 %diffusion coefficient from mode analysis
                
                %get current particles' information
                diffCoef2SegmentTime = diffCoef2Segment(indxTime);
                trajXCoordTime = trajXCoord(indxTime,:);
                trajYCoordTime = trajYCoord(indxTime,:);
                
                %trajectories without a diffusion coefficient i.e.
                %unclassified trajectories
                indx = find(isnan(diffCoef2SegmentTime));
                if ~isempty(indx)
                    plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                        'o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerSize',1.5);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %go over the different diff. coef. segments
                for iSegment = 1 : numSegments
                    indx = find(diffCoef2SegmentTime==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                            'o','Color',segmentColor(iSegment,:),'LineWidth',2,'MarkerSize',1.5);
                    end
                end
                
            case 8 %number of neighbors
                
                %get current particles' information
                smDensitySegmentTime = smDensitySegment(indxTime);
                trajXCoordTime = trajXCoord(indxTime,:);
                trajYCoordTime = trajYCoord(indxTime,:);
                
                %go over the different density segments
                for iSegment = 1 : numSegments
                    indx = find(smDensitySegmentTime==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoordTime(indx,iPos),trajYCoordTime(indx,iPos),...
                            'o','Color',segmentColor(iSegment,:),'LineWidth',2,'MarkerSize',1.5);
                    end
                end
                
        end
        
        %overlay masks
        maskBounds = bwboundaries(maskStack1);
        cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
        maskBounds = bwboundaries(maskStack2);
        cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',2)),maskBounds);
        maskBounds = bwboundaries(maskStack1);
        cellfun(@(x)(plot(x(:,2),x(:,1),'g','LineWidth',2)),maskBounds);
        
    end
    
    %     %add subplot title
    %     title([plottedProperty{iProperty} ' vs. ' plottedPosition{iPos}])
    
    %write color legend
    %     legend1 = legend(legendText);
    %     set(legend1,'Box','off','Fontsize',8);
    
%     text(128,imageRange(2,1)+...
%         textDeltaCoord+2,'Immobile','Color','black'); %,'FontSize',20);
%     text(190,imageRange(2,1)+...
%         textDeltaCoord+2,'Slow','Color','blue'); %,'FontSize',20);
%     text(128,imageRange(2,1)+...
%         textDeltaCoord+17,'Fast','Color','cyan'); %,'FontSize',20);
%     text(170,imageRange(2,1)+...
%         textDeltaCoord+17,'Very fast','Color','magenta'); %,'FontSize',20);
%     text(128,imageRange(2,1)+...
%         textDeltaCoord+32,'Unclassified','Color',[0.7 0.7 0.7]); %,'FontSize',20);
    
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

% %hold off
% hold off
% 
% % SUBPLOT 2: Distribution of values %
% 
% subplot(1,3,3)
% hold on
% 
% %label axes
% xlabel(plottedProperty{iProperty});
% ylabel('Fraction');
% 
% %initialize figure legend text
% legendText = [];
% 
% switch iProperty
%     
%     case 1 %classification
%         
%         %plot the bars with different colors
%         bar([1 2],[fracTrajClass(1) 0],'BarWidth',1,'FaceColor',...
%             'b','EdgeColor','none');
%         bar([2 3],[fracTrajClass(2) 0],'BarWidth',1,'FaceColor',...
%             'c','EdgeColor','none');
%         bar([3 4],[fracTrajClass(3) 0],'BarWidth',1,'FaceColor',...
%             'm','EdgeColor','none');
%         legendText = {'confined','free','directed'};
%         
%     case 2 %diffusion coefficient MSS
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesDC,2);
%         segmentWidth = segmentEdgesDC(:,2) - segmentEdgesDC(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsDC(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesDC(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesDC(iSegment,2))]; %#ok<AGROW>
%         end
%         
%     case 3 %confinement radius
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesCR,2);
%         segmentWidth = segmentEdgesCR(:,2) - segmentEdgesCR(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsCR(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesCR(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesCR(iSegment,2))]; %#ok<AGROW>
%         end
%         
%     case 4 %lifetime
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesLft,2);
%         segmentWidth = segmentEdgesLft(:,2) - segmentEdgesLft(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsLft(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesLft(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesLft(iSegment,2))]; %#ok<AGROW>
%         end
%         
%     case 5 %frame-to-frame displacement
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesFFD,2);
%         segmentWidth = segmentEdgesFFD(:,2) - segmentEdgesFFD(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsFFD(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesFFD(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesFFD(iSegment,2))]; %#ok<AGROW>
%         end
%         
%     case 6 %diffusion mode
%         
%         %plot the bars with different colors
%         bar([1 2],[fracDiffMode(1) 0],'BarWidth',1,'FaceColor',...
%             'k','EdgeColor','none');
%         bar([2 3],[fracDiffMode(2) 0],'BarWidth',1,'FaceColor',...
%             'b','EdgeColor','none');
%         bar([3 4],[fracDiffMode(3) 0],'BarWidth',1,'FaceColor',...
%             'c','EdgeColor','none');
%         bar([4 5],[fracDiffMode(4) 0],'BarWidth',1,'FaceColor',...
%             'm','EdgeColor','none');
%         legendText = {'mode 1','mode 2','mode 3','mode 4'};
%         
%     case 7 %diffusion coefficient from mode analysis
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesDC2,2);
%         segmentWidth = segmentEdgesDC2(:,2) - segmentEdgesDC2(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsDC2(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesDC2(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesDC2(iSegment,2))]; %#ok<AGROW>
%         end
%         
%     case 8 %number of neighbors
%         
%         %get the center and width of each segment
%         segmentCenter = mean(segmentEdgesDensity,2);
%         segmentWidth = segmentEdgesDensity(:,2) - segmentEdgesDensity(:,1);
%         
%         %plot the bars with different colors
%         for iSegment = 1 : numSegments
%             bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
%                 segmentWidth(iSegment)],[fracInSegmentsDensity(iSegment) 0],...
%                 'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
%                 'EdgeColor','none');
%             legendText{end+1} = [num2str(segmentEdgesDensity(iSegment,1)) ...
%                 ' - ' num2str(segmentEdgesDensity(iSegment,2))]; %#ok<AGROW>
%         end
%         
% end
% 
% %add subplot title
% title(['Distribution of ' plottedProperty{iProperty} ' values'])
% 
% %write color legend - this won't include unclassified/unconfined
% legend(legendText)
% 
% %hold off
% hold off
% 
% %save figure
% saveas(h,[figureName '_' plottedProperty{iProperty} '_vs_' plottedPosition{iPos}],'fig');
% 
% 



%% subfunctions

function [valVectorSegment,segmentEdges,fracInSegment] = divideRangeIntoSegments(...
    valVector,numSegments,fixedSegLength)

%find minimum and maximum values
minVal = min(valVector);
maxVal = max(valVector);

if fixedSegLength %if fixed segment length ...
    
    %divide the range of values into equal segments
    valIncrement = (maxVal - minVal) / numSegments;
    
    %define the segment upper edges
    segmentUpperEdges = minVal + (1 : numSegments) * valIncrement;
    
else %if equal distribution of elements among segments ...
    
    %determine segment upper edges based on percentiles
    segmentUpperEdges = prctile(valVector,(1:numSegments)*100/numSegments);
    
end

%determine which segment each value falls into
valVectorSegment = NaN(size(valVector));
for iSegment = numSegments : -1 : 1
    valVectorSegment(valVector<=segmentUpperEdges(iSegment)) = iSegment;
end

%store the segment edges
segmentEdges = [[minVal; segmentUpperEdges(1:end-1)'] segmentUpperEdges'];

%calculate the fraction of values in each segment
fracInSegment = hist(valVectorSegment,1:numSegments);
fracInSegment = fracInSegment/sum(fracInSegment);


% function segmentColor = distributeSegmentColors(numSegments)
%
% %define the color regimes - fixed to 5 for now
% numColors = 5;
% regimeEdgeColor = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1]; %bcgyrm
%
% %generate vector of segment indices
% segmentIndx = 1 : numSegments;
% 
% %initialize the vector of segment colors
% segmentColor = zeros(numSegments,3);
% 
% %go over color regimes and assign segments to each color regime
% for iColor = 1 : numColors
% 
%     %divide number of remaining segments by number of remaining colors
%     numSegments4thisColor = length(segmentIndx) / (numColors - iColor + 1);
%     
%     %assign segments to this color regime
%     segments4thisColor = segmentIndx(1:numSegments4thisColor);
%     
%     %assign colors to the segments of this color regime
%     for iSegment = 1 : numSegments4thisColor
%         segmentColor(segments4thisColor(iSegment),:) = ...
%             regimeEdgeColor(iColor,:) * ...
%             (numSegments4thisColor-iSegment+1)/numSegments4thisColor ...
%             + regimeEdgeColor(iColor+1,:) * (iSegment-1)/numSegments4thisColor;
%     end
% 
%     %update vector of segment indices by removing the indices already
%     %assigned to this color
%     segmentIndx = segmentIndx(numSegments4thisColor+1:end);
%     
% end

