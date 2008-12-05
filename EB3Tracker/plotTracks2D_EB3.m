function [selectedTracks] = plotTracks2D_EB3(trackedFeatureInfo,timeRange,...
    newFigure,img,flipXY,ask4sel,plotCurrentOnly,roiYX,movieInfo)
%PLOTTRACKS2D plots a group of tracks in 2D and allows user to click on them and extract track information
%
%SYNOPSIS plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure,img,flipXY,ask4sel)
%
%INPUT  trackedFeatureInfo: either tracksFinal structure,
%                           or trackedFeatureInfo matrix
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row
%                           consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%       timeRange         : 2-element row vector indicating time range to plot.
%                           Optional. Default: whole movie.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       img                : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%       flipXY            : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%       ask4sel           : 1 if user should be asked to select tracks in
%                           plot in order to show track information.
%                           Optional. Default: 1.
%       plotCurrentOnly   : number of a frame, will only plot those tracks
%                           between timeRange that exist during that frame
%       roiYX             : coordinates of a region-of-interest (closed
%                           polygon), the rectangle circumscribing the ROI
%                           will be used to limit plotting to within the
%                           region
%
%OUTPUT The plot.
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(trackedFeatureInfo)
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(trackedFeatureInfo);
    clear trackedFeatureIndx
end

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;

selectedTracks=[];

%check whether a time range for plotting was input
if nargin < 2 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        error('--plotTracks2D: Wrong time range for plotting!');
    end
end

%check whether newFigure was input
if nargin < 3 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        error('--plotTracks2D: newFigure should be 0 or 1!');
    end
end

%check whether user supplied an image
if nargin < 4 || isempty(img)
    img = [];
end

%check whether detection coordinates should be flipped for plotting
if nargin < 5 || isempty(flipXY)
    flipXY = false;
end

%check whether user wants to get track data by selection
if nargin < 6 || isempty(ask4sel)
    ask4sel = true;
end

%check whether user wants to only plot "active" tracks at
%plotCurrentOnly-frame or all tracks in frame range
if nargin<7 || isempty(plotCurrentOnly)
    plotCurrentOnly=0; %default: all frames in timeRange
end

%check for coordinates of a ROI in which results should be plotted
if nargin<8 || isempty(roiYX) || isequal(roiYX,0)
    [imL imW]=size(img);
    minY=1;
    maxY=imL;
    minX=1;
    maxX=imW;
else
    minY=floor(min(roiYX(:,1)));
    maxY=ceil(max(roiYX(:,1)));
    minX=floor(min(roiYX(:,2)));
    maxX=ceil(max(roiYX(:,2)));
end

if nargin<9 || ~isstruct(movieInfo)
    %no detection results to plot
    movieInfo=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% trackInfo contains gap start/end info, while trackVelocities has
% per-frame or average-per-segment track veloctiies
[trackedFeatureInfoInterp,trackInfo,trackVelocities] = getVelocitiesFromMat(trackedFeatureInfo,3);

%get the x,y-coordinates of features in all tracks
if flipXY
    tracksY = trackedFeatureInfo(:,1:8:end)';
    tracksX = trackedFeatureInfo(:,2:8:end)';
    tracksYInterp = trackedFeatureInfoInterp(:,1:8:end)';
    tracksXInterp = trackedFeatureInfoInterp(:,2:8:end)';
else
    tracksX = trackedFeatureInfo(:,1:8:end)';
    tracksY = trackedFeatureInfo(:,2:8:end)';
    tracksXInterp = trackedFeatureInfoInterp(:,1:8:end)';
    tracksYInterp = trackedFeatureInfoInterp(:,2:8:end)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%open new figure window
figure(gcf)

if ~isempty(img)
    img=img(minY:maxY,minX:maxX);
    imagesc(img)
    colormap gray
else
    %find coordinate limits
    maxXCoord =  ceil(nanmax(tracksX(:)));
    maxYCoord =  ceil(nanmax(tracksY(:)));
    img = repmat(0.75*ones(maxYCoord,maxXCoord),[1,1,3]);
    imagesc(img); %plot an empty image
    colormap gray
end

%show coordinates on axes
ah = gca;
set(ah,'visible','on');

%label axes
if flipXY
    xlabel('y-coordinate (pixels)');
    ylabel('x-coordinate (pixels)');
else
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');
end


%hold on figure
hold on


%concatenate all segments and gaps into n x 4 matrices
%each has format [trackNum startFrame endFrame velocity]
segs=(vertcat(trackInfo.seg));
fgaps=(vertcat(trackInfo.fgap));
bgaps=(vertcat(trackInfo.bgap));
ugaps=(vertcat(trackInfo.ugap));

if isempty(segs)
    segs=zeros(1,4);
end
if isempty(fgaps)
    fgaps=zeros(1,4);
end
if isempty(bgaps)
    bgaps=zeros(1,4);
end
if isempty(ugaps)
    ugaps=zeros(1,4);
end


%get list of all tracks that have begun before plotCurrentOnly-frame but
%finish after
if plotCurrentOnly~=0
    plotCurrentTracksIdx=unique([segs(segs(:,2)<=plotCurrentOnly & segs(:,3)>=plotCurrentOnly,1);...
        fgaps(fgaps(:,2)<=plotCurrentOnly & fgaps(:,3)>=plotCurrentOnly,1);...
        bgaps(bgaps(:,2)<=plotCurrentOnly & bgaps(:,3)>=plotCurrentOnly,1);...
        ugaps(ugaps(:,2)<=plotCurrentOnly & ugaps(:,3)>=plotCurrentOnly,1)]);
else
    plotCurrentTracksIdx=1:numTracks;
end
%these ones don't meet that criterion
removeTheseTracks=setdiff(1:numTracks,plotCurrentTracksIdx);

%fill matrices with coordinates for segments
segMatX = NaN(numTimePoints,size(segs,1));
segMatY = NaN(numTimePoints,size(segs,1));
for i=1:size(segs,1)
    if isempty(intersect(segs(i,1),removeTheseTracks))
        segMatX(segs(i,2):segs(i,3),i) = tracksXInterp(segs(i,2):segs(i,3),segs(i,1));
        segMatY(segs(i,2):segs(i,3),i) = tracksYInterp(segs(i,2):segs(i,3),segs(i,1));
    end
end

%fill matrices with coordinates for forward gaps
fgapMatX = NaN(numTimePoints,size(fgaps,1));
fgapMatY = NaN(numTimePoints,size(fgaps,1));
for i=1:size(fgaps,1)
    if isempty(intersect(fgaps(i,1),removeTheseTracks)) && sum(fgaps(:))~=0
        fgapMatX(fgaps(i,2):fgaps(i,3),i) = tracksXInterp(fgaps(i,2):fgaps(i,3),fgaps(i,1));
        fgapMatY(fgaps(i,2):fgaps(i,3),i) = tracksYInterp(fgaps(i,2):fgaps(i,3),fgaps(i,1));
    end
end

%fill matrices with coordinates for backward gaps
bgapMatX = NaN(numTimePoints,size(bgaps,1));
bgapMatY = NaN(numTimePoints,size(bgaps,1));
for i=1:size(bgaps,1)
    if isempty(intersect(bgaps(i,1),removeTheseTracks)) && sum(bgaps(:))~=0
        bgapMatX(bgaps(i,2):bgaps(i,3),i) = tracksXInterp(bgaps(i,2):bgaps(i,3),bgaps(i,1));
        bgapMatY(bgaps(i,2):bgaps(i,3),i) = tracksYInterp(bgaps(i,2):bgaps(i,3),bgaps(i,1));
    end
end

%fill matrices with coordinates for unclassified gaps
ugapMatX = NaN(numTimePoints,size(ugaps,1));
ugapMatY = NaN(numTimePoints,size(ugaps,1));
for i=1:size(ugaps,1)
    if isempty(intersect(ugaps(i,1),removeTheseTracks))  && sum(ugaps(:))~=0
        ugapMatX(ugaps(i,2):ugaps(i,3),i) = tracksXInterp(ugaps(i,2):ugaps(i,3),ugaps(i,1));
        ugapMatY(ugaps(i,2):ugaps(i,3),i) = tracksYInterp(ugaps(i,2):ugaps(i,3),ugaps(i,1));
    end
end

%if ROI limits were given, make null all coordinates outside the ROI and
%transform coordinates to roi-frame of reference. if a segment or gap goes
%outside the roi but comes back in, part of the track will be lost visually
%since the plot function doesn't deal with nans.  this could be fixed by
%looping through to remove missing points as khuloud does in plotTracks2D,
%but here we don't do that to achieve much faster plotting

outOfRangeIdx=find(segMatX(:)<minX | segMatX(:)>maxX | segMatY(:)<minY | segMatY(:)>maxY);
segMatX(outOfRangeIdx)=nan;
segMatY(outOfRangeIdx)=nan;
segMatX=segMatX-minX+1;
segMatY=segMatY-minY+1;

outOfRangeIdx=find(fgapMatX(:)<minX | fgapMatX(:)>maxX | fgapMatY(:)<minY | fgapMatY(:)>maxY);
fgapMatX(outOfRangeIdx)=nan;
fgapMatY(outOfRangeIdx)=nan;
fgapMatX=fgapMatX-minX+1;
fgapMatY=fgapMatY-minY+1;

outOfRangeIdx=find(bgapMatX(:)<minX | bgapMatX(:)>maxX | bgapMatY(:)<minY | bgapMatY(:)>maxY);
bgapMatX(outOfRangeIdx)=nan;
bgapMatY(outOfRangeIdx)=nan;
bgapMatX=bgapMatX-minX+1;
bgapMatY=bgapMatY-minY+1;

outOfRangeIdx=find(ugapMatX(:)<minX | ugapMatX(:)>maxX | ugapMatY(:)<minY | ugapMatY(:)>maxY);
ugapMatX(outOfRangeIdx)=nan;
ugapMatY(outOfRangeIdx)=nan;
ugapMatX=ugapMatX-minX+1;
ugapMatY=ugapMatY-minY+1;




hold on
plot(ugapMatX(timeRange(1):timeRange(2),:),ugapMatY(timeRange(1):timeRange(2),:),'m:','LineWidth',2)
plot(fgapMatX(timeRange(1):timeRange(2),:),fgapMatY(timeRange(1):timeRange(2),:),'c:','LineWidth',2)
plot(bgapMatX(timeRange(1):timeRange(2),:),bgapMatY(timeRange(1):timeRange(2),:),'y:','LineWidth',2)
plot(segMatX(timeRange(1):timeRange(2),:),segMatY(timeRange(1):timeRange(2),:),'r','LineWidth',2)

if ~isempty(movieInfo)
    colorOverTime = jet(timeRange(2)-timeRange(1)+1);
    frmCount1=1;
    for j=timeRange(1):timeRange(2)
        xCoord = vertcat(movieInfo(j).xCoord); xCoord = xCoord(:,1);
        yCoord = vertcat(movieInfo(j).yCoord); yCoord = yCoord(:,1);
        outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
        xCoord(outOfRangeIdx) = [];
        yCoord(outOfRangeIdx) = [];
        xCoord=xCoord-minX+1;
        yCoord=yCoord-minY+1;
        hold on
        scatter(xCoord,yCoord,'Marker','.','cData',repmat(colorOverTime(frmCount1,:),[length(xCoord),1]));
        frmCount1=frmCount1+1;
    end
end


if ask4sel

    %extract the portion of tracksX and tracksY that is of interest
    tracksXPInterp = tracksXInterp(timeRange(1):timeRange(2),:)-minX+1;
    tracksYPInterp = tracksYInterp(timeRange(1):timeRange(2),:)-minY+1;

    %ask the user whether to click on figure and get frame information
    userEntry = input('select points in figure? y/n ','s');
    if strcmp(userEntry,'y')
        disp('Use normal button clicks to add points.')
        disp('Double-click to add a final point and end selection.')
        disp('Press Enter to end selection without adding a final point.')
        disp('Press Backspace or Delete to remove previously selected point.')
        disp('--------------------------------------------------------------')
    end
    count = 1;
    while strcmp(userEntry,'y')

        %let the user choose the points of interest
        [x,y] = getpts;

        %find the time points of the indicated points
        for i=1:length(x)

            %find the distances between those points and the tracks
            distTrack2Point = (tracksXPInterp-x(i)).^2+(tracksYPInterp-y(i)).^2;

            %determine the minimum distance for each chosen point
            [frameChosen,rowChosen] = find(distTrack2Point==min(distTrack2Point(:)));

            %go over all chosen rows
            for j = 1 : length(rowChosen)

                %find the track corresponding to each minimum distance
                trackChosen = find(1:numTracks <= rowChosen(j),1,'last');

                tSeg  =  segs( segs(:,1)==trackChosen,:);
                tFgap = fgaps(fgaps(:,1)==trackChosen,:);
                tBgap = bgaps(bgaps(:,1)==trackChosen,:);
                tUgap = ugaps(ugaps(:,1)==trackChosen,:);


                disp(['Track: ' num2str(trackChosen) ...
                    '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
                    '   Coordinates: ' num2str(tracksXPInterp(frameChosen(j),rowChosen(j))) ...
                    ' ' num2str(tracksYPInterp(frameChosen(j),rowChosen(j)))]);

                prof = [tSeg 1*ones(size(tSeg,1),1); ...
                    tFgap 2*ones(size(tFgap,1),1); ...
                    tBgap 3*ones(size(tBgap,1),1); ...
                    tUgap 4*ones(size(tUgap,1),1);];
                trackProfile = sortrows(prof,2)
                selectedTracks{count,1} = trackProfile;
                
                text(tracksXPInterp(frameChosen(j),rowChosen(j)),tracksYPInterp(frameChosen(j),rowChosen(j)),...
                    ['\leftarrow' num2str(count)],'Color','y','FontWeight','bold')
                count = count+1;
            end

        end

        %ask the user again whether to click on figure and get frame information
        userEntry = input('select points again? y/n ','s');

    end

end

%%%%% ~~ the end ~~ %%%%%

