function [selectedTracks] = plusTipPlotTracks(projData,subIdx,timeRange,img,ask4sel,plotCurrentOnly,roiYX,movieInfo,rawToo)
%plusTipPlotTracks overlays tracks on an image and allows user to select individual tracks
%
%INPUT  projData         : structure containing track info (stored in meta
%                          folder); if not given, user will be asked to
%                          select the file
%       subIdx           : indices for subtracks that should be considered
%                          for plotting, corresponds to rows of projData
%                          .nTrack_start_end_velMicPerMin_class_lifetime
%                          if given as [], all subtracks will be
%                          considered
%       timeRange        : row vector of the form [startFrame endFrame]
%                          indicating time range to plot. if not given or
%                          given as [], tracks from the whole movie will
%                          be displayed
%       img              : image for overlay. if not given or given as [],
%                          user will be asked to select the image.
%       ask4sel          : 1 (default) if user should be asked to select
%                          tracks in plot in order to show track
%                          information, 0 otherwise
%       plotCurrentOnly  : a particular frame number within timeRange
%                          if given, will only plot those tracks
%                          that exist during that frame.
%       roiYX            : coordinates of a region-of-interest (closed
%                          polygon) in which tracks should be plotted, or
%                          a logical mask the size of the image.  the
%                          rectangle circumscribing the ROI will be used
%                          to limit plotting to within the
%                          region. if not given or given as [], whole
%                          image will be used.
%       movieInfo        : structure containing detection results. can be
%                          loaded from "feat" directory. if given, ALL
%                          detected features are plotted with color
%                          corresponding to time (startFrame in dark blue
%                          to endFrame in deep red). if not given or given
%                          as [], detected features will not be plotted.
%                          note that all features, not just the ones
%                          accepted by the tracking, will be plotted.
%       rawToo            : 1  to show raw image on left of overlay
%                          {0} to make overlay without dual panel
%
%OUTPUT The plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectedTracks=[]; % initialize output

% load projData
if nargin<1 || isempty(projData)
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    if isequal(fileName,0)
        return
    end
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
end
projData.anDir=formatPath(projData.anDir);
projData.imDir=formatPath(projData.imDir);

% get number of time points
numFrames=projData.numFrames;

% check whether specific subtrack indices were given
if nargin<2 || isempty(subIdx)
    subIdx=[];
end

% check whether a time range for plotting was input
if nargin<3 || isempty(timeRange)
    timeRange=[1 numFrames];
else
    if timeRange(1)<1
        timeRange(1)=1;
    end
    if timeRange(2)>numFrames
        timeRange(2)=numFrames;
    end
end

% check whether user supplied an image
if nargin<4 || isempty(img)
    % let user pick image
    homeDir=pwd;
    cd(projData.imDir)
    [fileName,pathName] = uigetfile('*.tif','Select image to use for track overlay');
    if isequal(fileName,0)
        return
    end
    img=double(imread([pathName filesep fileName]));
    cd(homeDir);
end

% check whether user wants to get track data by selection
if nargin<5 || isempty(ask4sel)
    ask4sel=true;
end

% check whether user wants to only plot "active" tracks at
% plotCurrentOnly-frame or all tracks in frame range
if nargin<6 || isempty(plotCurrentOnly)
    plotCurrentOnly=0; % default: all frames in timeRange
    figure
else
    if plotCurrentOnly<timeRange(1) || plotCurrentOnly>timeRange(2)
        error('--plusTipPlotTracks: plotCurrentOnly must be within timeRange');
    end
end

% check input for ROI coordinates
if nargin<7 || isempty(roiYX)
    [r c]=find(ones(size(img)));
    roiYX=[min(r) min(c); max(r) max(c)];
elseif islogical(roiYX)
    [r c]=find(roiYX);
    roiYX=[min(r) min(c); max(r) max(c)];
else
    if size(roiYX,2)~=2
        error('--plusTipPlotTracks: roiYX should be nx2 matrix of coordinates or logical mask')
    end
end

if nargin<8 || isempty(movieInfo)
    movieInfo=[];
end

if nargin<9 || isempty(rawToo)
    rawToo=0;
end

% get track information
trackData=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

if isempty(subIdx)
    fullIdx=trackData(:,1);
    trackType=trackData(:,5);
    sF=trackData(:,2);
else
    fullIdx=trackData(subIdx,1);
    trackType=trackData(subIdx,5);
    sF=trackData(subIdx,2);
end

% find bounds of the image to use
minY=floor(min(roiYX(:,1)));
maxY=ceil(max(roiYX(:,1)));
minX=floor(min(roiYX(:,2)));
maxX=ceil(max(roiYX(:,2)));
img=img(minY:maxY,minX:maxX);

if rawToo==1
    img=[img img];
end

% nSubtrack x nFrames matrices of coordinates for all subtracks in subIdx
[xMat,yMat]=plusTipGetSubtrackCoords(projData,subIdx);

% if a subtrack doesn't begin before and end after plotCurrentOnly frame,
% then it doesn't get plotted
if plotCurrentOnly~=0
    correspFullIdx=fullIdx(~isnan(xMat(:,plotCurrentOnly)));
    if ~isempty(correspFullIdx)
        subtracks2keep=find(ismember(fullIdx,correspFullIdx));
        xMat=xMat(subtracks2keep,:);
        yMat=yMat(subtracks2keep,:);
        fullIdx=fullIdx(subtracks2keep);
        trackType=trackType(subtracks2keep);
        sF=sF(subtracks2keep);
    else
        xMat=[];
        yMat=[];
    end
end

%open new figure window
gcf;
imagesc(img);
colormap gray
%axis equal

%show coordinates on axes
ah = gca;
set(ah,'visible','on');

%label axes
xlabel('x-coordinate (pixels)');
ylabel('y-coordinate (pixels)');

if isempty(xMat)
    return
end


% any coordinates outside the ROI's rectangle don't get plotted. if a
% subtrack goes outside the rectangle but comes back in, part of the track
% will be lost visually since the plot function doesn't deal with nans.
xMat(xMat(:)<minX | xMat(:)>maxX)=nan;
yMat(yMat(:)<minY | yMat(:)>maxY)=nan;
% transform coords to rectangle frame of reference
xMat=xMat-minX+1;
yMat=yMat-minY+1;

% only the parts of the subtracks within timeRange will be plotted
frames2Plot=nan(size(xMat));
frames2Plot(:,timeRange(1):timeRange(2))=1;
xMat=xMat.*frames2Plot;
yMat=yMat.*frames2Plot;

if rawToo==1
    xMat=xMat+length(img(1,:))/2;
    % y values don't change
end

%hold on figure
hold on

for t=1:5
    switch t
        case 1 % growth
            c1='r';
            c2='r';
        case 2 % forward gaps - pause
            c1='c:';
            c2='c';
        case 3 % backward gaps - shrinkage
            c1='y:';
            c2='y';
        case 4 % unclassified gaps
            c1='m:';
            c2='m';
        case 5 % forward gaps - reclassified as growth
            c1='g';
            c2='g';
        case 6 % backward gaps - reclassified as pause
            c1='b:';
            c2='b';
      
    end

    plot(xMat(trackType==t,timeRange(1):timeRange(2))',...
        yMat(trackType==t,timeRange(1):timeRange(2))',c1,'LineWidth',1)

    % plot big circle around transition events (red circle when start of
    % growth, yellow circle when start of shrinkage, etc.)
    if plotCurrentOnly~=0
        cTTIdx=find(trackType==t); % current track type indices
        currenStartingSubtrackIdx=cTTIdx(sF(cTTIdx)==plotCurrentOnly);
        if ~isempty(currenStartingSubtrackIdx)
            coordIdx=sub2ind(size(xMat),currenStartingSubtrackIdx,plotCurrentOnly*ones(length(currenStartingSubtrackIdx),1));
            scatter(xMat(coordIdx),yMat(coordIdx),'MarkerEdgeColor',c2,'SizeData',(72/3)^2)
        end
    end

end


% if movieInfo given then ALL detected features will be shown with jet
% colormap (not just the ones selected by the tracking)
if ~isempty(movieInfo)
    colorOverTime = jet(timeRange(2)-timeRange(1)+1);
    nFr=length(timeRange(1):timeRange(2));
    xCoord=arrayfun(@(i) vertcat(movieInfo(i).xCoord(:,1)),[timeRange(1):timeRange(2)]','UniformOutput',0);
    yCoord=arrayfun(@(i) vertcat(movieInfo(i).yCoord(:,1)),[timeRange(1):timeRange(2)]','UniformOutput',0);
    l=arrayfun(@(i) length(xCoord{i}),[1:nFr]');
    col=cell2mat(arrayfun(@(i) repmat(colorOverTime(i,:),[l(i) 1]), [1:nFr]','UniformOutput',0));
    xCoord=cell2mat(xCoord);
    yCoord=cell2mat(yCoord);

    outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
    xCoord(outOfRangeIdx) = [];
    yCoord(outOfRangeIdx) = [];
    col(outOfRangeIdx,:) = [];
    xCoord=xCoord-minX+1;
    yCoord=yCoord-minY+1;
    hold on
    scatter(xCoord,yCoord,'Marker','.','cData',col);

end

% get user selections and plot on image
if ask4sel
    if isempty(subIdx)
        subIdx=[1:length(trackType)]';
    end

    userEntry = 'yes';
    h=msgbox('Left-click to add points. Double-click to add final point and end selection, or press Enter to end selection without adding a final point. Press Backspace or Delete to remove previously selected point.','Help for track selection','help');
    uiwait(h);

    count = 1;
    while strcmpi(userEntry,'yes')

        %let the user choose the points of interest
        [x,y] = getpts;

        %find the time points of the indicated points
        for i=1:length(x)

            %find the distances between those points and the tracks
            distTrack2Point = sqrt((xMat-x(i)).^2+(yMat-y(i)).^2);

            %determine the minimum distance for each chosen point
            [rowChosen,colChosen] = find(distTrack2Point==min(distTrack2Point(:)));
            if ~isempty(rowChosen)
                subtrackChosen=rowChosen(1);
                trackChosen=trackData(subIdx(subtrackChosen),1);
                frameChosen=colChosen(1);

                % get whole track profile
                trackProfile=trackData(trackData(:,1)==trackChosen,:);
                selectedTracks{count,1} = trackProfile;

                disp(['Track: ' num2str(trackChosen) ',  Frame: ' num2str(frameChosen)]);
                disp(trackProfile)

                text(xMat(subtrackChosen,frameChosen),yMat(subtrackChosen,frameChosen),...
                    ['\leftarrow' num2str(trackProfile(1,1))],'Color','y','FontWeight','bold')
                count = count+1;
            else
                h=errordlg('No tracks found.');
                uiwait(h)
            end

        end

        %ask the user again whether to click on figure and get frame information
        userEntry = questdlg('Do you want to select more points?');

    end

end

%%%%% ~~ the end ~~ %%%%%

