function plotTracks2D_EB3(trackedFeatureInfo,timeRange,colorTime,markerType,...
    indicateSE,newFigure,image,flipXY,ask4sel,velocityMatrix)
%PLOTTRACKS2D plots a group of tracks in 2D and allows user to click on them and extract track information
%
%SYNOPSIS plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure,image,flipXY,ask4sel)
%
%INPUT  trackedFeatureInfo: -- EITHER --
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row
%                           consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%                           -- OR --
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       timeRange         : 2-element row vector indicating time range to plot.
%                           Optional. Default: whole movie.
%       colorTime         : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                           blue, red, etc.
%                           -'2' if tracks are colored by cycling through
%                           the plot's default color order (except black)
%                           Optional. Default: 'r'.
%                           -'vel' if tracks should be colored by speed. In
%                           this case make sure to add the optional
%                           parameter "trackInfo"
%       markerType        : String indicating marker type for plotting.
%                           Only used if colorTime is not '1'.
%                           Optional. Default: 'none'.
%       indicateSE        : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%       flipXY            : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%       ask4sel           : 1 if user should be asked to select tracks in
%                           plot in order to show track information.
%                           Optional. Default: 1.
%       velocityMatrix    : nTrack x nFrames-1 matrix containing the
%                           velocities between timepoints (e.g. estimated
%                           by getVelocitiesFromMat (opt, used if colorTime
%                           = 'vel')
%
%OUTPUT The plot.
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get number of tracks and number of time points

[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;


errFlag = 0;

%check whether a time range for plotting was input
if nargin < 2 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracks2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether colorTime was input
if nargin < 3 || isempty(colorTime)
    colorTime = 'r';
end
% make sure colorTime 1,2 are strings
if isnumeric(colorTime) && isscalar(colorTime)
    colorTime = num2str(colorTime);
end
if strmatch(colorTime,'vel')
    if nargin<10
        error('if plotting by velocity, input must include trackInfo')
    end
end

%check whether markerType was input
if nargin < 4 || isempty(markerType)
    markerType = 'none';
end

%check whether indicateSE was input
if nargin < 5 || isempty(indicateSE)
    indicateSE = 1;
else
    if indicateSE ~= 0 && indicateSE ~= 1
        disp('plotTracks2D: indicateSE should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 6 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        disp('--plotTracks2D: newFigure should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether user supplied an image
if nargin < 7 || isempty(image)
    image = [];
end

if nargin < 8 || isempty(flipXY)
    flipXY = false;
end
if nargin < 9 || isempty(ask4sel)
    ask4sel = true;
end

%exit if there are problem in input variables
if errFlag
    disp('--plotTracks2D: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% tracks should be in matrix format



%indicate that each track consists of one segment
numSegments = ones(numTracks,1);

%locate the row of the first track of each compound track in the
%big matrix of all tracks
%in this case of course every compound track is simply one track
%without branches
trackStartRow = (1:numTracks)';


%get the x,y-coordinates of features in all tracks
if flipXY
    tracksY = trackedFeatureInfo(:,1:8:end)';
    tracksX = trackedFeatureInfo(:,2:8:end)';
else
    tracksX = trackedFeatureInfo(:,1:8:end)';
    tracksY = trackedFeatureInfo(:,2:8:end)';
end

%find x-coordinate limits
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
maxYCoord =  ceil(max(tracksY(:)));

%calculate the number of time points to be plotted
numTimePlot = timeRange(2) - timeRange(1) + 1;

%define colors to loop through in case colorTime = '2'
colorLoop = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: r,g,b,y,m,c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%if the user wants to plot in a new figure window

if newFigure

    %open new figure window
    figure

    if ~isempty(image) & colorTime~='vel'
        imshow(image,[]); %plot the image
    elseif ~isempty(image)%if user did not supply an image
        img = repmat(0.75*ones(size(image)),[1,1,3]);
        imshow(img); %plot an empty image    
    else
        img = repmat(0.75*ones(maxYCoord,maxXCoord),[1,1,3]);
        imshow(img); %plot an empty image    

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

end

%hold on figure
hold on

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

switch colorTime

    case '1' %if user wants to color-code time

        %plot tracks ignoring missing points
        %gaps are depicted as a dotted black line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
        end

        %get the overall color per time interval
        colorOverTime = timeColormap(numTimePlot);

        %overlay tracks with color coding wherever a feature has been detected
        for i=1:numTimePlot-1
            plot(tracksXP(i:i+1,:),tracksYP(i:i+1,:),'color',colorOverTime(i,:));
        end

    case '2' %no time color-coding, loop through series of colors to color tracks

        %plot tracks by looping through colors
        %missing intervals are indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
            plot(tracksXP(:,i),tracksYP(:,i),'color',colorLoop(mod(i-1,6)+1,:),...
                'marker',markerType);
        end

    case 'vel' % color code according to velocity
        cMap = jet;
        %plot tracks ignoring missing points
        %gaps are depicted as a dotted white line
        for i = 1:numTracks
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'w:','linewidth',1);
        end
        
        maxVel = max(abs(velocityMatrix(:)));
        scaleVec = linspace(-maxVel,maxVel,size(jet,1))';
        D = abs(createDistanceMatrix(velocityMatrix(:),scaleVec));
        [B,idx] = sort(D,2); % idx is index of colorMap
        
        % each value of velocityMatrix now has corresponding colorMap
        % index, cIdx
        cIdx = reshape(idx(:,1),size(velocityMatrix));
        % cut down to frame range we care about
        cIdx = cIdx(:,timeRange(1):timeRange(2)-1);
        
        for iTrack = 1:numTracks
            for iTimePoint=1:numTimePlot-1
            plot(tracksXP(iTimePoint:iTimePoint+1,iTrack),tracksYP(iTimePoint:iTimePoint+1,iTrack),'color',cMap(cIdx(iTrack,iTimePoint),:),...
                'marker',markerType,'linewidth',1);
           end
        end
        colormap(cMap);
        colorbar('YTick',round(linspace(1,size(cMap,1),7)),'YTickLabel',round(1000.*linspace(-maxVel,maxVel,7))./1000)
      

    otherwise %no time color-coding, all tracks same color

        %plot tracks with the line color indicated
        %missing intervals are indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
            plot(tracksXP(:,i),tracksYP(:,i),colorTime,'marker',markerType);
        end

end %(switch colorTime)


if indicateSE %if user wants to indicate starts and ends

    %find the beginning and end of each track
    for i=numTracks:-1:1
        timePoint = find(~isnan(tracksX(:,i)));
        startInfo(i,:) = [tracksX(timePoint(1),i) ...
            tracksY(timePoint(1),i) timePoint(1)];
        endInfo(i,:) = [tracksX(timePoint(end),i) ...
            tracksY(timePoint(end),i) timePoint(end)];
    end

    %place circles at track starts and squares at track ends if they happen to
    %be in the plotting region of interest
    switch colorTime
        case '1'
            indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
            plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o');
            indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
            plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square');
        case '2'
            indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
            plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o','MarkerEdgeColor','w');
            indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
            plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square','MarkerEdgeColor','w');
        case 'vel'
            indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
            plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o','MarkerEdgeColor','w');
            indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
            plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square','MarkerEdgeColor','w');
        
        otherwise
            indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
            plot(startInfo(indx,1),startInfo(indx,2),colorTime,...
                'LineStyle','none','marker','o');
            indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
            plot(endInfo(indx,1),endInfo(indx,2),colorTime,...
                'LineStyle','none','marker','square');
    end


end %(if indicateSE)

if ask4sel

    %ask the user whether to click on figure and get frame information
    userEntry = input('select points in figure? y/n ','s');

    while strcmp(userEntry,'y')

        %let the user choose the points of interest
        [x,y] = getpts;

        %find the time points of the indicated points
        for i=1:length(x)

            %find the distances between those points and the tracks
            distTrack2Point = (tracksXP-x(i)).^2+(tracksYP-y(i)).^2;

            %determine the minimum distance for each chosen point
            [frameChosen,rowChosen] = find(distTrack2Point==min(distTrack2Point(:)));

            %go over all chosen rows
            for j = 1 : length(rowChosen)

                %find the track corresponding to each minimum distance
                trackChosen = find(trackStartRow <= rowChosen(j),1,'last');
                segmentChosen = rowChosen(j) - trackStartRow(trackChosen) + 1;

                disp(['Track: ' num2str(trackChosen) ...
                    '   Segment: ' num2str(segmentChosen) ...
                    '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
                    '   Coordinates: ' num2str(tracksXP(frameChosen(j),rowChosen(j))) ...
                    ' ' num2str(tracksYP(frameChosen(j),rowChosen(j)))  ]);

            end

        end

        %ask the user again whether to click on figure and get frame information
        userEntry = input('select points again? y/n ','s');

    end

end

%%%%% ~~ the end ~~ %%%%%

