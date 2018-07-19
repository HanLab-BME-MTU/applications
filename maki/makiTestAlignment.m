function makiTestAlignment(dataStruct,what2plot,whichRotation)
%TESTROTATION ...
%
% SYNOPSIS: testRotation(dataStruct,what2plot,whichRotation)
%
% INPUT dataStruct : maki data structure (see makiMakeDataStruct for
%                    details).
%       what2plot  : 'p' to plot positions, 'd' to plot displacements.
%       whichRotation: 1 to plot rotatedCoord from planeFit, 2 to plot
%                      alignedCoord from frameAlignment.
%
% created by: kjaqaman
% DATE: 31-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 

% check input for empty
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

nTimePoints = dataStruct.dataProperties.movieSize(end);
rotationInfo.initCoord = dataStruct.initCoord;
rotationInfo.tracks = dataStruct.tracks;
rotationInfo.projectName = dataStruct.projectName;
rotationInfo.what2plot = what2plot;
rotationInfo.whichRotation = whichRotation;
switch whichRotation
    case 1
        rotationInfo.planeFit = dataStruct.planeFit;
    case 2
        rotationInfo.frameAlignment = dataStruct.frameAlignment;
end

% initialize plot window
figH = figure('Name',sprintf('Frame alignment %s (%2d / %2d)',rotationInfo.projectName,1, nTimePoints));
set(figH,'KeyPressFcn',@figure_keyPress);
% user data of the figure stores the current time point and the number of frames
rotationInfo.timeInfo = struct('nTimePoints',nTimePoints,...
    'currentTimePoint',1);
set(figH,'UserData',rotationInfo);

if strcmp(what2plot,'p')
    plotPoints(figH,1);
elseif strcmp(what2plot,'d')
    plotDisp(figH,1);
end


%% LOCAL FUNCTIONS

%% figure_keyPress

function figure_keyPress(src,event)

rotationInfo = get(src,'UserData');
timeInfo = rotationInfo.timeInfo;
what2plot = rotationInfo.what2plot;

switch event.Key
    case 'uparrow'
        timeInfo.currentTimePoint = timeInfo.currentTimePoint + 1;
    case 'downarrow'
        timeInfo.currentTimePoint = timeInfo.currentTimePoint - 1;
end

if timeInfo.currentTimePoint > timeInfo.nTimePoints
    timeInfo.currentTimePoint = 1;
end
if timeInfo.currentTimePoint < 1
    timeInfo.currentTimePoint = timeInfo.nTimePoints;
end

rotationInfo.timeInfo = timeInfo;

set(src,'UserData',rotationInfo);

if strcmp(what2plot,'p')
    plotPoints(src,0);
elseif strcmp(what2plot,'d')
    plotDisp(src,1);
end
    
% end of function figure_keyPress

%% plotPoints

function plotPoints(figH,initialPlot)

cameraPropertyNames = {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'};

% figure is still open
    figure(figH)
    cameraProps = get(gca,cameraPropertyNames);
    
    
    rotationInfo = get(figH,'UserData');
    timeInfo = rotationInfo.timeInfo;
    currentTimePoint = timeInfo.currentTimePoint;
    whichRotation = rotationInfo.whichRotation;
    originalCoord = rotationInfo.initCoord(currentTimePoint).allCoord(:,1:3);
    switch whichRotation
        case 1
            alignedCoord = rotationInfo.planeFit(currentTimePoint).rotatedCoord(:,1:3);
        case 2
            alignedCoord = rotationInfo.frameAlignment(currentTimePoint).alignedCoord(:,1:3);
    end            
    
%     centerMassOriginal = mean(originalCoord);
%     centerMassAligned = mean(alignedCoord);
%     alignedCoord = alignedCoord + repmat(centerMassOriginal-centerMassAligned,size(alignedCoord,1),1);
    
    plot3(originalCoord(:,1),originalCoord(:,2),originalCoord(:,3),'o')
    hold on;
    plot3(alignedCoord(:,1),alignedCoord(:,2),alignedCoord(:,3),'r+')
    axis([-10 30 -10 30 -10 30]);
    grid on;
    
    if ~initialPlot
        set(gca,cameraPropertyNames,cameraProps);
    end
    
    set(figH,'Name',sprintf('Frame alignment %s (%2d / %2d)',rotationInfo.projectName,...
        timeInfo.currentTimePoint,timeInfo.nTimePoints));
    hold off;

% end of function plotPlanes

%% plotDisp

function plotDisp(figH,initialPlot)

cameraPropertyNames = {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'};

% figure is still open
    figure(figH)
    cameraProps = get(gca,cameraPropertyNames);

    rotationInfo = get(figH,'UserData');
    timeInfo = rotationInfo.timeInfo;
    whichRotation = rotationInfo.whichRotation;

    %put tracks into matrix format
    [tracksInfo,tracksIndx] = convStruct2MatNoMS(rotationInfo.tracks);
    numTracks = size(tracksInfo,1);

    currentTimePoint = timeInfo.currentTimePoint;
    if currentTimePoint == timeInfo.nTimePoints
        nextTimePoint = 1;
    else
        nextTimePoint = currentTimePoint + 1;
    end
    
    originalCoord = tracksInfo(:,[((currentTimePoint-1)*8+1:(currentTimePoint-1)*8+3) ...
        ((nextTimePoint-1)*8+1:(nextTimePoint-1)*8+3)]);
    
    goodPoints = find(~isnan(sum(originalCoord,2)))';
    
    alignedCoord = originalCoord;

    switch whichRotation
        case 1
            for iTrack = 1 : numTracks
                iFeature = tracksIndx(iTrack,currentTimePoint);
                if iFeature ~= 0
                    alignedCoord(iTrack,1:3) = rotationInfo.planeFit(currentTimePoint).rotatedCoord(iFeature,1:3);
                end
                iFeature = tracksIndx(iTrack,nextTimePoint);
                if iFeature ~= 0
                    alignedCoord(iTrack,4:6) = rotationInfo.planeFit(nextTimePoint).rotatedCoord(iFeature,1:3);
                end
            end
        case 2
            for iTrack = 1 : numTracks
                iFeature = tracksIndx(iTrack,currentTimePoint);
                if iFeature ~= 0
                    alignedCoord(iTrack,1:3) = rotationInfo.frameAlignment(currentTimePoint).alignedCoord(iFeature,1:3);
                end
                iFeature = tracksIndx(iTrack,nextTimePoint);
                if iFeature ~= 0
                    alignedCoord(iTrack,4:6) = rotationInfo.frameAlignment(nextTimePoint).alignedCoord(iFeature,1:3);
                end
            end
    end

    for iTrack = goodPoints
        plot3(originalCoord(iTrack,[1 4])',originalCoord(iTrack,[2 5]'),originalCoord(iTrack,[3 6]),'LineWidth',3)
        hold on;
    end
    for iTrack = goodPoints
        plot3(alignedCoord(iTrack,[1 4])',alignedCoord(iTrack,[2 5]'),alignedCoord(iTrack,[3 6]),'r','LineWidth',3)
        hold on;
    end
    axis([-10 30 -10 30 -10 30]);
    grid on;
    
    if ~initialPlot
        set(gca,cameraPropertyNames,cameraProps);
    end
    
    set(figH,'Name',sprintf('Frame alignment %s (%2d / %2d)',rotationInfo.projectName,...
        timeInfo.currentTimePoint,timeInfo.nTimePoints));
    hold off;

% end of function plotDisp
