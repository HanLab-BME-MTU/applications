function makiPlotRotatingPlanes(dataStruct)
%MAKIPLOTROTATINGPLANES plots the plane rotation calculated by makiFitPlane
%
% SYNOPSIS: makiPlotRotatingPlanes(dataStruct)
%
% INPUT dataStruct : maki data structure (see makiMakeDataStruct for details)
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: gdanuser
% DATE: 26-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 

% check input for empty
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
% check whether analysis has been done
if isempty(dataStruct.planeFit)
    dataStruct = makiFitPlane(dataStruct,0);
end

nTimePoints = dataStruct.dataProperties.movieSize(end);
planeInfo.planeFit = dataStruct.planeFit;
planeInfo.projectName = dataStruct.projectName;

% check whether any of the plane is defined
defined = 0;
timePoint = 1;
while (timePoint <= nTimePoints) && ~defined
    defined = ~isempty(dataStruct.planeFit(timePoint).planeVectors);
    timePoint = timePoint + 1;
end
if ~defined
    warning('None of the time points contains a valid planeVector field');
    return;
end

% initialize plot window
figH = figure('Name',sprintf('Plane Orientation %s (%2d / %2d)',dataStruct.projectName,1, nTimePoints));
set(figH,'KeyPressFcn',@figure_keyPress);
% user data of the figure stores the current time point and the number of frames
planeInfo.timeInfo = struct('nTimePoints',nTimePoints,...
    'currentTimePoint',1);
set(figH,'UserData',planeInfo);
plotPlanes(figH,1);


%% LOCAL FUNCTIONS

function figure_keyPress(src,event)

planeInfo = get(src,'UserData');
timeInfo = planeInfo.timeInfo;

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

planeInfo.timeInfo = timeInfo;

set(src,'UserData',planeInfo);
plotPlanes(src,0,event.Key);
    
% end of function figure_keyPress

function plotPlanes(figH,initialPlot,key)

cameraPropertyNames = {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'};

if nargin == 3
    fwd = strmatch(key,'uparrow');
else
    fwd = 1;
end

% figure is still open
figure(figH)
cameraProps = get(gca,cameraPropertyNames);


planeInfo = get(figH,'UserData');
timeInfo = planeInfo.timeInfo;
currentTimePoint = timeInfo.currentTimePoint;
planeData = planeInfo.planeFit(currentTimePoint);

if ~isempty(planeData.planeVectors)
    % plot normal vector in red
    [dx,dy,dz] = eigenVec2plotVec(planeData.planeVectors(:,1));
    if planeData.planeVectorClassifier
        % the plane is defined by the eigenvectors
        plot3(dx,dy,dz,'r-');
        hold on
    
        % which column in the eigenvector matrix defined the normal ?
        normalIndx = find(planeData.planeVectors(1,1)==planeData.eigenVectors(1,:));
        eigenVecRatio = planeData.eigenValues(normalIndx)/...
            mean(planeData.eigenValues(setdiff([1 2 3],normalIndx)));
        % replot normal vector multiplied by the eigen vector ratio in much
        % larger width
        lineH = plot3(eigenVecRatio*dx,eigenVecRatio*dy,eigenVecRatio*dz,'r-');
        set(lineH,'LineWidth',5);
    else 
        plot3(dx,dy,dz,'r--');
        hold on
    end
        
    % plot in-plane vectors
    for i =2:3
        [dx,dy,dz] = eigenVec2plotVec(planeData.planeVectors(:,i));
        if planeData.planeVectorClassifier
            plot3(dx,dy,dz,'k-');
        else
            plot3(dx,dy,dz,'k--');
        end
    end
    
    if ~initialPlot
        set(gca,cameraPropertyNames,cameraProps);
    end
    axis([-1 1 -1 1 -1 1]);
    
    set(figH,'Name',sprintf('Plane Orientation %s (%2d / %2d)',planeInfo.projectName,...
        timeInfo.currentTimePoint,timeInfo.nTimePoints));
    hold off;
else
    % increase or deacreas the currrent time point by 1 and write it back into the
    % figure handle's user data
    if fwd
        timeInfo.currentTimePoint = timeInfo.currentTimePoint + 1;
    else 
        timeInfo.currentTimePoint = timeInfo.currentTimePoint - 1;
    end
    
    if timeInfo.currentTimePoint > timeInfo.nTimePoints
        timeInfo.currentTimePoint = 1;
    end
    if timeInfo.currentTimePoint < 1
        timeInfo.currentTimePoint = timeInfo.nTimePoints;
    end
    
    planeInfo.timeInfo = timeInfo;
    set(figH,'UserData', planeInfo);
    % recursive call of plot plane; the recursion will end as soon as
    % one time point with a valid plane is found
    plotPlanes(figH,0,key);
end

% end of function plotPlanes
    
function [dx,dy,dz]=eigenVec2plotVec(vec)
% service function of makiPlotRotatingPlanes to rearrange the eigenvector
% into a format for line plotting

dx = [0 vec(1)];
dy = [0 vec(2)];
dz = [0 vec(3)];

% end of function eigenVec2plotVec