function imarisApplication = imarisPlot3(plotData)
%IMARISPLOT3 generates a 3D data plot
%
% SYNOPSIS imarisApplication = imarisPlot3(plotData)
%
% INPUT    plotData: Structure nSpotGroups-by-1 array containing the data
%                       for plotting with fields:
%               - XYZ: n-by-3 arrays of coordinates
%               - spotRadius (opt): size of plot spots. Can be set
%                    individually for every data point. Default: [0.5]
%               - color (opt): [R/G/B/Opacity] for the group of spots.
%                    Opacity of 0 = opaque. Default: [1,0,0,0.5]
%               - name (opt) : name of group of spots. Default: data_#,
%                    where # is the place of the set in the structure
%
% OUTPUT  imarisApplication: Handle to the imaris Application
%
% REMARKS: (1) Imaris will close once all its handles are deleted. Even if
%              no output argument is assigned, a variable called
%              imarisApplication will be created in the workspace
%          (2) Imaris can not handle doubles. All values will be converted
%              to singles
%
% c: jonas, 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%======================
% TEST INPUT
%======================

% defaults
def_spotRadius = 0.5;
def_color = [1,0,0,0.5]; % R/G/B/opacity
def_nameStub = 'data_'; % stub to which number will be added

% input arguments
if nargin == 0 || isempty(plotData) || ~isstruct(plotData)
    error('please specify a plotData-Structure as input arguments')
end

% check for xyz
if ~isfield(plotData,'XYZ')
    error('please specify at least a coordinate field!')
end

% loop through structure and make sure everything's allright
nGroups = length(plotData);

for i = 1:nGroups
    if isempty(plotData(i).XYZ)
        error('every entry into plotData needs at least a nonempty field XYZ')
    end

    % test xyz
    plotData(i).XYZ = returnRightVector(plotData(i).XYZ,3);

    % test spotSize
    if ~isfield(plotData,'spotRadius') || isempty(plotData(i).spotRadius)
        % fill in spotRadius
        plotData(i).spotRadius = ...
            repmat(def_spotRadius,size(plotData(i).XYZ,1),1);

        % make sure there are enough radii
    elseif length(plotData(i).spotRadius) == 1
        plotData(i).spotRadius = ...
            repmat(plotData(i).spotRadius,size(plotData(i).XYZ,1),1);

    else
        plotData(i).spotRadius = returnRightVector(plotData(i).spotRadius);
        if length(plotData(i).spotRadius) ~= size(plotData(i).XYZ,1)
            error('if a list of radii is given, it has to match the number of coordinates!')
        end

    end

    % color
    if ~isfield(plotData,'color') || isempty(plotData(i).color)
        plotData(i).color = def_color;
    else
        plotData(i).color = returnRightVector(plotData(i).color,1,'r');
        if length(plotData(i).color) ~= 4
            error('please specify color as [R,G,B,opacity]')
        end
    end

    % name
    if ~isfield(plotData,'name') || isempty(plotData(i).name)
        plotData(i).name = [def_nameStub,num2str(i)];
    elseif ~isstr(plotData(i).name)
        error('name has to be a string')
    end
end



%=========================================


%=========================================
% PREPARE DATA
%=========================================

% Here, we have to bring all the data into the right shape to be able to
% pass it to Imaris
% 1) Find extent of data, calculate limits
% 2) Calculate coordinate transformation
% 3) Later: join same colors, make 3D-histogram

% find extent
allData = cat(1,plotData.XYZ);
dataExtent = [min(allData);max(allData)];
dataRange = dataExtent(2,:) - dataExtent(1,:);


% to transform coordinates: subtract origin, divide by range to make data
% going from 0 to 10
origin = dataExtent(1,:);
divideRange = dataRange/10;

%========================================


%========================================
% PUT SPOTS INTO IMARIS
%========================================

% 1) start Imaris, set 10-by-10-by-10 matrix (later used for histo)
% 2) loop through plotData, put points into Surpass

% start new imaris
imaApp = imarisStartNew;
% insert a pause. 0.2 seconds works, too, but I want to avoid a crash on
% slower machines.
pause(0.5) 




% make dataSet and put into imaris
imaDataSet = imaApp.mFactory.CreateDataSet;
imaDataSet.SetData(single(zeros(12,12,12)));
imaDataSet.mExtendMinX = -1;
imaDataSet.mExtendMinY = -1;
imaDataSet.mExtendMinZ = -1;
imaApp.mDataSet = imaDataSet;

% make top-level surpass scene 
imaSurpassScene = imaApp.mFactory.CreateDataContainer;

% fill surpass scene with light and frame
imaLight = imaApp.mFactory.CreateLightSource;
imaSurpassScene.AddChild(imaLight);
imaFrame = imaApp.mFactory.CreateFrame;
imaSurpassScene.AddChild(imaFrame);


% add surpass scene and set view
imaApp.mSurpassScene = imaSurpassScene;
imaApp.mViewer = 'eViewerSurpass';

% add surpass scene for spots
imaSSSpots = imaApp.mFactory.CreateDataContainer;
imaSSSpots.mName = 'Spots';
imaSurpassScene.AddChild(imaSSSpots);

% loop through data and add spots
for iGroup = 1:nGroups
    
    % create spot object
    imaSpots = imaApp.mFactory.CreateSpots;
    
    % set spots
    nCoords = size(plotData(iGroup).XYZ,1);
    coords = single((plotData(iGroup).XYZ - repmat(origin,nCoords,1))...
        ./repmat(divideRange,nCoords,1));
    imaSpots.Set(coords, single(zeros(nCoords,1)),...
        single(plotData(iGroup).spotRadius));
    
    % set color
    color = single(plotData(iGroup).color);
    imaSpots.SetColor(color(1),color(2),color(3),color(4));
    
    % set name
    imaSpots.mName = plotData(iGroup).name;
    
    imaSSSpots.AddChild(imaSpots);
    
end

%===========================


%===========================
% SHOW DATA ABOUT AXES
%===========================

% since it is not possible in Imaris to set axes labels, we give the
% information in matlab
xStep = dataRange(1)/10;
yStep = dataRange(2)/10;
zStep = dataRange(3)/10;
xString = sprintf('Xmin: %f, Xmax: %f, Xstep: %f',...
    dataExtent(1,1)-1.5*xStep,...
    dataExtent(2,1)+1.5*xStep,...
    xStep);
yString = sprintf('Ymin: %f, Ymax: %f, Ystep: %f',...
    dataExtent(1,2)-1.5*yStep,...
    dataExtent(2,2)+1.5*yStep,...
    yStep);
zString = sprintf('Zmin: %f, Zmax: %f, Zstep: %f',...
    dataExtent(1,3)-1.5*zStep,...
    dataExtent(2,3)+1.5*zStep,...
    zStep);
%plot to command line
sprintf([xString '\n' yString '\n' zString])