function imarisApplication = imarisPlot3(plotData,aspectRatio)
%IMARISPLOT3 generates a 3D data plot
%
% SYNOPSIS imarisApplication = imarisPlot3(plotData)
%
% INPUT    plotData: Structure nSpotGroups-by-1 array containing the data
%                       for plotting with fields:
%               - XYZ: n-by-3 arrays of coordinates
%               - spotRadius (opt): size of plot spots. Can be set
%                    individually for every data point. Default: [0.25]
%                    If aspect ratio is set to other than 0, radius 1 is
%                    about 1/10 of the plot box. Otherwise, it is measured
%                    in the units of the input.
%               - color (opt): [R/G/B/Opacity] for the group of spots.
%                    Opacity of 0 = opaque. Default: [1,0,0,0]
%               - name (opt) : name of group of spots. Default: data_#,
%                    where # is the place of the set in the structure
%               - time (opt) : timepoint where the group should be plotted.
%                    Can be set individually for each spot (like
%                    spotRadius). Default: 1.
%               - class (opt): name (string!) of the class the group
%                    belongs to. Spot groups belonging to the same class
%                    will be collected in a folder in Imaris.
%                    If classes are used, ALL groups have to belong to a
%                    class.
%           aspectRatio (opt): Aspect ratio of the data. With the default,
%                              [1,1,1], the plot box has the shape of a
%                              cube. Specify [0,0,0] if a unit step should
%                              be the same in every direction.
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
def_spotRadius = 0.25;
def_color = [1,0,0,0]; % R/G/B/opacity
def_nameStub = 'data_'; % stub to which number will be added
def_aspectRatio = [1,1,1]; % aspect ratio
def_maxImSize = 1000; % maximum number of pixels for a single image

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
nCoords = zeros(nGroups,1);

for i = 1:nGroups
    if isempty(plotData(i).XYZ)
        error('every entry into plotData needs at least a nonempty field XYZ')
    end

    % test xyz
    plotData(i).XYZ = returnRightVector(plotData(i).XYZ,3);

    nCoords(i) = size(plotData(i).XYZ,1);

    % test spotSize
    if ~isfield(plotData,'spotRadius') || isempty(plotData(i).spotRadius)
        % fill in spotRadius
        plotData(i).spotRadius = ...
            repmat(def_spotRadius,nCoords(i),1);

        % make sure there are enough radii
    elseif length(plotData(i).spotRadius) == 1
        plotData(i).spotRadius = ...
            repmat(plotData(i).spotRadius,nCoords(i),1);

    else
        plotData(i).spotRadius = returnRightVector(plotData(i).spotRadius);
        if length(plotData(i).spotRadius) ~= nCoords(i)
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
    elseif ~ischar(plotData(i).name)
        error('name has to be a string')
    end

    % time
    if ~isfield(plotData,'time') || isempty(plotData(i).time)
        plotData(i).time = 1;
    end
    if length(plotData(i).time) == 1
        % adjust number of time
        plotData(i).time = repmat(plotData(i).time,[nCoords(i), 1]);
    else
        % make sure the time vector has the correct length
        plotData(i).time = returnRightVector(plotData(i).time);
        if length(plotData(i).time) ~= nCoords(i)
            error('if a list of timepoints is given, it has to match the number of coordinates!')
        end

    end

    % class
    if ~isfield(plotData,'class')
        plotData(i).class = '';

        % bad if there is something in class which is not a char
    elseif ~isempty(plotData(i).class) && ~ischar(plotData(i).class)
        error('class names have to be character arrays!')
    end
end

% collect classes
classNames = strvcat(plotData.class);
if isempty(classNames)
    nClasses = 0;
else
    % check whether there is a class assigned for every group
    if ~isequal(size(classNames,1),nGroups)
        error('every group needs to belong to a class if classes are assigned!')
    else
        % find unique classes
        [classList,dummy,classIdx] = unique(classNames,'rows');
        nClasses = size(classList,1);
    end
end


% aspect ratio
if nargin < 2 || isempty(aspectRatio)
    aspectRatio = def_aspectRatio;
end

%=========================================


%=========================================
% PREPARE DATA
%=========================================

% Here, we have to bring all the data into the right shape to be able to
% pass it to Imaris
% 1) Find extent of data, calculate limits
% 2) Calculate coordinate transformation
%       This has to be done, in case the aspect ratio is fixed. Otherwise,
%       we just adjust the image extent.
% 3) Later: join same colors, make 3D-histogram

% find extent. make sure we have enough space to place a sphere
allData = cat(1,plotData.XYZ);
allRadii = cat(1,plotData.spotRadius);
maxRadius = max(allRadii);
dataExtent = [min(allData) - repmat(maxRadius,[1,3]);...
    max(allData) + repmat(maxRadius,[1,3])];
dataRange = dataExtent(2,:) - dataExtent(1,:);

% extent in time
allTime = cat(1,plotData.time);
maxTime = max(allTime);


% imageDelta takes into account that the pixel coordinates are at the
% center of the pixel. Therefore, if we want to go from 0 to 10, we need 11
% pixels.
imageDelta = [1, 1, 1];

% select size of image - consider aspect ratio
if any(aspectRatio) == 0
    % use dataRange as image size

    % we do not need to make any kind of coordinate transformation
    doTransform = 0;

    % test for max image size
    bigExtent = ceil(dataRange);
    if prod(bigExtent) < def_maxImSize


        % assign imageSize - correct with imageDelta
        imSize = bigExtent + imageDelta;

        % find origin of image (coord - radius - some round-off delta)
        deltaExtent = bigExtent - dataRange;
        origin = dataExtent(1,:) - deltaExtent/2;

    else
        % we need to adjust image size. However, by setting the image
        % extent, we will get the correct aspect ratio.
        multImg = repmat((def_maxImSize/prod(bigExtent))^(1/3),[1,3]);

        % assign imSize
        imSize = ceil(bigExtent .* multImg) + imageDelta;
        % this is a reasonable approximation
        deltaExtent = (bigExtent - dataRange);
        % origin does not need to be adjusted - this is done below
        origin = dataExtent(1,:) - deltaExtent/2;

    end

    % set extends of image
    extendMin = (origin - ((imageDelta-1)/2));
    extendMax = (bigExtent + origin + ((imageDelta-1)/2));

    % store data for axis labels
    stepSize = (extendMax - extendMin)./(imSize-imageDelta); % one less steps that pix
    plotBoxData = [extendMin-stepSize/2; extendMax+stepSize/2; stepSize]';



else % use aspect ratio

    % we need to transform the ranges onto [0...10].*aspectRatio
    doTransform = 1;

    % to transform coordinates: subtract origin, divide by range to make data
    % going from 0 to 10
    origin = dataExtent(1,:);

    % "inner image": 10x10x10*aspectRatio
    imSize = ([10,10,10]).* aspectRatio + imageDelta;


    % we will divide the coordinates by divideRange -> multiply by imSize
    divideRange = dataRange ./ (imSize - imageDelta);

    % the extends of the image will be imSize - 1;
    extendMin = zeros(1,3);
    extendMax = imSize - 1;

    % store data for axis labels
    stepSize = divideRange;
    plotBoxData = [origin - divideRange/2;...
        origin + dataExtent(2,:) + divideRange/2; stepSize]';


end

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
imaDataSet.SetData(single((zeros([imSize,1,maxTime]))));
% image starts at imageDelta/2 pixels left of origin - divide to ensure
% aspect ratio!
imaDataSet.mExtendMinX =  extendMin(1);
imaDataSet.mExtendMinY =  extendMin(2);
imaDataSet.mExtendMinZ =  extendMin(3);
% image ends at imageDelta/2 pixels right of dataRange + origin - divide to
% ensure aspect ratio!

imaDataSet.mExtendMaxX = extendMax(1);
imaDataSet.mExtendMaxY = extendMax(2);
imaDataSet.mExtendMaxZ = extendMax(3);
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

switch nClasses

    case 0

        % loop through data and add spots
        for iGroup = 1:nGroups

            % create spot object
            imaSpots = imaApp.mFactory.CreateSpots;

            % set spots

            if doTransform
                coords = single((plotData(iGroup).XYZ - repmat(origin,[nCoords(iGroup),1]) )...
                    ./ repmat(divideRange,[nCoords(iGroup),1]));
            else
                coords = single(plotData(iGroup).XYZ);
            end
            imaSpots.Set(coords, single(plotData(iGroup).time-1),...
                single(plotData(iGroup).spotRadius));


            % set color
            color = single(plotData(iGroup).color);
            imaSpots.SetColor(color(1),color(2),color(3),color(4));

            % set name
            imaSpots.mName = plotData(iGroup).name;

            imaSSSpots.AddChild(imaSpots);

        end

    otherwise

        % loop through classes, and add spots to each class

        for iClass = 1:nClasses

            % make new subdir
            imaClass = imaApp.mFactory.CreateDataContainer;
            % name it
            imaClass.mName = classList(iClass,:);
            % and add to spot folder
            imaSSSpots.AddChild(imaClass);

            % find spot groups belonging to class
            % groupIdx points to entry in plotData
            groupIdx = find(classIdx == iClass);

            for iGroup = groupIdx'

                % create spot object
                imaSpots = imaApp.mFactory.CreateSpots;

                % set spots

                if doTransform
                    coords = single((plotData(iGroup).XYZ - repmat(origin,[nCoords(iGroup),1]) )...
                        ./ repmat(divideRange,[nCoords(iGroup),1]));
                else
                    coords = single(plotData(iGroup).XYZ);
                end
                imaSpots.Set(coords, single(plotData(iGroup).time-1),...
                    single(plotData(iGroup).spotRadius));


                % set color
                color = single(plotData(iGroup).color);
                imaSpots.SetColor(color(1),color(2),color(3),color(4));

                % set name
                imaSpots.mName = plotData(iGroup).name;

                imaClass.AddChild(imaSpots);
            end % for iGroup = groupdx'
        end % for iClass = 1:nClasses

end % switch

%===========================


%===========================
% SHOW DATA ABOUT AXES
%===========================

% since it is not possible in Imaris to set axes labels, we give the
% information in matlab

% the center of the pixel is the actual coordinate. However, Imaris plots
% the frame around the pixel. Therefore, we need to add/subtract the
% pixelsize
% We do this already when calculating plotBoxData
xString = sprintf('Xmin: %f, Xmax: %f, Xstep: %f',...
    plotBoxData(1,:));
yString = sprintf('Ymin: %f, Ymax: %f, Ystep: %f',...
    plotBoxData(2,:));
zString = sprintf('Zmin: %f, Zmax: %f, Zstep: %f',...
    plotBoxData(3,:));
%plot to command line
disp(sprintf(['axis limits and step sizes:\n' xString '\n' yString '\n' zString]))

%=========================

%=========================
% ASSIGN OUTPUT
%=========================

if nargout > 0
    imarisApplication = imaApp;
end