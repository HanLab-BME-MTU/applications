function [windows]  = calculateWindowInformation(data,mask);
% calculateWindowInformation divides a cell mask into windows and
% calculates pairCorrelation, nucleationDensity, and pit lifetimes for
% each window
%
%inputs:
%       data
%
%       mask (optional): Give mask to use. Otherwise program will have you
%       hand select one and then cleanup using detections
%
%outputs:
%       windows: structure with the pairCorrelation, nucleationDensity,
%       list of pit nucleation sites, and lifetimes for a given window

%LOAD MOVIES
if nargin < 1 || isempty(data)
    data = loadConditionData;
end

%MAKE CELL MASK

%BY HAND
if nargin < 2 || isempty(mask)
    [filename, pathname] = uigetfile('.tif','choose image to segment by hand');
    image = imread([pathname filesep filename]);
    imshow(image,[]);
    %make mask
    maskHandCut = roipoly;
    
    %REFINE USING DETECTIONS
    %area mask parameters
    closureRadius = 50;
    dilationRadius = 5;
    doFill = 1;
    %Load Lifetime Information
    lftInfo = load([data.source filesep 'LifetimeInfo' filesep 'lftInfo']);
    lftInfo = lftInfo.lftInfo;
    matX = lftInfo.Mat_xcoord;
    % y-coordinate matrix
    matY = lftInfo.Mat_ycoord;
    % image size
    imsize  = data.imagesize;
    %MAKE MASK
    imsizS = [imsize(2) imsize(1)];
    maskDetections = makeCellMaskDetections([matX(:),matY(:)],closureRadius,dilationRadius,doFill,imsize,1,[]);
    mask = zeros(size(maskDetections));
    mask = maskDetections & maskHandCut;
end

%MAKE WINDOWS
goodWindows = 0;
perpSize = 200;
paraSize = 500;
while ~goodWindows
    %make polygons
    windowPolygons = getMaskWindows(mask,perpSize,paraSize);
    %visualize polygons
    plotWindows(windowPolygons)
    goodWindows = input('Is the window size good (1 for yes 0 for no)?');
    if ~goodWindows
        perpSize = input('Enter window dimension perpendicular to membrane');
        paraSize = input('Enter window dimension parallel to the membrane');
    end
end

%CALCULATE PAIR CORRELATION, LIFETIMES, AND NUCLEATION DENSITY FOR EACH WINDOW
windows = plotPairCorrelation_SpatialWindows(data,[],[],[],windowPolygons);
