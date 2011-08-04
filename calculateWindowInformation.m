function [windows,windowPolygons]  = calculateWindowInformation(data,mask);
% calculateWindowInformation divides a cell mask into windows and
% calculates pairCorrelation, nucleationDensity, and pit lifetimes for
% each window
%
%inputs:
%       data: The typical structure, but only one movie!!!
%
%       mask (optional): Give mask to use. Otherwise program will have you
%       hand select one and then cleanup using detections
%
%outputs:
%       windows: structure with the pairCorrelation, nucleationDensity,
%       list of pit nucleation sites, and lifetimes for a given window


%MAKE CELL MASK

%BY HAND
if nargin < 2 || isempty(mask)
    [filename, pathname] = uigetfile('.tif','choose image to segment by hand');
    image = imread([pathname filesep filename]);
    handle = figure;
    imshow(image,[]);
    %make mask
    maskHandCut = roipoly;
    close(handle)
    
    %REFINE USING DETECTIONS
    %area mask parameters
    closureRadius = 50;
    dilationRadius = 5;
    doFill = 1;
    imsize  = data.imagesize;
    %Load Lifetime Information
    try load([data.source filesep 'Tracking' filesep 'trackAnalysis.mat'])
        
        %positions used to calculate mask
        maskPositionsX = arrayfun(@(t) t.x(1),tracks)';
        maskPositionsY = arrayfun(@(t) t.y(1),tracks)';
        maskPositions = [maskPositionsX, maskPositionsY];
        
    catch ME
        
        lftInfo = load([data.source filesep 'LifetimeInfo' filesep 'lftInfo']);
        lftInfo = lftInfo.lftInfo;
        matX = lftInfo.Mat_xcoord;
        % y-coordinate matrix
        matY = lftInfo.Mat_ycoord;
        %positions used to calculate mask
        maskPositions = [matX(~isnan(matX)),matY(~isnan(matY))];
        
    end
    
    %MAKE MASK
    imsizS = [imsize(2) imsize(1)];
    maskDetections = makeCellMaskDetections(maskPositions,closureRadius,dilationRadius,doFill,imsize,0,[]);
    mask = zeros(size(maskDetections));
    mask = maskDetections & maskHandCut;
end

%MAKE WINDOWS
goodWindows = 0;
perpSize = 100;
paraSize = 100;
while ~goodWindows
    %make polygons
    windowPolygons = getMaskWindows(mask,perpSize,paraSize);
    %visualize polygons
    handle = figure;
    plotWindows(windowPolygons)
    goodWindows = input('Is the window size good (1 for yes 0 for no)?');
    if ~goodWindows
        perpSize = input('Enter window dimension perpendicular to membrane');
        paraSize = input('Enter window dimension parallel to the membrane');
    end
    close(handle)
end

%CALCULATE PAIR CORRELATION, LIFETIMES, AND NUCLEATION DENSITY FOR EACH WINDOW
windows = plotPairCorrelation_SpatialWindows(data,[],[],[],windowPolygons);
