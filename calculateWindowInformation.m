function [windows,windowPolygons]  = calculateWindowInformation(data,varargin);
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

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(var) isstruct(var)&numel(var) == 1);
ip.addParamValue('mask', [], @isnumeric);
ip.addParamValue('chooseImageToSegment', false, @islogical);
ip.parse(data, varargin{:});
%MAKE CELL MASK

%BY HAND
if isempty(ip.Results.mask)
    if ip.Results.chooseImageToSegment
        [filename, pathname] = uigetfile('.tif','choose image to segment by hand');
        image = imread([pathname filesep filename]);
    else
        image = double(imread(data.framePaths{1}{1}));
    end
    handle = figure;
    imshow(image,[]);
    %make mask
    maskHandCut = roipoly;
    close(handle)
    
    %refine using intensity
    
    areamask = getCellMask(data);
    mask = zeros(size(areamask));
    mask = areamask & maskHandCut;
else
    image = double(imread(data.framePaths{1}{1}));
    mask = ip.Results.mask;
end

%MAKE WINDOWS
goodWindows = 0;
perpSize = 100;
paraSize = 100;
while ~goodWindows
    %make polygons
    windowPolygons = getMaskWindows(mask,perpSize,paraSize);
    %visualize polygons
    handle  = imshow(image,[]);
    hold on
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
