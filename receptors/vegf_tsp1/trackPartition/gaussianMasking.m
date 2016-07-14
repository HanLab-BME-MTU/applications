function mask = gaussianMasking(movieInfo,MD,threshold,minSize,upscale)
%GAUSSIANMASKING Creates an image mask based on Gaussian detection
%results
%   mask = gaussianMasking(movieInfo,MD,threshold,minSize)
%   This function utilizes the Gaussian fitting information from the object
%   detection step to create a binary mask marking the location of each
%   object. A Gaussian is placed at each object location with the s.d.
%   obtained during fitting, and the result is thresholded. Specify a
%   minimum mask diameter in the fourth input argument (default 2 pixels).
%   Inputs:
%       movieInfo:      detection result from Gaussian fitting
%
%       MD:             Movie Data object 
%
%       threshold:      mask has value 'true' where the value of the
%                       Gaussian exceeds the given threshold. Set to 0 to
%                       use constant diameter masks only (specified by
%                       minSize) or use constantSizeMasking()
%
%       minSize:        minimum mask diameter in nm (default 100 nm)
%
%       upscale:        factor by which to upscale the mask from the
%                       size and resolution of the movie. This allows a
%                       more circular mask when dealing with mask diameters
%                       that are only a couple of pixels in the original
%                       movie image size. (default 1)
%
%   Output:
%       mask:           binary mask in a 3D matrix of the same size as the
%                       movie
%
%Kevin Nguyen, July 2016

% Default params
if nargin < 5
    upscale = 1;
end
if nargin < 4
    minSize = 100;
end

% Open up the movieInfo
nFramesDetected = size(movieInfo,1);
movieInfoCell = struct2cell(movieInfo)';
xData = movieInfoCell(:,1);
yData = movieInfoCell(:,2);
% ampData = movieInfoCell(:,3);

% Get information from movie data
imSize = MD.imSize_;
nFramesAll = MD.nFrames_;
if ~isempty(MD.pixelSize_)
    pixelSize = MD.pixelSize_;
else 
    pixelSize = 90; % default 90 nm pixel size
end

% Convert nm to pixels
minSizePx = round(minSize/pixelSize);

% Upscale
imSize = imSize*upscale;
minSizePx = minSizePx*upscale;

% Initialize empty matrix
mask = false(imSize(1),imSize(2),nFramesAll);

[xGrid,yGrid] = meshgrid(1:imSize(1),1:imSize(2));
% xGrid = xGrid(:);
% yGrid = yGrid(:);

parfor f = 1:nFramesDetected
    xList = xData{f}*upscale;
    yList = yData{f}*upscale;
    nObjects = size(xList,1);
    
    for i = 1:nObjects
        if threshold ~= 0
            maskNew = placeGaussian(xGrid,xList(i,1),xList(i,2),yGrid,yList(i,1),yList(i,2),threshold);
            minMask = placeConstantSize(xGrid,xList(i,1),yGrid,yList(i,1),minSizePx);
            maskNew = maskNew|minMask;
        else
            maskNew = placeConstantSize(xGrid,xList(i,1),yGrid,yList(i,1),minSizePx);
        end
        mask(:,:,f) = mask(:,:,f)|maskNew;
    end
end

end

function gMask = placeGaussian(xGrid,xC,xS,yGrid,yC,yS,threshold)
    gMask = (exp(-((xGrid-xC).^2)/(2*xS^2)-((yGrid-yC).^2)/(2*yS^2)))>threshold;
end

function cMask = placeConstantSize(xGrid,xC,yGrid,yC,diameter)
    r = round(diameter/2);
    cMask = (abs(xGrid-xC) <= r).*(abs(yGrid-yC) <= r);
end