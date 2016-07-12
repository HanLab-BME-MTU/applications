function mask = gaussianMaskingInner(movieInfo,MD,threshold,minSize)
%GAUSSIANMASKINGINNER Creates an image mask based on Gaussian detection
%results
%   mask = gaussianMaskingInner(movieInfo,MD,threshold,minSize)
%   This function utilizes the Gaussian fitting information from the object
%   detection step to create a binary mask marking the location of each
%   object. A Gaussian is placed at each object location with the s.d.
%   obtained during fitting, and the result is thresholded. Specify a
%   minimum mask diameter in the fourth input argument (default 2 pixels).
%   Inputs:
%       movieInfo:      detection result from Gaussian fitting
%       MD:             Movie Data object 
%       threshold:      mask has value 'true' where the value of the
%                       Gaussian exceeds the given threshold. Set to 0 to
%                       use constant diameter masks only (specified by
%                       minSize)
%       minSize:        minimum mask diameter (default 2 pixels)
%   Output:
%       mask:           binary mask in a 3D matrix of the same size as the
%                       movie

% Open up the movieInfo
nFramesDetected = size(movieInfo,1);
movieInfoCell = struct2cell(movieInfo)';
xData = movieInfoCell(:,1);
yData = movieInfoCell(:,2);
% ampData = movieInfoCell(:,3);

% Get information from movie data
imSize = MD.imSize_;
nFramesAll = MD.nFrames_;

% Initialize empty matrix
mask = false(imSize(1),imSize(2),nFramesAll);

[xGrid,yGrid] = meshgrid(1:imSize(1),1:imSize(2));
% xGrid = xGrid(:);
% yGrid = yGrid(:);

for f = 1:nFramesDetected
    xList = xData{f};
    yList = yData{f};
    nObjects = size(xList,1);
    assert(nObjects == size(yList,1),'xCoord and yCoord have lengths in frame %g',f)
    
    for i = 1:nObjects
        if threshold ~= 0
            maskNew = placeGaussian(xGrid,xList(i,1),xList(i,2),yGrid,yList(i,1),yList(i,2),threshold);
            minMask = placeConstantSize(xGrid,xList(i,1),yGrid,yList(i,1),minSize);
            maskNew = maskNew|minMask;
        else
            maskNew = placeConstantSize(xGrid,xList(i,1),yGrid,yList(i,1),minSize);
        end
        mask(:,:,f) = mask(:,:,f)|maskNew;
    end
end
end

function gMask = placeGaussian(xGrid,xC,xS,yGrid,yC,yS,threshold)
%     xMask = (xGrid > (xC-threshold*xS)).*(xGrid < (xC+threshold*xS));
%     yMask = (yGrid > (yC-threshold*yS)).*(yGrid < (yC+threshold*yS));
%     gMask = xMask.*yMask;
    gMask = (exp(-((xGrid-xC).^2)/(2*xS^2)-((yGrid-yC).^2)/(2*yS^2)))>threshold;
end

function cMask = placeConstantSize(xGrid,xC,yGrid,yC,diameter)
    cMask = ((xGrid-xC).^2 <= diameter).*((yGrid-yC).^2 <= diameter);
end