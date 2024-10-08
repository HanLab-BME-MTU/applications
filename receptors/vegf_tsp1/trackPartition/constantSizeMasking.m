function mask = constantSizeMasking(movieInfo,MD,diameter,upscale)
%CONSTANTSIZEMASKING Creates an image mask based on detection results
%   mask = constantSizeMasking(movieInfo,MD,diameter,upscale)
%   This function utilizes object detection results to create a binary mask 
%   marking the location of each object. 
%   Inputs:
%       movieInfo:      detection result from Gaussian fitting
%       MD:             Movie Data object 
%       diameter:       diameter of each mask in nm (default 100 nm)
%       upscale:        factor by which to upscale the mask from the
%                       size and resolution of the movie. This allows a
%                       more circular mask when dealing with mask diameters
%                       that are only a couple of pixels in the original
%                       movie image size. (default 1)
%   Output:
%       mask:           binary mask in a 3D matrix of the same size as the
%                       movie
% 
%Kevin Nguyen, July 2016

% Default params
if nargin < 4
    upscale = 1;
end
if nargin < 3
    diameter = 100;
end

% Open up the movieInfo
nFramesDetected = size(movieInfo,1);
movieInfoCell = struct2cell(movieInfo)';
xData = movieInfoCell(:,1);
yData = movieInfoCell(:,2);

% Get information from movie data
imSize = MD.imSize_;
nFramesAll = MD.nFrames_;
if ~isempty(MD.pixelSize_)
    pixelSize = MD.pixelSize_;
else 
    pixelSize = 90; % default 90 nm pixel size
end

% Convert nm to pixels
diameterPx = round(diameter/pixelSize);

% Upscale
imSize = imSize*upscale;
diameterPx = diameterPx*upscale;

% Initialize empty matrix
mask = false(imSize(1),imSize(2),nFramesAll);

[xGrid,yGrid] = meshgrid(1:imSize(1),1:imSize(2));

parfor f = 1:nFramesDetected
    xList = xData{f}*upscale;
    yList = yData{f}*upscale;
    nObjects = size(xList,1);
    assert(nObjects == size(yList,1),'xCoord and yCoord have different lengths in frame %g',f)
    
    % Iterate through each object and create masks
    for i = 1:nObjects
        maskNew = placeConstantSize(xGrid,xList(i,1),yGrid,yList(i,1),diameterPx);
%         maskNew = placeConstantSize2(xList(i,1),yList(i,1),diameterPx,imSize);
        mask(:,:,f) = mask(:,:,f)|maskNew;
    end
end
end


function cMask = placeConstantSize(xGrid,xC,yGrid,yC,diameter)
    r = round(diameter/2);
    cMask = sqrt(((xGrid-xC).^2)+((yGrid-yC).^2)) <= r;
end

% Use FFT for faster computation? Does not seem to improve much
function cMask = placeConstantSize2(xC,yC,diameter,imSize)
objects = zeros(imSize);
objects(round(yC),round(xC)) = 1;


objectsFFT = fft2(objects);

[kx,ky] = meshgrid(1:size(objects,2),1:size(objects,1));
radius = round(diameter/2);
center = round(size(objects)/2);
kernel = sqrt((kx-center(2)).^2+(ky-center(1)).^2) <= radius;
kernelFFT = fft2(kernel);
cMask = ifftshift(ifft2(kernelFFT.*objectsFFT));
end