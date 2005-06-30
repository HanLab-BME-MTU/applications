function [avgSNR,varargout] = fsmCalSNR(img,bitDepth,sDN,mask,fSize,roi,nGridY,nGridX)
%calSNR: Calculate the signal to noise ratio of an image.
%
% SYNOPSIS:
%    avgSNR = fsmCalSNR(img,bitDepth,sDN,mask);
%    [avgSNR,snr,y,x] = fsmCalSNR(img,bitDepth,sDN,mask,fSize,roi,nGridY,nGridX);
%
% INPUT:
%    img      : The input image.
%    bitDepth : The bit depth of the camera.
%    sDN      : Calibrated background noise from FSM experiment database.
%    mask     : Size of a mask for finding local maximum. It is
%               in the form [yLen,xLen] where 'yLen' and 'xLen' needs to be
%               odd. Otherwise, it is set to the next odd number. 
%               Default: [3 3]. Pass [] for default.
%    fSize    : Length scale of small features for calculating signal 
%               contrast. It needs to be odd. Otherwise, it is set to the 
%               next odd number. Default: 7 pixels. Pass [] for default.
%    roi      : The regeion of interest for calculating the average SNR.
%               It is either a string of value 'crop' (or 'whole') or 
%                a 2-column matrix in the form [xi yi].
%    nGridY
%    nGridX   : Number of grid points in Y and X direction for generating 
%               snr map.
%

if nargin < 4
   mask = [3 3];
end

if nargin < 5
   fSize = 7;
end

if nargin < 6
    roi = 'whole';
end

if nargin < 7
    nGridY = 10;
    nGridX = 10;
end

if nargout > 4
    error('Too many output arguments. See help calSNR.');
end

if isempty(mask)
    mask = [3 3];
end
if isempty(fSize)
    fSize = 7;
end

if mod(fSize,2) == 0
   fSize = fSize+1;
end

fSizeH = (fSize-1)/2;

figure; imshow(img,[]); hold on;

if ischar(roi)
   if strcmp(roi,'crop')
      [bw,xi,yi] = roipoly;
      roi = [xi yi];
   elseif strcmp(roi,'whole')
      roi = [];
   else
      error('ROI is not recogonized.');
   end
elseif ~isempty(roi) & isnumeric(roi)
   xi = roi(:,1);
   yi = roi(:,2);
elseif ~isempty(roi)
   error('ROI is not correctly defined. See help calSNR.');
end

if ~isempty(roi)
    plot(xi,yi,'r-.');
end

%Scale the image.
scaleI = double(img)./(2^bitDepth-1);
maxI = max(scaleI(:));
minI = min(scaleI(:));

%Smoth the image with a Gaussian kernal to get rid of dark noise and use the
%smoothed image for the calculation of local maximum.
sigma = 1;
smthI = imInterp(img,[],'Gaussian',sigma);

%Since 'smthI' is in the range [0 1], we need to scale it back to the range of
% 'scaleI' which is used in the calculation of 'sDN' by fsm calibration. 
smthI = smthI*(maxI-minI)+minI;

%Calculate local min and loc max intensity with mask size 'mask'.
%locMinI = locmin2d(smthI,mask);
locMaxI = locmax2d(smthI,mask);

%The non-zero values contain the original value of the image at that place.
%[pMinY,pMinX] = find(locMinI~=0);
[pMaxY,pMaxX] = find(locMaxI~=0);

% %Triangulization with [pMinX pMinY];
% triMin = delaunay(pMinX,pMinY);
% 
% %Get triangles for each local maximum;
% triIndex = tsearch(pMinX,pMinY,triMin,pMaxX,pMaxY);
% 
% %Get rid of local maximum that do not belong to any triangle.
% ind = find(isnan(triIndex));
% triIndex(ind) = [];
% pMaxX(ind)    = [];
% pMaxY(ind)    = [];
% 
% %Get the coordinates of the triangle for each maximum.
% triX = pMinX([triMin(triIndex,:),triMin(triIndex,1)]);
% triY = pMinY([triMin(triIndex,:),triMin(triIndex,1)]);

[m,n] = size(smthI);

%The local base intensity for each maximum is the average of the 3 closest
% local min intensity.
%locI0 = (smthI(m*(triX(:,1)-1)+triY(:,1))+smthI(m*(triX(:,2)-1)+triY(:,2))+ ...
%   smthI(m*(triX(:,3)-1)+triY(:,3)))/3;

%We use the local minimum in the box of size 'mask' around each local maximum
% as the base intensity.
locI0 = zeros(size(pMaxX));
for k = 1:length(pMaxX)
   locI0(k) = min(min(smthI(max(1,pMaxY(k)-fSizeH):min(m,pMaxY(k)+fSizeH), ...
      max(1,pMaxX(k)-fSizeH):min(n,pMaxX(k)+fSizeH))));
end

%The local signal is difference between the local maximum and local base
% intensity.
locMaxS = smthI(m*(pMaxX-1)+pMaxY);
locS = locMaxS-locI0;

%The local signal to noise ratio.
locSNR = locS/sDN;

%Get local maximum that are inside 'roi'.
if ~isempty(roi)
    in = inpolygon(pMaxX,pMaxY,xi,yi);
    inROI = find(in==1);
    avgSNR = mean(locSNR(inROI));
else
    avgSNR = mean(locSNR);
end

if nargout > 1
    %Interpolate to grid points.
    x = linspace(1,n,nGridX);
    y = linspace(1,m,nGridY);
    [gridX,gridY] = meshgrid(x,y);

    gridSNR = griddata(pMaxX,pMaxY,locSNR,gridX,gridY);

    snr = reshape(gridSNR,nGridY,nGridX);
    
    if nargout >= 2
        varargout{1} = snr;
    end
    
    if nargout == 4
        varargout{2} = y;
        varargout{3} = x;
    end
end
