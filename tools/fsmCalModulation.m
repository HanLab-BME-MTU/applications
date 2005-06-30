function [modS,varargout] = fsmCalModulation(img,varargin)
%calModulation: Calculate the average local feature modulation of an image.
%
% SYNOPSIS:
%    modS = calModulation(img,mask);
%    modS = calModulation(img,mask,fSize);
%    [modS,roi] = calModulation(img,mask,fSize,roi);
% 
% INPUT:
%    img : The image to be analyzed.
%
%    Optional input:
%    mask  : A 2-element vector in the form [yLen xLen] that is used in 
%            finding local maxima. Default: [3 3] (pixels). 
%            See help locMax2d.
%    fSize : Local feature size. Default: 7 pixels. It is base on a typical
%            speckle size in FSM images.
%    roi   : region of interest. It can be either string of value 'crop' 
%            or it is a 2-columns vector in the form [xi yi] that 
%            specifies the polygon. Default value: 'crop'.

mask  = [3 3]; % pixels.
fSize = 7; % pixels.

roi = 'whole';

if nargin > 4
   error('Too many input arguments.');
end

if nargout > 2
   error('Too many output arguments.');
end

if nargin > 1
   mask = varargin{1};
end

if nargin > 2
   fSize = varargin{2};
end

if nargin > 3
   roi = varargin{3};
end

%You can also pass [] for the default.
if isempty(mask)
    mask = [3 3];
end
if isempty(fSize)
    fSize = 7;
end

if mod(fSize,2) == 0
    fSize = fSize+1;
end

%Half of feature size.
fSizeH = (fSize-1)/2;

figure;
imshow(img,[]); hold on;

if ischar(roi)
   if strcmp(roi,'crop')
      [bw,xi,yi] = roipoly;
      roi = [xi yi];
   elseif strcmp(roi,'whole')
      roi = [];
   else
      error('ROI is not recogonized.');
   end
elseif isnumeric(roi)
   xi = roi(:,1);
   yi = roi(:,2);
else
   error('ROI is not correctly defined. See help calModulation.');
end

if ~isempty(roi)
   plot(xi,yi,'r-.');
end

if nargout == 2
   varargout{1} = roi;
end

%Smoth the image with a Gaussian kernal to get rid of dark noise and use the
%smoothed image for the calculation of local maximum.
sigma = 1;
smthI = imInterp(img,[],'Gaussian',sigma);

%Calculate local min and loc max intensity map with mask size 'mask'.
locMaxMap = locmax2d(smthI,mask);

%The non-zero values contain the original value of the image at that place.
[pMaxY,pMaxX] = find(locMaxMap~=0);

%Get local maximum that are inside 'roi'.
if ~isempty(roi)
   in = inpolygon(pMaxX,pMaxY,xi,yi);
   pMaxX = pMaxX(find(in==1));
   pMaxY = pMaxY(find(in==1));
end
[m,n] = size(smthI);

%We use the local minimum in the box of size 'mask' around each local 
% maximum point as the local base intensity, 'locI0'. And, we use the 
% difference between the local maximum intensity and local minimum 
% intensity in a box of size double of 'mask' around each point as the 
% reference variation for that point.
quatMinI = zeros(1,4);
locI0 = zeros(size(pMaxX));
refI = ones(size(pMaxX));
%refIV = ones(size(pMaxX));
for k = 1:length(pMaxX)
   yInd1 = max(1,pMaxY(k)-fSizeH);
   yInd2 = min(m,pMaxY(k)+fSizeH);
   xInd1 = max(1,pMaxX(k)-fSizeH);
   xInd2 = min(n,pMaxX(k)+fSizeH);
   quatMinI(1) = min(min(smthI(yInd1:pMaxY(k),pMaxX(k):xInd2)));
   quatMinI(2) = min(min(smthI(pMaxY(k):yInd2,pMaxX(k):xInd2)));
   quatMinI(3) = min(min(smthI(pMaxY(k):yInd2,xInd1:pMaxX(k))));
   quatMinI(4) = min(min(smthI(yInd1:pMaxY(k),xInd1:pMaxX(k))));
   locI0(k) = mean(quatMinI);
   
   %The reference intensity is calulated as the minimum intensity in a
   % local box region of size 5*fSize around each local maximum.
   yInd1 = max(1,pMaxY(k)-5*fSizeH);
   yInd2 = min(m,pMaxY(k)+5*fSizeH);
   xInd1 = max(1,pMaxX(k)-5*fSizeH);
   xInd2 = min(n,pMaxX(k)+5*fSizeH);
   refI(k) = min(min(smthI(yInd1:yInd2,xInd1:xInd2)));
   
%    refIV(k) = max(max(smthI(max(1,pMaxY(k)-3*mask(1)):min(m,pMaxY(k)+3*mask(1)), ...
%       max(1,pMaxX(k)-3*mask(2)):min(n,pMaxX(k)+3*mask(2))))) - ...
%       min(min(smthI(max(1,pMaxY(k)-3*mask(1)):min(m,pMaxY(k)+3*mask(1)), ...
%       max(1,pMaxX(k)-3*mask(2)):min(n,pMaxX(k)+3*mask(2)))));
end

%The local signal variation is defined by the difference between the local maximum and local base
% intensity.
locMaxI = smthI(m*(pMaxX-1)+pMaxY);
locIV = locMaxI-locI0;

%And, the averaged local modulation is defined by
%modS = sum(locIV./locMaxI)/length(locIV);
modS = sum(locIV./(locMaxI-refI))/length(locIV);

