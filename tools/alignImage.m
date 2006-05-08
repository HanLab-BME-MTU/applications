function alignShift = alignImage(refImgFile,alignImgPath,firstAlignImgFile,varargin)
%alignImage: This function aligns shifted images to the reference image according to markers in a
%            cropped region.
%
% SYNOPSIS: 
%    alignShift = alignImage(refImgFile,alignImgPath,firstAlignImgFile)
%    alignShift = alignImage(refImgFile,alignImgPath,firstAlignImgFile, ...
%       'modelChannel',modelChannel,'markerROI',markerROI)
%
% INPUT:
%    refImgFile       : The full path to the reference image file.
%    alignImgPath     : A cell array of directories where multi-channel images to be aligned are 
%                       stored. It can also be a string for single-channel images.
%    firstAlignImgFile: A cell array of the file names of the first images in each channel.
%
%    OPTIONAL:
%      par                           value
%    ---------------------------------------------------------------------------------------------
%    modelChannel: The image channel used for calculating the alignment shift. The calculated
%                  shift will be used for aligning the other channels. If it is empty, the
%                  first channel will be used as the model channel. If it is 0, then each
%                  channel will be aligned by iteself.
%    alignFrame : The list of frames to be aligned. Default, 0 for the whole stack.
%    markerROI   : An mx2 matrices that stores the (x,y) coordinates of the vertices of the
%                       polygon that defines the region of markers.
%    maxXShift   : Maximum shift in x-direction. Default, 50 pixels.
%    maxYShift   : Maximum shift in y-direction. Default, 50 pixels.
%
% OUTPUT: 
%    alignShift : A vector that stores the aligning shift of each image. If a model channel is used,
%                 it is a Nx2 matrix where N is the number of images. If each channel align by
%                 itself, it is a cell array of Nx2 matrix. The length of the cell is m where m is 
%                 the number of image channels. The aligned images will be stored in directories 
%                 whose names are composed by adding a suffix '_align' to the directory names 
%                 in 'alignImgPath'. For example, if the 'alignImgPath' is '/foo/imgChannel1', 
%                 then the output directory for aligned images is '/foo/imgChannel1_align'.
%
% Author: Lin Ji, May 04, 2006.

%Default parameters.
modelChannel = 0; %Each image path align by itself.
alignFrame   = 0;
markerROI    = []; %The whole image is used for alignment.
maxXShift    = 50; %Unit: pixels.
maxYShift    = 50; %Unit: pixels.

if nargin > 3
   for k = 1:2:nargin-3
      switch varargin{k}
         case 'modelChannel'
            modelChannel = varargin{k+1};
         case 'alignFrame'
            alignFrame = varargin{k+1};
         case 'markerROI'
            markerROI = varargin{k+1};
         case 'maxXShift'
            maxXShift = varargin{k+1};
         case 'maxYShift'
            maxYShift = varargin{k+1};
      end
   end
end
%Read the reference image.
refImg = imread(refImgFile);

if ~iscell(alignImgPath)
   alignImgPath = {alignImgPath};
end
numImgChannels = length(alignImgPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the image file list in each 'alignImgPath'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alignImgFileList = cell(1,numImgChannels);
numImages        = Inf;
for k = 1:numImgChannels
   alignImgFileList{k} = getFileStackNames([alignImgPath{k} filesep firstAlignImgFile{k}]);
   numImages = min(numImages,length(alignImgFileList{k}));
end
if alignFrame == 0
   alignFrame = 1:numImages;
end
numAlignFrames = length(alignFrame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align images in the model channel of 'alignImgPath' list to the reference image 'refImg'. The
% calculated shift for alignment is recorded for the other channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(num2str(numAlignFrames));
strForm = sprintf('%%.%dd',L);
fprintf(1,['Aligning frame (total: ' strForm '): '],numAlignFrames);
if modelChannel == 0
   %Each channel align by itself in this case.
   alignShift = zeros(numAlignFrames,2,numImgChannels);
   for k = 1:numImgChannels
      alignShift{k} = zeros(numAlignFrames,2);
      for ii = 1:numAlignFrames
         jj = alignFrame(ii);
         %Read each image in one channel and align it
         imgToAlign = imread(alignImgFileList{k}{jj});
         alignShift(ii,:,k) = calAlignShift(refImg,imgToAlign,markerROI, ...
            maxXShift,maxYShift);
      end
   end
else
   alignShift = zeros(numAlignFrames,2);
   for ii = 1:numAlignFrames
      jj = alignFrame(ii);
      procStatusStr = sprintf(strForm,jj);
      fprintf(1,procStatusStr);

      imgToAlign = imread(alignImgFileList{modelChannel}{jj});
      alignShift(ii,:) = calAlignShift(refImg,imgToAlign,markerROI, ...
         maxXShift,maxYShift);

      for k = 1:length(procStatusStr)
         fprintf(1,'\b');
      end
   end
   fprintf(1,[strForm '\n'],alignFrame(end));
end


function alignShift = calAlignShift(refImg,imgToAlign,markerROI,maxXShift,maxYShift)
%calAlignShift: This function calculate the shift required to align 'imgToAlign' to the
%               reference image 'refImg'.
% INPUT:
%    refImg: Reference image.
%    imgToAling: Image to be aligned to 'refImg'.
%    markerROI : The marker region of interest used in the correlation for calculating the shift.
%    maxXShift : The maximum shift in x direction. Default, 50 pixels.
%    maxYShift : The maximum shift in y direction. Default, 50 pixels.

[imgHeight,imgWidth] = size(refImg);

[pixX,pixY] = meshgrid(1:imgWidth,1:imgHeight);

if ~isempty(markerROI)
   IN = inpolygon(pixX,pixY,markerROI(:,1),markerROI(:,2));
   roiIN = find(IN==1);
else
   roiIN = 1:imgHeight*imgWidth;
end
refImgROI = double(refImg(roiIN));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Debugging:
%markerMask        = zeros(size(refImg));
%markerMask(roiIN) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

score = zeros(2*maxYShift+1,2*maxXShift+1);
for xShift = -maxXShift:maxXShift
   for yShift = -maxYShift:maxYShift
      %After shift, 'markerROI' may be outside of the image. Identify pixel index of 'markerROI' 
      % that are inside the image.
      in = find(pixX(roiIN)+xShift>0 & pixX(roiIN)+xShift<=imgWidth & ...
         pixY(roiIN)+yShift > 1 & pixY(roiIN)+yShift <= imgHeight);
      roiShiftIN = sub2ind(size(refImg),pixY(roiIN(in))+yShift,pixX(roiIN(in))+xShift);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Debugging:
      %imgToAlignMask = imgToAlign;
      %imgToAlignMask(roiShiftIN) = 0;
      %[xShift yShift]
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      imgToAlignROI = double(imgToAlign(roiShiftIN));

      refImgROINorm     = sqrt(sum(refImgROI(in).^2));
      imgToAlignROINorm = sqrt(sum(imgToAlignROI.^2));

      ii = yShift + maxYShift + 1;
      jj = xShift + maxXShift + 1;
      score(ii,jj) = sum(refImgROI(in).*imgToAlignROI)/refImgROINorm/imgToAlignROINorm;
   end
end
[maxS,maxIJ] = max(score(:));
[maxI,maxJ]  = ind2sub(size(score),maxIJ);

alignShift = [maxJ-maxXShift-1 maxI-maxYShift-1];

return;
