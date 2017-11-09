function [SegmentationCoords] = segmentationCoordinates(originalImage,maskofSegments,blobOfInterest,...
              sizeThreshold,intensityThreshold,plot, varargin)
%segmentationCoordinates: get the coordinates of the selected segmentations
%      
%
% synopsis function [maskBlobs,SegmentationCoords] = classify(OriginalImage,...
%   image,BlobOfInterest,SizeThreshold,IntensityThreshold,labels,imageNorm)
%
% Inputs:
%         originalImage: image used to create segmentations.
%        maskOfSegments: mask from segmentations.
%        blobOfInterest: 'all' for all segmentations  
%                        'big' for segementations bigger than sizeThreshold
%                        or intensityThreshold
%                        'small' for segmentations smaller than
%                        sizeThreshold or intensityThreshold
%                        Optional. Default: ;'all'
%         sizeThreshold: size threshold for segmentations.
%                        Optional: Default:100
%    intensityThreshold: intensity threshold for segmentations.                
%                 plot : to plot selected segmentations  1. 
%                        no plotting                     0.
%                        Optional. Default: 0
%   varargin   (labels): labels from segementations. 
%                        Optional. Needed when blobOfInterest is 'big' or 'small' 
%           (imageNorm): normalized original image.
%                        Optional. Only needed when using integrated
%                        intensity
%
% Outputs:
%         SegmentationCoords: coordinates of the selected segmentations
%
% Jesus Vega July 2017

%% Input

%check inputs
if nargin < 1
    disp('Enter original image');
end

if nargin < 2
    disp('Enter mask from segmentations');
    return
end

% blobOfInterest
if nargin < 3 || isempty(blobOfInterest)
    blobOfInterest = 'all';
end

% if both thresholds are empty use a predetermined sizeThreshold
if nargin < 4 || isempty(sizeThreshold) && isempty(intensityThreshold)
    sizeThreshold = 100;
end

%if both thresholds are full use sizeThreshold
if ~isempty(sizeThreshold) && ~isempty(intensityThreshold)
    intensityThreshold = [];
end

% intensityThreshold
if ~isempty(intensityThreshold)
    
    %check that you have labels and normalized image
    if isempty(varargin) || numel(varargin) < 2
         disp('Enter normalized image and/or labels');
         return
    else 
        
      stats = regionprops(maskofSegments,varargin{2},'PixelValues');
        
      Integrated_Int = zeros(size(stats,1),1);
       
      %calculate integrated intensity
       for k = 1:size(stats,1)
           Integrated_Int(k) = sum(stats(k).PixelValues);
       end 
    end
end

if nargin < 7 || isempty(plot)
    plot = 0;
end

%% Separation

%separate segmentations
switch blobOfInterest
    
    case 'all'
        maskBlobs = maskofSegments;
        
    case 'big'    
        if ~isempty(sizeThreshold)
            
            stats = regionprops(maskofSegments,'Area');
            idx = find([stats.Area] > sizeThreshold);
            maskBlobs = ismember(varargin{1}, idx);
       
        else
            idx = find(Integrated_Int > intensityThreshold);
            maskBlobs = ismember(varargin{1}, idx);
        end
        
        
    case 'small'
        if ~isempty(sizeThreshold)
            
            stats = regionprops(maskofSegments,'Area');
            idx = find([stats.Area] <= sizeThreshold);
            maskBlobs = ismember(varargin{1}, idx);
        
        else
            idx = find(Integrated_Int <= intensityThreshold);
            maskBlobs = ismember(varargin{1}, idx);
        end
end

%Number of segmentations 
% NumOfSegments = numel(idx);
% disp(['Number of Segementations: ' num2str(NumOfSegments)]); 

%% Getting Coordinates 

    SegIdx = regionprops(maskBlobs,'PixelIdxList');
    
    % Concatenate Idxs from position
    SegIdxList = vertcat(SegIdx.PixelIdxList);
    
    % Convert indeces into coordinates
    [SegRow, SegCol] = ind2sub(size(maskofSegments),SegIdxList);
    SegmentationCoords = [SegRow SegCol];
    
%% Plot

 if plot
   originalImage = double(originalImage);

   imageScaled = (originalImage - prctile(originalImage(:),1)) / (prctile(originalImage(:),99) - prctile(originalImage(:),1));
    imageScaled(imageScaled<0) = 0;
    imageScaled(imageScaled>1) = 1;
    
    %get the blob edges from the final blob mask
    SE = strel('square',3);
    edgesBlobs = imdilate(maskBlobs,SE) - maskBlobs;
    
    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;
    
    
    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
         
    figure; imshow(image3Color,[]);
    
 end
end


