function [coord, imgHalo]= ptFindHalos (inputImage, erosionDiskSize, haloLevel, method)
% ptFindHalos detects light areas and tries to fit calls into them
%
% SYNOPSIS       [coord,imgHalo] = ptFindHalos (inputImage, erosionDiskSize, haloLevel, method)
%
% INPUT          inputImage    : either original image or segmented image (depends on method)
%                erosionDiskSize : size of disk to erode with
%                haloLevel     : threshold
%                method        : 1 or 2. Says if clustering or image segmentation has been applied applied
%                                (changes what ptFindHalos actually does)
%
% OUTPUT         coord       : coordinates that have been found
%                imgHalo : binary image giving the areas of our halos
%
% DEPENDENCIES   ptFindHalos uses { nothing }
%                                  
%                ptFindHalos is used by { ptTrackCells
%                                         ptInitializeJob }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Do some input checking first to determine the segmentation method used
if method == 1          % Clustering
   % The input image is already segmented, so take the pixels with value 3 
   imgHalo = inputImage == 3;

elseif method == 2 || method == 3     % Segmentation
   % Prepare storage space
   index = [];
   imgHaloLabeled = [];
   erodedHaloArea = [];
   uniqueErodedHaloArea = [];
   imgHalo = zeros (size (inputImage, 1), size (inputImage, 2));

   % Look for pixels that are above a certain value provided by haloLevel
   haloIndex = find (inputImage >= haloLevel);
   imgHalo (haloIndex) = 1;
else
   fprintf (1, 'ptFindHalos: invalid method. Please use method 1 or 2.\n');
end 

% Let's do some morphological ops now to clean the image up
imgHalo = imfill (imgHalo,'holes');
imgHalo = imerode (imgHalo, strel ('disk', 5));
imgHalo = imdilate (imgHalo, strel ('disk', 5));
imgHalo = bwareaopen (imgHalo, 50);

% Close small gaps using a disk with radius 3
%imgHalo = imclose (imgHalo, strel ('disk', 3));

% Label groups
imgHaloLabeled = bwlabel (imgHalo);

% Erode the labeled image. In this way even big groups are reduced to nothing, if
% their area is not coherent
%erodedHaloArea = imerode (imgHaloLabeled, strel ('disk',erosionDiskSize));

% See what groups are still present
%uniqueErodedHaloArea = unique (erodedHaloArea);
%uniqueErodedHaloArea (1,:) = [];

% Calculate nucloi properties using the regionprops function from the image toolbox
%haloProperties = regionprops (erodedHaloArea, 'Area', 'PixelList', 'MajorAxisLength', 'Centroid', 'Eccentricity');
haloProperties = regionprops (imgHaloLabeled, 'Area', 'PixelList', 'MajorAxisLength', 'Centroid', 'Eccentricity');

% Get the centroid coordinates. Note that the centroids are of the eroded picture.
% Generally properties of the eroded picture are NOT representative!
%uniqueHaloProperties = haloProperties (uniqueErodedHaloArea);
%coord = round (cat (1, uniqueHaloProperties.Centroid));
coord = round (cat (1, haloProperties.Centroid));
