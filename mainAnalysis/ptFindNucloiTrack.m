function [coord, imgNucloiArea] = ptFindNucloiTrack (inputImage, nucloiLevel, nucloiMinSize, nucloiMaxSize, method)
% ptFindNucloiTrack detects dark areas and tries to fit cells into them
%
% SYNOPSIS       [coord, imgNucloiArea] = ptFindNucloiTrack (inputImage, nucloiLevel, nucloiMinSize, nucloiMaxSize, method)
%
% INPUT          inputImage    : either original image or segmented image (depends on method)
%                nucloiLevel   : level used for minima detection
%                nucloiMinSize : minimal size for nuclei
%                nucloiMaxSize : maximal size for nuclei
%                method        : 1 or 2. Says if clustering or image segmentation has been applied 
%                                (changes what ptFindNucloiTrack actually does)
%
% OUTPUT         coord         : coordinates found
%                imgNucloiArea : binary image showing the areas of nuclei
%
% DEPENDENCIES   ptFindNucloiTrack uses { nothing }
%                                  
%                ptFindNucloiTrack is used by { ptTrackCells
%                                               ptInitializeJob } 
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Do some input checking first to determine the segmentation method used
if method == 1          % Clustering
   % The input image is already segmented, so take the pixels with value 3
   imgNucloiArea = inputImage == 1;
    
elseif method == 2      % Segmentation
   % Calculate the binary image based on the nucloi intensity level
   imgNucloiArea = imextendedmin (inputImage, nucloiLevel);
    
else
   fprintf (1, 'ptFindNucloiTrack: invalid method. Please use method 1 or 2.\n');
end 

% Label the image
imgNucloiAreaLabeled = bwlabel (imgNucloiArea);

% Calculate nucloi properties using the regionprops function from the image toolbox
nucloiProperties = regionprops (imgNucloiAreaLabeled, 'Area', 'PixelList');

% Find really small and big regions and temporarily store this info
nucloiArea = [nucloiProperties.Area];
nucloiMinMaxArea = find ((nucloiArea < round (nucloiMinSize / 4 * 3)) | (nucloiArea > nucloiMaxSize));
tempProperties = nucloiProperties (nucloiMinMaxArea);

% Get all pixels within those regions
minMaxCoord = cat (1, tempProperties(:).PixelList);

% Set all of those pixels to zero so that we can forget about them
for iCount = 1 : size (minMaxCoord, 1)
   imgNucloiArea (minMaxCoord(iCount, 2), minMaxCoord(iCount, 1)) = 0;
end

% Close gaps which are smaller than a disk with radius 3
imgNucloiArea = imclose(imgNucloiArea, strel('disk',3));

% And dilate the picture as well (same disk size)
imgNucloiArea = imdilate(imgNucloiArea, strel('disk',3));

% Label results of the processed image
imgNucloiAreaLabeledBefore = bwlabel (imgNucloiArea);
lastNucloiGroup = max (max (imgNucloiAreaLabeledBefore)) ;

% More morphological operations, but this time the regions are already labeled.
% So if one region gets seperated thus becoming two or more, it still has the same numbers.
imgNucloiArea = imerode (imgNucloiArea, strel('disk',7));
imgNucloiArea = imdilate (imgNucloiArea, strel('disk',6));

% Now we label the image again. So we have one image labeled before the
% second morph and one labeled after
imgNucloiAreaLabeledAfter = bwlabel (imgNucloiArea);

counter = 0;
membership = [];

% For every group of imgNucloiAreaLabeledBefore (labeled before second morph oper) we allow
% one labeled group in imgNucloiAreaLabeledAfter (labeled after...). If there are more then 
% one region we selected the biggest.
for group = 1 : lastNucloiGroup
   counter = counter + 1;
   nucloiGroup = imgNucloiAreaLabeledAfter (find (imgNucloiAreaLabeledBefore == group));
   nucloiGroup = nucloiGroup (find (nucloiGroup));
   uniqueNucloiGroup = unique (nucloiGroup);
   if ~isempty (uniqueNucloiGroup)
      if length (uniqueNucloiGroup) > 1.1  
         [uniqueEntries, numberOfOccurences] = countEntries (nucloiGroup);
         [dummy, groupIndex] = max (numberOfOccurences);
         membership (group, 1) = uniqueEntries (groupIndex);
      else 
         membership (counter, 1) = uniqueNucloiGroup;
      end
   else
      counter = counter - 1;
   end
end

% In membership we have at most one group (belonging to imgNucloiAreaLabeledAfter) for
% every each group belonging to imgNucloiAreaLabeledBefore
% imgNucloiArea is overwritten now only including the groups (of imgNucloiAreaLabeledAfter) given in membership
if ~isempty (membership)
   % find all the non-zero groups
   membership = membership (find (membership));
   imgNucloiArea = ismember (imgNucloiAreaLabeledAfter, membership);

   % Label the new image again
   clear imgNucloiAreaLabeled;
   imgNucloiAreaLabeled = bwlabel (imgNucloiArea);
end

% Calculate nucloi properties using the regionprops function from the image toolbox
nucloiProperties = regionprops(imgNucloiAreaLabeled, 'Area', 'PixelList', 'MajorAxisLength', 'Centroid', 'Eccentricity');

% Remove really small and big areas
clear nucloiArea;
clear tempProperties;
nucloiArea = [nucloiProperties.Area];
nucloiMinArea = find (nucloiArea > nucloiMinSize);
tempProperties = nucloiProperties (nucloiMinArea);

clear nucloiArea;
nucloiArea = [tempProperties.Area];
nucloiMaxArea = find(nucloiArea < nucloiMaxSize);
tempProperties2 = tempProperties (nucloiMaxArea);

coord = round (cat (1, tempProperties2.Centroid));
