function [cellProps, imgCellArea, imgLabeledCellArea] = ptCalculateCellArea (inputImage, coord, imgNucleiArea, imgHaloArea, distanceToCellArea, method)       
% ptCalculateCellArea  determines what areas of an image are occupied by cells and
%                      calculates image properties
%
% SYNOPSIS     [cellProps,imgCellArea,imgLabeledCellArea] = ptCalculateCellArea(inputImage,coord,imgNucleiArea,imgHaloArea,distanceToCellArea,method)
%
% INPUT        inputImage         : either original image (segmentation) or segmented image (clustering)
%              coord              : set of coordinates
%              imgNucleiArea      : binary image showing the areas of nuclei
%              imgHaloArea        : binary image showing the areas of halos
%              distanceToCellArea : distance a set of coordinates may have to an cell area and still belong to it
%              method             : Says if clustering (1) or image segmentation (2) has been applied
%                                   (which changes what ptCalculateCellArea actually does)
%
% OUTPUT       cellProps : 
%                 cellProps (:,1) = coord (:,1);
%	          cellProps (:,2) = coord (:,2);
%	          cellProps (:,3) = belongsToCluster (:);  (number of cluster - label)
%	          cellProps (:,4) = numberOfOccurences (:);  (how many cells in the cluster this cell is in)
%	          cellProps (:,5) = bodyCount (:);  (area of the cluster with the number given in belongsToCluster)
%	          cellProps (:,6) = perimeterDividedByArea (:);  (cluster)
%              imgCellArea : is the binary image of the areas occupied by cells                
%              imgLabeledCellArea : b/w labeled imgCellArea
%
%
% DEPENDENCIES ptCalculateCellArea uses {nothing}
%                                  
%              ptCalculateCellArea is used by { ptTrackCells }
%
% CHANGE REVISION
%
% Name                  Date            Comment
% --------------------- ----------      -----------------------------------------------
% Colin Glass           Feb 04          Initial version
% Andre Kerstens        03/31/04        Changed name from body to ptCalculateCellArea
%                                       and general cleanup of code

bodyDiskSize = 35;
background = []; 
I2 = []; 
I3 = [];
bw = []; 
imgCellArea = [];
imgLabeledCellArea = [];
belongsToCluster = [];

% Copy the coordinates into a temporary storage space
tempCoord = round (coord);

% Depending on which method of image analysis we used, there are two different approaches
if method == 1            % Clustering (which uses the segmented image)
   % The binary image is an or-function of the input image where pixels are 1 and the 
   % input image where pixels are 3
   imgCellArea = (inputImage==1 | inputImage==3);

elseif method == 2        % Segmentation (which uses a normal image)
   % Here we have to do a bit more since we're starting from a normal image
   % First we find the background by filtering out all the smaller 
   % objects (smaller than a disk with radius 35)
   %background = imopen (inputImage, strel ('disk', bodyDiskSize));

   % Then we subtract the background from the original image
   %I2 = imsubtract (inputImage,background); 

   % Normalize all the intensities to the range (0..1) ; stretchlim finds the maximum
   % intensities of the image which is needed as an input to imadjust
   %I3 = imadjust (I2, stretchlim(I2), [0 1]);

   % Find a global image threshold (between 0 and 1) using Otsu's method
   imgThreshold = graythresh (I3);
   
   % And finally generate a binary image out of all of this and combine (or function)
   % it with the binary images of the nuclei and the halos to show full cell areas
   imgCellArea = im2bw (I3, imgThreshold) | imgNucleiArea | imgHaloArea; 

else
   % The method number is unknown: print an error message
   printf (1, 'ptCalculateCellArea: invalid method %d. Please use 1 for clustering and 2 for segmentation.\n', method);
end

% Fill out any gaps by dilating and eroding the image a couple of times.
% AK: what are the disk numbers based on except for making the disk smaller and smaller?
imgCellArea = imdilate (imgCellArea, strel ('disk', 8));
imgCellArea = imerode  (imgCellArea, strel ('disk', 10));
imgCellArea = imdilate (imgCellArea, strel ('disk', 4));
imgCellArea = imerode  (imgCellArea, strel ('disk', 3));
imgCellArea = imdilate (imgCellArea, strel ('disk', 2));
imgCellArea = imerode  (imgCellArea, strel ('disk', 1));
imgCellArea = imdilate (imgCellArea, strel ('disk', 1));
imgCellArea = imdilate (imgCellArea, strel ('disk', 1));

% Another way of doing this could be the function imfill with parameter 'holes' maybe combined
% with dilation and erosion
%imgCellArea = imfill(bw,'holes');

% Label the objects in imgCellArea: 0 = background, 1 = first object, 2 = second object, etc
imgLabeledCellArea = bwlabel (imgCellArea);

% Prepare a matrix for the grouping of coordinates to objects
belongsToCluster = zeros (length (tempCoord), 1);

% Determine to which group each set of coordinates belongs
for i = 1 : length (tempCoord)          % process all cols in tempCoord (that's what length returns)
   if imgLabeledCellArea (tempCoord (i,2), tempCoord (i,1)) ~= 0
      % If the coordinates fall into a labeled area, that cluster label is stored in a matrix
      belongsToCluster(i) = imgLabeledCellArea (tempCoord(i,2), tempCoord(i,1));

   else
      % Initialize some variables
      helper = [];
      where = [];
      img_h = [];       % Image height
      img_w = [];       % Image width
      x_1 = [];
      y_1 = []; 
      x_2 = [];
      y_2 = [];
               
      % Get the size of the input image
      [img_h, img_w] = size (inputImage);

      % Calculate the coordinates of the point before the cell area
      x_1 = round (tempCoord (i,2) - distanceToCellArea);
      y_1 = round (tempCoord (i,1) - distanceToCellArea);

      % Calculate the coordinates of the point after the cell area
      x_2 = round (tempCoord (i,2) + distanceToCellArea);
      y_2 = round (tempCoord (i,1) + distanceToCellArea);
                   
      % Make sure all of the calculated coordinates are within image boundaries
      if x_1 < 1
         x_1 = 1;
      end 

      if x_2 > img_h
         x_2 = img_h;    
      end 

      if y_1 < 1
         y_2 = y_2 - y_1 + 1;
         y_1 = 1;
      end 
				
      if y_2 > img_w
         y_2 = img_w;
      end 

      % Now get the area that we have defined above
      where = imgLabeledCellArea (x_1:x_2, y_1:y_2);

      % And find all the entries that are different from 0
      helper = find (where);

      % In case there is something in the found part of the image do the following:
      if ~isempty (helper) 
         %index = [];
         %value = [];
         %uniqueIndex = [];
         %uniqueHelper = [];
         %numberOfOcc = [];

         % Sort the non-zero entries and find the unique values
         % uniqueIndex returns the last occurence of the respective unique entry
         helper = sort (helper);
         [uniqueHelper, uniqueIndex] = unique (helper);

         % Having sorted m before, we can now count the number of occurences
         if size (uniqueHelper, 1) > size (uniqueHelper, 2);
            uniqueIndex = [0 ; uniqueIndex];      % col vector
	 else
            uniqueIndex = [0 , uniqueIndex];      % row vector
	 end

	 numberOfOcc = diff (uniqueIndex); 
         [value, index] = max (numberOfOcc);
         belongsToCluster(i) = uniqueHelper (index);

         clear index;
         clear value;
         clear uniqueIndex;
         clear uniqueHelper;
         clear numberOfOcc;
      end

      clear helper;
      clear where;
      clear img_h
      clear img_w
      clear x_1
      clear x_2
      clear y_1
      clear y_2
   end
end

% Find all the pixels that belong to the background (aka not to an area of cells)
% and mark these empty; also unset the coordinates
onBackGround = find (belongsToCluster == 0);
if ~isempty (onBackGround)
   coord (onBackGround,:) = [];
   belongsToCluster (onBackGround) = [];
end
 
% Now we determine how many times a set of coordinates fall into the same
% labeled area: aka how many nuclei per area
sortedBelongsToCluster = sort (belongsToCluster);
[uniqueEntries, uniqueIdx] = unique (sortedBelongsToCluster);

% uniqueIdx returns the last occurence of the respective unique entry
% having sorted m before, we can now count the number of occurences
if size (uniqueEntries,1) > size (uniqueEntries,2);
   uniqueIdx = [0 ; uniqueIdx];
else
   uniqueIdx = [0 , uniqueIdx];
end 
numberOfOccurences = diff (uniqueIdx); 

noSensibleInfo = [];

bodyCount = zeros (length (uniqueEntries), 1);
for iCount = 1 : length (uniqueEntries)
   % Area of a group devided by the amount of nucloi in the group. Average
   onebody = [];
   onebody = imgLabeledCellArea == uniqueEntries(iCount);
   if length (find(onebody)) > 10
      areabod = length (find (onebody));
      perimeterDividedByArea (iCount) = length (find (bwperim (onebody))) / areabod;
      bodyCount (iCount) = areabod;
   else 
      noSensibleInfo (end+1, 1) = iCount;
      perimeterDividedByArea (iCount) = 0;
      bodyCount (iCount) = 0;
   end
end 

perimeterDividedByArea = sum (perimeterDividedByArea) / length (perimeterDividedByArea);

% One row of areaGroup is equivalent to the properties of one group (area)
areaGroup = zeros (length (uniqueEntries), 4);
areaGroup (:,1) = uniqueEntries (:);
areaGroup (:,2) = numberOfOccurences (:);
areaGroup (:,3) = bodyCount (:);
areaGroup (:,4) = perimeterDividedByArea (:);
areaGroup (noSensibleInfo,:) = [];
    
%one row of cellProps gives all information for one set of coordinates
cellProps = zeros (length (coord), 6);
cellProps(:,1) = coord(:,1);
cellProps(:,2) = coord(:,2);
cellProps(:,3) = belongsToCluster(:);

noSensibleProp = [];
% fill the right information into the right spots of cellProps, via belongsToCluster(now cellProps(:,3)) (cellProps<->areaGroup)
for i = 1 : length (coord)
   f = [];
   f = find (cellProps (i,3) == areaGroup (:,1));
   if ~isempty (f)
      cellProps (i,4) = areaGroup (f,2);
      cellProps (i,5) = areaGroup (f,3);
      cellProps (i,6) = areaGroup (f,4);
   else
      noSensibleProp (end+1,1) = i;
   end
end 

cellProps(noSensibleProp,:) = [];
clear f;

%avarea is the average area of all cells
%avarea = sum (bodyCount.*numberOfOccurences) / length (coord);
