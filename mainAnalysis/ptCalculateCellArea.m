function [cellProps, clusterProps, imgCellArea, imgLabeledCellArea] = ptCalculateCellArea (inputImage, coord, imgNucleiArea, imgHaloArea, distanceToCellArea, minSizeNucleus, maxSizeNucleus)
% ptCalculateCellArea  determines what areas of an image are occupied by cells and
%                      calculates image properties
%
% SYNOPSIS     [newCoord, cellProps,imgCellArea,imgLabeledCellArea] = ptCalculateCellArea (inputImage, 
%              coord, imgNucleiArea, imgHaloArea, distanceToCellArea, minSizeNucleus, maxSizeNucleus)
%
% INPUT    inputImage     : either original image (segmentation) or segmented image (clustering)
%          coord          : set of coordinates
%          imgNucleiArea  : binary image showing the areas of nuclei
%          imgHaloArea    : binary image showing the areas of halos
%          distanceToCellArea : distance a set of coordinates may have to an cell area and still belong to it
%	   minSizeNucleus : everything smaller than this size will be removed from the binary image
%	   maxSizeNucleus : everything bigger than this size will be removed from the binary image
%
% OUTPUT       cellProps :
%                 cellProps (:,1) = coord (:,1);
%	          cellProps (:,2) = coord (:,2);
%	          cellProps (:,3) = clusterMember (:);  (number of cluster - label)
%	          cellProps (:,4) = numberOfOccurences (:);  (how many cells in the cluster this cell is in)
%	          cellProps (:,5) = clusterArea (:);  (area of the cluster with the number given in clusterMember)
%	          cellProps (:,6) = perimeterDividedByArea (:);  (cluster)
%              clusterProps :
%                 clusterProps (:,1) = uniqClusterMembers (:);          (which cells are in this cluster)
%                 clusterProps (:,2) = numberOfOccurences (:);          (how many cells in this cluster)
%                 clusterProps (:,3) = clusterArea (:);                 (the area of this cluster)
%                 clusterProps (:,4) = perimeterDividedByArea (:);      (a value for perimeter / area)
%                 clusterProps (:,5) = clusterPerimeter (:);            (the length of the perimeter)
%              imgCellArea : is the binary image of the areas occupied by cells
%              imgLabeledCellArea : b/w labeled imgCellArea
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
% Andre Kerstens        May 04          Changed name from body to ptCalculateCellArea
%                                       and general cleanup of code

% The binary image is an or-function of the input image where pixels are 1 (nuclei) and the
% input image where pixels are 3 (halos)
imgCellArea = (inputImage == 1 | inputImage == 3);

% Fill out any gaps by dilating and eroding the image a couple of times.
% AK: what are the disk numbers based on except for making the disk smaller and smaller?
imgCellArea = imdilate (imgCellArea, strel ('disk', 8));
imgCellArea = imerode  (imgCellArea, strel ('disk', 8));
imgCellArea = imdilate (imgCellArea, strel ('disk', 4));
imgCellArea = imerode  (imgCellArea, strel ('disk', 3));
imgCellArea = imdilate (imgCellArea, strel ('disk', 2));
imgCellArea = imerode  (imgCellArea, strel ('disk', 1));
imgCellArea = imdilate (imgCellArea, strel ('disk', 1));

% Get rid of all the really small objects (smaller than a nucleus)
imgCellArea = bwareaopen (imgCellArea, minSizeNucleus);

% Process the binary image to find new coordinates. This is basically done
% by subtracting circular areas that are centered around the already found
% coordinates and see if any areas are left that are big enough to be cells themselves
%newNucCoord = ptSubtractCellArea (imgCellArea, coord, maxSizeNucleus);

% Add the newly found coordinates to the ones we already have
%newCoord = cat (1, coord, newNucCoord);

% Copy the coordinates into a temporary storage space
newCoord = coord;
tempCoord = newCoord;

% Label the objects in imgCellArea: 0 = background, 1 = first object, 2 = second object, etc
imgLabeledCellArea = bwlabel (imgCellArea);

% Prepare a matrix for the grouping of coordinates to objects
clusterMember = zeros (length (tempCoord), 1);

% Determine to which group each set of coordinates belongs
for i = 1 : length (tempCoord)          % process all rows in tempCoord

   % Do the coordinates fall in a labeled area at all?
   if imgLabeledCellArea (tempCoord(i,2), tempCoord(i,1)) ~= 0
      % If the coordinates fall into a labeled area, that cluster label is stored in a matrix
      clusterMember(i) = imgLabeledCellArea (tempCoord(i,2), tempCoord(i,1));
   else
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
      searchArea = imgLabeledCellArea (x_1:x_2, y_1:y_2);

      % And find all the entries that are different from 0
      labelArea = find (searchArea);

      % In case there is something in the found part of the image do the following:
      if ~isempty (labelArea) 
         % Sort the non-zero entries and find the unique values
         % uniqLabelIndex returns the last occurence of the respective unique entry
         labelArea = sort (labelArea);
         [uniqLabelArea, uniqLabelIndex] = unique (labelArea);

         % Having sorted labelArea before, we can now count the number of occurences
         if size (uniqLabelArea, 1) > size (uniqLabelArea, 2);
            uniqLabelIndex = [0 ; uniqLabelIndex];      % col vector
	 else
            uniqLabelIndex = [0 , uniqLabelIndex];      % row vector
	 end

         % Get the maximum label value which will the cluster the cell belongs to
	 numberOfOcc = diff (uniqLabelIndex); 
         [value, index] = max (numberOfOcc);
         clusterMember(i) = uniqLabelArea (index);
      end
   end
end

% Find all the cell coordinates that are part of the background (aka not of an area of cells)
% and mark these empty; also unset the coordinates
onBackGround = find (clusterMember == 0);
if ~isempty (onBackGround)
   coord (onBackGround,:) = [];
   clusterMember (onBackGround) = [];
end
 
% Determine how many times a set of coordinates falls into the same
% labeled area (aka how many nuclei per area)
clusterMemberSorted = sort (clusterMember);
[uniqClusterMembers, uniqClusterIndex] = unique (clusterMemberSorted);

% uniqClusterIndex returns the last occurence of the respective unique entry
% having sorted clusterMember before, we can now count the number of occurences
if size (uniqClusterMembers,1) > size (uniqClusterMembers,2);
   uniqClusterIndex = [0 ; uniqClusterIndex];
else
   uniqClusterIndex = [0 , uniqClusterIndex];
end 

numberOfOccurences = diff (uniqClusterIndex); 

% Initialize the counter for disregarded clusters
disregardCount = 1;
disregardedClusters = [];

% How many clusters do we have?
clusterArea = zeros (length (uniqClusterMembers), 1);

% Calculate cluster properties
for iCount = 1 : length (uniqClusterMembers)
   % Get the area of the next cluster in the labeled image
   cluster = imgLabeledCellArea == uniqClusterMembers (iCount);

   % Calculate the area of this cluster
   clusterAreaTemp = length (find (cluster));
   
   % Is the area at least as big as a nucleus (1 cell)? If no disregard it.
   if clusterAreaTemp > minSizeNucleus
      % Calculate the perimeter and perimeter area of this cluster. The function bwperim will take holes
      % into account as well
      clusterPerimeter (iCount) = length (find (bwperim (cluster)));
      perimeterDividedByArea (iCount) = clusterPerimeter (iCount) / clusterAreaTemp;
      clusterArea (iCount, 1) = clusterAreaTemp;
   else 
      % The cluster is too small to be of interest, so set all values to 0
      disregardedClusters (disregardCount) = iCount;
      disregardCount = disregardCount + 1;

      clusterPerimeter (iCount) = 0;
      perimeterDividedByArea (iCount) = 0;
      clusterArea (iCount, 1) = 0;
   end
end 

% Calculate the average perimeter length and area of the clusters
if length (perimeterDividedByArea) > 0
   perimeterDividedByArea = sum (perimeterDividedByArea) / length (perimeterDividedByArea);
else
   perimeterDividedByArea = 0;
end
   
% One row of clusterProps is equivalent to the properties of one cluster
clusterProps = zeros (length (uniqClusterMembers), 4);
clusterProps (:,1) = uniqClusterMembers (:);
clusterProps (:,2) = numberOfOccurences (:);
clusterProps (:,3) = clusterArea (:);
clusterProps (:,4) = perimeterDividedByArea (:);
clusterProps (:,5) = clusterPerimeter (:);
clusterProps (disregardedClusters,:) = [];
    
% One row of cellProps gives all information for one set of coordinates
cellProps = zeros (length (newCoord), 6);
cellProps(:,1) = newCoord(:,1);
cellProps(:,2) = newCoord(:,2);
cellProps(:,3) = clusterMember(:);

% Initialize the counter for disregarded cells
disregardCount = 1;
disregardedCells = [];

% Put the right information into the right rows of cellProps (via clusterMember, now cellProps(:,3)) 
% (cellProps <-> clusterProps)
for jCount = 1 : length (coord)
   index = find (cellProps (jCount,3) == clusterProps (:,1));
   if ~isempty (index)
      cellProps (jCount, 4) = clusterProps (index, 2);
      cellProps (jCount, 5) = clusterProps (index, 3);
      cellProps (jCount, 6) = clusterProps (index, 4);
   else
      disregardedCells (disregardCount, 1) = jCount;
      disregardCount = disregardCount + 1;
   end
end

% Disregard the cells which are not member of a cluster
cellProps (disregardedCells,:) = [];
