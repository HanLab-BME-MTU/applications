function [cellProps, clusterProps] = ptCalculateCellAreaUsingVariance (edgeImageLabeled, coord, distanceToCellArea, minSizeNucleus, edgeKernelSize)
% ptCalculateCellAreaUsingVariance  determines what areas of an image are occupied by cells and
%                                   calculates image properties
%
% SYNOPSIS     [cellProps, clusterProps] = ptCalculateCellAreaUsingVariance (edgeImage, 
%                                          coord, distanceToCellArea, minSizeNucleus, edgeKernelSize)
%
% INPUT    edgeImageLabeled : a labeled binary image showing the area of cells and clusters
%          coord          : set of coordinates
%          distanceToCellArea : distance a set of coordinates may have to an cell area and still belong to it
%	       minSizeNucleus : everything smaller than this size will be removed from the binary image
%          edgeKernelSize : size of the kernel that was used to create the edge variance image
%
% OUTPUT       cellProps :
%                 cellProps (:,1) = coord (:,1);
%	              cellProps (:,2) = coord (:,2);
%	              cellProps (:,3) = clusterNr (:);  (number of cluster - label)
%
%              clusterProps :
%                 clusterProps (:,1) = uniqClusterNr (:);          (which cells are in this cluster)
%                 clusterProps (:,2) = numberOfCells (:);          (how many cells in this cluster)
%                 clusterProps (:,3) = clusterArea (:);            (the area of this cluster)
%                 clusterProps (:,4) = clusterPerimeter (:);       (the length of the perimeter)
%                 clusterProps (:,5) = clusterPerimeterElements (:); (how many elements does the perimeter exist of)
%
% DEPENDENCIES ptCalculateCellAreaUsingVariance uses {nothing}
%
%              ptCalculateCellAreaUsingVariance is used by { ptTrackCells }
%
% CHANGE REVISION
%
% Name                  Date            Comment
% --------------------- ----------      -----------------------------------------------
% Colin Glass           Feb 04          Initial version
% Andre Kerstens        May 04          Changed name from body to ptCalculateCellArea
%                                       and general cleanup of code
% Andre Kerstens        Jun 04          Calculations are done based on edges found in variance image

% % Process the edge image so that we end up with a labeled image of areas
% % First remove the white edge that was created during convolution
% edgeImageNoEdge = edgeImage(edgeKernelSize:end-edgeKernelSize, edgeKernelSize:end-edgeKernelSize);
% 
% % Some morphological operations to increase cluster quality
% edgeImageOpened = bwareaopen (edgeImageNoEdge, minSizeNucleus);
% edgeImageClosed = imfill (edgeImageOpened, 'holes');
% edgeImageBridged = bwmorph (edgeImageClosed,'bridge');
% edgeImageFilled = imfill (edgeImageBridged, 'holes');
% 
% % Label the image
% edgeImageLabeled = bwlabel (edgeImageFilled);

% Prepare a matrix for the grouping of coordinates to objects
clusterNr = zeros (length (coord), 1);

% Determine to which group each set of coordinates belongs
for iCount = 1 : length (coord)          % process all rows in coord

   % Store the coordinates
   x = coord(iCount,2);
   y = coord(iCount,1);
   
   % Do the coordinates fall in a labeled area at all?
   if edgeImageLabeled (x,y) ~= 0
      % If the coordinates fall into a labeled area, that cluster label is stored in a matrix
      clusterNr(iCount) = edgeImageLabeled (x,y);
   else
      % We also take cells into account that are near enough a cell area
      
      % Get the size of the input image
      [img_h, img_w] = size (edgeImageLabeled);

      % Calculate the coordinates of the point before the cell area
      x_1 = round (coord (iCount,2) - distanceToCellArea);
      y_1 = round (coord (iCount,1) - distanceToCellArea);

      % Calculate the coordinates of the point after the cell area
      x_2 = round (coord (iCount,2) + distanceToCellArea);
      y_2 = round (coord (iCount,1) + distanceToCellArea);
                   
      % Make sure all of the calculated coordinates are within image
      % boundaries (minus the edge lost by convolving)
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
      searchArea = edgeImageLabeled (x_1:x_2, y_1:y_2);

      % And find all the entries that are different from 0
      labelArea = searchArea(find (searchArea));

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

         % Get the maximum label value which the cell belongs to
	     numberOfOcc = diff (uniqLabelIndex); 
         [value, index] = max (numberOfOcc);
         clusterNr(iCount) = uniqLabelArea (index);
      end   % if ~isempty (labelArea)
   end   % if edgeImageLabeled
end   % for iCount = 1

% % Find all the cell coordinates that are part of the background (aka not of an area of cells)
% % and mark these empty; also unset the coordinates
% onBackGround = find (clusterNr == 0);
% if ~isempty (onBackGround)
%    coord (onBackGround,:) = [];
%    clusterNr (onBackGround) = [];
% end
 
% Determine how many times a set of coordinates falls into the same
% labeled area (aka how many nuclei per area)
clusterNrSorted = sort (clusterNr);
[uniqClusterNr, uniqClusterIndex] = unique (clusterNrSorted);

% uniqClusterIndex returns the last occurence of the respective unique entry
% having sorted clusterNr before, we can now count the number of cells per cluster
% Before doing this we have to add a 0 entry to be able to calculate diffs later on
if size (uniqClusterNr,1) > size (uniqClusterNr,2);
   uniqClusterIndex = [0 ; uniqClusterIndex];
else
   uniqClusterIndex = [0 , uniqClusterIndex];
end 

% Get the number of cells in every cluster by cleverly using the diff
% function to calculate the differenc between indexes
numberOfCells = diff (uniqClusterIndex); 

% Calculate cluster properties
for iCount = 1 : length (uniqClusterNr)
   % Get the area of the next cluster in the labeled image
   cluster = edgeImageLabeled == uniqClusterNr (iCount);

   % Calculate the area of this cluster
   clusterArea (iCount) = length (find (cluster));
   
   % Calculate the perimeter and perimeter area of this cluster. The function bwperim will take holes
   % into account as well
   clusterPerimeter (iCount) = length (find (bwperim (cluster)));
   
   % Inverse the binary cluster image with the object of interest
   clusterInv = ~cluster;
   
   % Label this image: background and holes in the object will get a number
   clusterInvLabel = bwlabel (clusterInv);
   
   % To find the number of perimeter elements of the object take the max
   % label nr
   clusterPerimeterElements (iCount) = max (max (clusterInvLabel));
end
   
% One row of clusterProps is equivalent to the properties of one cluster
clusterProps = zeros (length (uniqClusterNr), 5);
clusterProps (:,1) = uniqClusterNr (:);
clusterProps (:,2) = numberOfCells (:);
clusterProps (:,3) = clusterArea (:);
clusterProps (:,4) = clusterPerimeter (:);
clusterProps (:,5) = clusterPerimeterElements (:);
    
% One row of cellProps gives all information for one set of coordinates
cellProps = zeros (length (coord), 3);
cellProps(:,1) = coord(:,1);
cellProps(:,2) = coord(:,2);
cellProps(:,3) = clusterNr(:);
