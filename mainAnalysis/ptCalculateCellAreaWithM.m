function [cellProps, clusterProps] = ptCalculateCellAreaWithM (M, distanceToCellArea, minSizeNucleus, clusterDirectory, startFrame, endFrame, increment)
% ptCalculateCellArea  determines what areas of an image are occupied by cells and
%                      calculates image properties
%
% SYNOPSIS     [cellProps, clusterProps] = ptCalculateCellAreaWithM (M, distanceToCellArea, minSizeNucleus,
%                                                                    clusterDirectory)
%
% INPUT    M : the magic position matrix
%          distanceToCellArea : distance a set of coordinates may have to an cell area and still belong to it
%	       minSizeNucleus : everything smaller than this size will be removed from the binary image
%          clusterDirectory : place where cluster images are stored
%          startFrame : frame the movie starts with (not necessarily 1)
%          endFrame : frame the movie ends with (not necessarily the last frame)
%          increment : what's the increment between frames
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

% Let the user know we're starting to calculate cell and cluster props
fprintf (1, '\n     Generating cell and cluster properties...\n');
fprintf (1, '     Processing frame: ');

% Use the size of the M matrix to determine how many frames we should process
nrOfFrames = size (M, 3) + 1;

% Initialize loop counter for M index
mCount = 0;

% Initialize empty cell and cluster rows
emptyCell            = zeros (1,3);
emptyCluster         = zeros (1,5);

% Initialize the coordinate matrix
coord = [];

for frameCount = startFrame : increment : endFrame

   % Where are we processing the movie?
   fprintf (1, '%d ', frameCount);
   
   % Increase M counter
   mCount = mCount + 1;
   
   % Load the cluster (binary) image of the cells (nuclei and halos combined)
   cd (clusterDirectory);
   formatStr = sprintf ('%%.%dd', 3);
   imageNr = sprintf (formatStr, frameCount);
   clusterFile = ['clusters' imageNr];
   load (clusterFile);
   
   % Label the cluster image
   imgLabeledCellArea = bwlabel (clusterImage);
   % Get the coordinates for this frame out of the M matrix. Since the
   % number of M entries is one less than the number of frames, the last
   % frame coordinates should be fetched differently.l
   if mCount < endFrame
      coord = M (:,1:2,mCount);
   else
      coord = M (:,3:4,mCount-1);
   end
    
   % Throw out the zeros that were likely to be in M
   coord (find (coord (:,1) == 0 & coord (:,2) == 0),:) = [];
    
   % Prepare a matrix for the grouping of coordinates to objects
   clusterNr = zeros (length (coord), 1);

   % Determine to which group each set of coordinates belongs
   for iCount = 1 : length (coord)          % process all rows in tempCoord

      % Do the coordinates fall in a labeled area at all?
      if imgLabeledCellArea (coord(iCount,2), coord(iCount,1)) ~= 0
         % If the coordinates fall into a labeled area, that cluster label is stored in a matrix
         clusterNr(iCount) = imgLabeledCellArea (coord(iCount,2), coord(iCount,1));
      else
         % We also take cells into account that are near enough a cell area
     
         % Get the size of the input image
         [img_h, img_w] = size (imgLabeledCellArea);

         % Calculate the coordinates of the point before the cell area
         x_1 = round (coord (iCount,2) - distanceToCellArea);
         y_1 = round (coord (iCount,1) - distanceToCellArea);

         % Calculate the coordinates of the point after the cell area
         x_2 = round (coord (iCount,2) + distanceToCellArea);
         y_2 = round (coord (iCount,1) + distanceToCellArea);
                   
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
      end   % if imgLabeledCellArea
   end   % for iCount = 1

   % Find all the cell coordinates that are part of the background (aka not of an area of cells)
   % and mark these empty; also unset the coordinates
   onBackGround = find (clusterNr == 0);
   if ~isempty (onBackGround)
      coord (onBackGround,:) = [];
      clusterNr (onBackGround) = [];
   end
 
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
      cluster = imgLabeledCellArea == uniqClusterNr (iCount);

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
   clusterProp = zeros (length (uniqClusterNr), 5);
   clusterProp (:,1) = uniqClusterNr (:);
   clusterProp (:,2) = numberOfCells (:);
   clusterProp (:,3) = clusterArea (:);
   clusterProp (:,4) = clusterPerimeter (:);
   clusterProp (:,5) = clusterPerimeterElements (:);
    
   % One row of cellProps gives all information for one set of coordinates
   cellProp = zeros (length (coord), 3);
   cellProp(:,1) = coord(:,1);
   cellProp(:,2) = coord(:,2);
   cellProp(:,3) = clusterNr(:);
 
   % Accumulate the cell properties
   cellProps (1 : size (emptyCell, 1), 1 : size (emptyCell, 2), mCount) = emptyCell;
   cellProps (1 : size (cellProp, 1), 1 : size (cellProp, 2), mCount) = cellProp;
   
   % Accumulate the cluster properties
   clusterProps (1 : size (emptyCluster, 1), 1 : size (emptyCluster, 2), mCount) = emptyCluster;
   clusterProps (1 : size (clusterProp, 1), 1 : size (clusterProp, 2), mCount) = clusterProp;
   
   % clear temporary variables
   clear onBackGround;
   clear cluster; clear clusterArea; clear clusterPerimeter; clear clusterInv;
   clear clusterInvLabel; clear clusterPerimeterElements;
   clear clusterProp; clear cellProp;
end

% Let the user know we've finished
fprintf (1, '\n     Finished generating cell and cluster properties.\n');
