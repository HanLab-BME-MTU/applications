function [newCoord, imgClusterArea, labeledClusterImage] = ptFindCoordFromClusters (inputImage, coord, minSizeNucleus)
% ptFindCoordFromClusters tries to find unfound single cells based on the average cell size
%       and clusters without any coordinates
%
% SYNOPSIS     [newCoord] = ptFindCoordFromClusters (labeledClusterImage, coord, avgCellArea)
%
% INPUT    inputImage : segmented image (1 = nuclei, 2 = background, 3 = halos)
%          coord : set of coordinates
%          minSizeNucleus : the minimum size a nucleus should have
%
% OUTPUT   newCoord : newly found coordinates
%          imgClusterImage: a binary image containing all of the clusters
%          labeledClusterImage: an image where all the clusters have unique numbers
%
% DEPENDENCIES ptFindCoordFromClusters uses {nothing}
%
%              ptFindCoordFromClusters is used by { ptTrackCells }
%
% CHANGE REVISION
%
% Name                  Date            Comment
% --------------------- ----------      -----------------------------------------------
% Andre Kerstens        May 04          Wrote a function to determine coordinates based
%                                       on average cell size

% The binary image is an or-function of the input image where pixels are 1 (nuclei) and the
% input image where pixels are 3 (halos)
imgClusterArea = (inputImage == 1 | inputImage == 3);

% Fill out any gaps by dilating and eroding the image a couple of times.
imgClusterArea = imdilate (imgClusterArea, strel ('disk', 8));
imgClusterArea = imerode  (imgClusterArea, strel ('disk', 8));
imgClusterArea = imdilate (imgClusterArea, strel ('disk', 4));
imgClusterArea = imerode  (imgClusterArea, strel ('disk', 3));
imgClusterArea = imdilate (imgClusterArea, strel ('disk', 2));
imgClusterArea = imerode  (imgClusterArea, strel ('disk', 1));
imgClusterArea = imdilate (imgClusterArea, strel ('disk', 1));

% Get rid of all the really small objects (smaller than a nucleus)
imgClusterArea = bwareaopen (imgClusterArea, minSizeNucleus);

% Label the objects in imgCellArea: 0 = background, 1 = first object, 2 = second object, etc
labeledClusterImage = bwlabel (imgClusterArea);

% Find the number of clusters in the image
maxClusterNr = max (max (labeledClusterImage));

% Get the area and centroids of all the clusters in the image
clusterProps = regionprops (labeledClusterImage, 'Centroid', 'Area');

% Initialize matrices to store cell area and clusters with no coords
cellArea = [];
clusterNoCoords = [];

% For every cluster determine whether it contains coordinates from coord
for iCount = 1 : maxClusterNr
   % Find the coordinates of the pixels in the cluster
   [pixelCoordY, pixelCoordX] = find (labeledClusterImage == iCount);

   % Match these with the coordinates in coord
   matchedCoords = [];
   for jCount = 1 : size (coord,1)
      coordIndex = find (pixelCoordY (:) == coord (jCount,2) & pixelCoordX (:) == coord (jCount,1));
      if ~isempty (coordIndex)
         matchedCoords(end+1) = coordIndex;
      end
   end

   % If we didn't find any, it means that this cluster doesn't contain any previously 
   % found nuclei so store its number
   if ~isempty (matchedCoords)
      % Calculate the average cell area using the cluster size and the number of cells
      % we found in it
      cellArea(end+1) = clusterProps(iCount).Area / length (matchedCoords);
   else
      % There are no coordinates in this cluster so keep it for later
      clusterNoCoords (end+1) = iCount;
   end
end 

% Determine the average cell area using the data of all the clusters
avgCellArea = round (sum (cellArea) / length (cellArea));

% Now that we have the average cell size compare it with all the coordinate-less
% clusters to see if it could be a cell
% Initialize the matrix for the new coordinates
newCoord = [];

% For all cells without any coordinates:
% Get the area and centroid of this cluster if it has the size of a single
% cell (size between half and twice the size of an average cell)
for iCount = 1 : length (clusterNoCoords)
   if (clusterProps(clusterNoCoords(iCount)).Area >= avgCellArea/2) & ...
      (clusterProps(clusterNoCoords(iCount)).Area < avgCellArea*2)
      % We found a new cell which wasn't detected before: store it
      newCoord (end+1,1) = round(clusterProps(clusterNoCoords(iCount)).Centroid(1));
      newCoord (end,2) = round(clusterProps(clusterNoCoords(iCount)).Centroid(2));
   end
end
