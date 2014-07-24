function newCoord = ptSubtractCellArea (imgIn, coord, cellMaxSize)
% ptSubtractCellArea subtracts a round cell area given by the center coordinates and 
% parameter from a binary image.
%
% SYNOPSIS       imgOut = ptSubtractCellArea (imgIn, coord, cellMaxSize)
%
% INPUT          imgIn     : a binary image where cell area has value 1 
%                coord     : the center coordinates of the circle to subtract
%                cellMaxSize: the maximum size a cell can have 
%
% OUTPUT         newCoord : coordinates found in the binary image after subtracting the old set
%                           of coordinates
%
% DEPENDENCIES   ptSubtractCellArea uses { nothing }
%                                  
%                ptSubtractCellArea is used by { ptCalculateCellArea }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        May 04          Initial Version

% Calculate the diameter of a circle with area maxSizeNucleus
maxDiameterCell = round (2 * sqrt ((cellMaxSize*2) / pi));

% Make sure this diameter is an odd number
if rem (maxDiameterCell, 2) == 0
   maxDiameterCell = maxDiameterCell + 1;
end

% Create a grid to get the circle mask
radius = round ((maxDiameterCell - 1) / 2);
[xx yy] = ndgrid ([-1 * radius : radius], [-1 * radius : radius]);

% From this grid create a distance matrix
dist = sqrt (xx.^2 + yy.^2);

% Get the size of the input binary image
[hIn, wIn] = size (imgIn);

% Loop through all the coordinates in coord
for iCount = 1 : size (coord, 1)
    
   % Initialize the circle mask
   circleMask = ones (maxDiameterCell);
  
   % And set the pixels in the circle area to 1
   circleMask (find (dist <= radius)) = 0;

   % Determine the min and max values of the area we want to subtract from
   xmin = round (coord (iCount, 1) - radius);
   xmax = round (coord (iCount, 1) + radius);
   ymin = round (coord (iCount, 2) - radius);
   ymax = round (coord (iCount, 2) + radius);

   % Make sure these values fall within the boundaries of the image
   if xmin < 1
      circleMask (:, 1:abs(xmin)+1) = [];
      xmin = 1;
   end

   if xmax > wIn
      circleMask (:, end-(abs(xmax)-wIn-1):end) = [];
      xmax = wIn;
   end

   if ymin < 1
      circleMask (1:abs(ymin)+1, :) = [];
      ymin = 1;
   end

   if ymax > hIn
      circleMask (end-(abs(ymax)-hIn-1):end, :) = [];
      ymax = hIn;
   end

   % Subtract the circle region from the binary in image
   imgIn (ymin:ymax, xmin:xmax) = imgIn (ymin:ymax, xmin:xmax) & circleMask;
   
end    % for loop

%figure,imshow(imgIn);
%hold on;

% Use the remainder of the input image to find any new coordinates
% First do some morphological operations
imgIn = imerode (imgIn, strel ('disk', 2));
imgIn = imdilate (imgIn, strel ('disk', 2));
imgIn = imopen (imgIn, strel ('disk', 3));
imgIn = imclose (imgIn, strel ('disk', 3));

% Label the image
imgInLabeled = bwlabel (imgIn);

% Calculate cell properties using the regionprops function from the image toolbox
cellProperties = regionprops (imgInLabeled, 'Area', 'PixelList');

% Find regions smaller than a cell so that we can remove these
cellArea = [cellProperties.Area];
cellRemoveArea = find (cellArea < cellMaxSize*2);
tempProperties = cellProperties (cellRemoveArea);

% Get all pixels within those regions
removeCoord = cat (1, tempProperties(:).PixelList);

% Set all of those pixels to zero so that we can forget about them
for iCount = 1 : size (removeCoord, 1)
   imgIn (removeCoord(iCount, 2), removeCoord(iCount, 1)) = 0;
end

% Now label the new image again
imgInLabeled = bwlabel (imgIn);

% Calculate cell properties using the regionprops function again
cellProperties = regionprops (imgInLabeled, 'Centroid');

% And get the coordinates of all of these
newCoord = round (cat (1, cellProperties.Centroid));

%plot (newCoord (:,1), newCoord (:,2), 'r.');
%plot (coord (:,1), coord (:,2), 'g.');
%hold off;

