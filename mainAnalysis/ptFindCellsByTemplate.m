function [newCoord, maxCorrelation] = ptFindCellsByTemplate (currentCoord, prevImage, currentImage, backgroundLevel,                          percentBackground, templateSize, searchAreaSize)
% ptFindCellsByTemplate finds cells in a subsequent frame using template matching (correlation)
%
% SYNOPSIS       [newCoord,maxCorrelation] = ptFindCellsByTemplate (currentCoord, prevImage, currentImage, backgroundLevel,
%                                                                                                    percentBackground, templateSize, searchAreaSize)
%
% INPUT 	 currentCoord          : coordinates of the cell we wish to find in the new image
% 		 prevImage          : the image in which that cell is
% 		 currentImage              : the image in which we wish to find the cell
% 		 backgroundLevel       : approximated level of the background
% 		 percentBackground     : threshold for ignoring areas within template and image
%                                        (1 +- percentBackground * backgroundLevel = thresh)
% 		 templateSize          : the size we want the template to have
% 		 searchAreaSize        : the size of the search area
%
% OUTPUT         newCoord       : found coordinates
%                maxCorrelation : maximum correlation
%
% DEPENDENCIES   ptFindCellsByTemplate uses { nothing }
%
%                ptFindCellsByTemplate is used by { ptTrackCells }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Determine the min and max values of the area we want to take as the template
xminTemplate = round (currentCoord (1, 1) - (templateSize - 1) / 2);
xmaxTemplate = round (currentCoord (1, 1) + (templateSize - 1) / 2);
yminTemplate = round (currentCoord (1, 2) - (templateSize - 1) / 2);
ymaxTemplate = round (currentCoord (1, 2) + (templateSize - 1) / 2);

% Get the size of the previous image
[prevImage_h, prevImage_w] = size (prevImage);

% Make sure these values fall within the boundaries of the image
if xminTemplate < 1
   xminTemplate = 1;
end

if xmaxTemplate > prevImage_w;
   xmaxTemplate = prevImage_w;
end

if yminTemplate < 1
   yminTemplate = 1;
end

if ymaxTemplate > prevImage_h
   ymaxTemplate = prevImage_h;
end

% Extract the template from the previous image
cellTemplate = prevImage (yminTemplate:ymaxTemplate, xminTemplate:xmaxTemplate);

% Now we replace the background with random noise, so that we don't correlate background to background
up = 1 + percentBackground;
down = 1 - percentBackground;
[backgroundRowIndex, backgroundColIndex] = find (cellTemplate < (backgroundLevel * up) & ...
                                                 cellTemplate > (backgroundLevel * down));

% replace the indexes of interest with a random noise value in between 0 and 1
for iCount = 1 : length (backgroundRowIndex)
   cellTemplate (backgroundRowIndex (iCount), backgroundColIndex (iCount)) = rand (1);
end

% Now we're starting the processing using the current image
% Determine the min and max values of the area we want to search in
xminSearch = round (currentCoord (1,1) - (searchAreaSize - 1) / 2);
xmaxSearch = round (currentCoord (1,1) + (searchAreaSize - 1) / 2);
yminSearch = round (currentCoord (1,2) - (searchAreaSize - 1) / 2);
ymaxSearch = round (currentCoord (1,2) + (searchAreaSize - 1) / 2);

% Get the size of the current image
[currentImage_h, currentImage_w] = size (currentImage);

% Make sure the min max values fall within the boundaries of the image
if xminSearch < 1
   xmaxSearch = xmaxSearch - xminSearch + 1;
   xminSearch = 1;
end

if xmaxSearch > currentImage_w
   xminSearch = xminSearch - xmaxSearch + currentImage_w;
   xmaxSearch = currentImage_w;
end

if yminSearch < 1
   ymaxSearch = ymaxSearch - yminSearch + 1;
   yminSearch = 1;
end

if ymaxSearch > currentImage_h
   yminSearch = yminSearch - ymaxSearch + currentImage_h;
   ymaxSearch = currentImage_h;
end

% Extract the search area part from the new image
searchArea = currentImage (yminSearch:ymaxSearch, xminSearch:xmaxSearch);

[backgroundRowIndex, backgroundColIndex] = find (searchArea < (backgroundLevel * up) & ...
                                                                               searchArea > (backgroundLevel * down));

% Replace the indexes of interest with a random noise value between 0 and 1
for jCount = 1 : length (backgroundRowIndex)
   searchArea (backgroundRowIndex (jCount), backgroundColIndex (jCount)) = rand(1);
end

% Perform a normalized crosscorrelation of the template in the search area
crossCorrelation = normxcorr2 (cellTemplate, searchArea);
crossCorSize = size (crossCorrelation);

% now we superimpose a distance criteria. First we have to build a matrix of
% the same size as crossCorrelation, then adjust it's values according to the
% distance from where the cell was in the last picture
distance = zeros (crossCorSize (1), crossCorSize (2));

for rowCount = 1 : crossCorSize (1,1)
   for colCount = 1 : crossCorSize (1,2)
      % Calculate the distance to the initial coordinates
      calcDist = sqrt ((colCount + xminSearch - round ((templateSize - 1) / 2) - currentCoord (1,1))^2 + ...
                              (rowCount + yminSearch - round ((templateSize - 1) / 2) - currentCoord (1,2))^2);
      if calcDist < 20
         % If the distance to the initial coordinates is smaller than 20 we make it 1
         distance (rowCount, colCount) = 1;
      else
         % It is larger than 20 thus the value slowly decreases
         distance (rowCount, colCount) = 1 - ((calcDist - 20) / 80);
      end

      % Negative values for the distance are not allowed so we make them 0
      if distance (rowCount, colCount) < 0
         distance (rowCount, colCount) = 0;
      end
   end
end

% Multiply crossCorrelation with distance, thus having introduced the distance criteria
distCriteria = crossCorrelation .* distance;

clear distance;

% Find the coordinates of the absolute intensity maxima in the image
[maxDistCrit distCritIndex] = max (distCriteria);
[maxCorrelation correlationIndex] = max (maxDistCrit);
lastIndex = distCritIndex (correlationIndex);

% Add together everything that's necessary to receive the coordinates of the
% found spot in the frame (not within crossCorrelation)
newCoord (1,1) = correlationIndex + xminSearch - round (size (cellTemplate, 1) / 2);
newCoord (1,2) = lastIndex + yminSearch - round (size (cellTemplate, 2) / 2);

% Again we have to make sure the new coordinates will fall within image boundaries
if newCoord (1,1) < 1
   newCoord (1,1) = 1;
end

if newCoord (1,1) > (currentImage_w - 1)
   newCoord (1,1) = currentImage_w - 1;
end

if newCoord (1,2) < 1
   newCoord (1,2) = 1;
end

if newCoord (1,2) > (currentImage_h - 1)
   newCoord (1,2) = currentImage_h - 1;
end

% Mark the new coordinates (using the fraction) as found by template
%newCoord (1,1) = newCoord (1,1) + oldCellTemplateMarker;
%newCoord (1,2) = newCoord (1,2) + oldCellTemplateMarker;

% Mark cells that were new cells as new cells propagated by templates
%realCurCoord = floor (currentCoord + 0.1);
%marker = round (10 * (currentCoord - realCurCoord));

% If the newly found coordinates belong to a cell found by template matching, mark it appropriately
% (which is done by setting the fraction to the correct marker value)
%if marker == newCell | marker == newCellTemplate
%  newCoord = floor (newCoord) + newCellTemplateMarker;
%end
