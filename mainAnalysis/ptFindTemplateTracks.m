function [newCoord, maxCorrelation] = ptFindTemplateTracks (currentCoord, currentImage, newImage, backgroundLevel, oldCellTemplateMarker, newCell, newCellTemplate, newCellTemplateMarker, percentBackground, templateSize, searchAreaSize)
% ptFindTemplateTracks finds cells in a subsequent frame using template matching (correlation)
%
% SYNOPSIS       [newCoord,maxCorrelation] = ptFindTemplateTracks (currentCoord, currentImage, newImage, backgroundLevel, 
%                                            oldCellTemplateMarker, newCell, newCellTemplate, newCellTemplateMarker, 
%                                            percentBackground, templateSize, searchAreaSize)
%
% INPUT 	 currentCoord          : coordinates of the cell we wish to find in the new image
% 		 currentImage          : the image in which that cell is 
% 		 newImage              : the image in which we wish to find the cell 
% 		 backgroundLevel       : approximated level of the background
% 		 oldCellTemplateMarker : marker for old cells found by templates
% 		 newCell               : identifier of new cells 
% 		 newCellTemplate       : identifier of new cells found by template
% 		 newCellTemplateMarker : marker for new cells found by templates
% 		 percentBackground     : threshold for ignoring areas within template and image 
%                                        (1 +- percentBackground * backgroundLevel = thresh)
% 		 templateSize          : the size we want the template to have
% 		 searchAreaSize        : the size of the search area
% 
% OUTPUT         newCoord       : found coordinates
%                maxCorrelation : maximum correlation
%
% DEPENDENCIES   ptFindTemplateTracks uses { nothing }
%                                  
%                ptFindTemplateTracks is used by { ptTrackCells }
%                                   
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Determine the min and max values of the area we want to take as the template
xminCur = round (currentCoord (1, 1) - (templateSize - 1) / 2);
xmaxCur = round (currentCoord (1, 1) + (templateSize - 1) / 2);
yminCur = round (currentCoord (1, 2) - (templateSize - 1) / 2);
ymaxCur = round (currentCoord (1, 2) + (templateSize - 1) / 2);

% Get the size of the current image
[curImage_h, curImage_w] = size (currentImage);

% Make sure these values fall within the boundaries of the image
if xminCur < 1
   xminCur = 1;
end

if xmaxCur > curImage_w;
   xmaxCur = curImage_w;
end

if yminCur < 1
   yminCur = 1;
end

if ymaxCur > curImage_h
   ymaxCur = curImage_h;
end
   
% Extract the template from the current image
cellTemplate = currentImage (yminCur:ymaxCur, xminCur:xmaxCur);

% Now we replace the background with random noise, so that we don't correlate background to background
up = 1 + percentBackground;
down = 1 - percentBackground;
[backgroundRowIndex, backgroundColIndex] = find (cellTemplate < (backgroundLevel * up) & ...
                                                 cellTemplate > (backgroundLevel * down));
% replace the indexes of interest with random noise (0 to 255)
for iCount = 1 : length (backgroundRowIndex)
   cellTemplate (backgroundRowIndex (iCount), backgroundColIndex (iCount)) = rand (1) * 255;
end

% Clear these values for later reuse
clear backgroundRowIndex;
clear backgroundColIndex;


% Here we're starting on the new image
% Determine the min and max values of the area we want to search in
xminNew = round (currentCoord (1,1) - (searchAreaSize - 1) / 2);
xmaxNew = round (currentCoord (1,1) + (searchAreaSize - 1) / 2);
yminNew = round (currentCoord (1,2) - (searchAreaSize - 1) / 2);
ymaxNew = round (currentCoord (1,2) + (searchAreaSize - 1) / 2);
   
% Get the size of the new image
[newImage_h, newImage_w] = size (newImage);

% Make sure the min max values fall within the boundaries of the image
if xminNew < 1
   xmaxNew = xmaxNew - xminNew + 1;
   xminNew = 1;
end

if xmaxNew > newImage_w
   xminNew = xminNew - xmaxNew + newImage_w;
   xmaxNew = newImage_w;    
end

if yminNew < 1
   ymaxNew = ymaxNew - yminNew + 1;
   yminNew = 1;
end

if ymaxNew > newImage_h
   yminNew = yminNew - ymaxNew + newImage_h;
   ymaxNew = newImage_h;
end

% Extract the search area part from the new image
searchArea = newImage (yminNew:ymaxNew, xminNew:xmaxNew);

[backgroundRowIndex, backgroundColIndex] = find (searchArea < (backgroundLevel * 1.2) & ...
                                                 searchArea > (backgroundLevel * 0.8));

% replace the indexes of interest with random noise (0 to 1) 
for jCount = 1 : length (backgroundRowIndex)
   searchArea (backgroundRowIndex (jCount), backgroundColIndex (jCount)) = rand(1);
end

% Perform a normalized crosscorrelation of the template in the search area
crossCorValue = normxcorr2 (cellTemplate, searchArea);
crossCorSize = size (crossCorValue);

% now we superimpose a distance criteria. First we have to build a matrix of
% the same size as crossCorValue, then adjust it's values according to the
% distance from where the cell was in the last picture
distance = zeros (crossCorSize (1,1), crossCorSize (1,2));

for rowCount = 1 : crossCorSize (1,1)
   for colCount = 1 : crossCorSize (1,2)
      % If the distance to the initial coordinates is smaller than 20, the value is 1
      if (sqrt ((colCount + xminNew - round ((templateSize - 1) / 2) - currentCoord (1,1))^2 + ...
                (rowCount + yminNew - round ((templateSize - 1) / 2) - currentCoord (1,2))^2)) < 20

         distance (rowCount, colCount) = 1;
            
      % If it is larger than 20, the value slowly decreases
      else 
         distance (rowCount, colCount) = 1 - ((sqrt ((colCount + xminNew - round ((templateSize - 1) / 2) - currentCoord (1,1))^2 + ...
                                                 (rowCount + yminNew - round ((templateSize - 1) / 2) - currentCoord (1,2))^2) - ...
                                                 20) / 80);
      end
      
      % Negative values become zero
      if distance (rowCount, colCount) < 0
         distance (rowCount, colCount) = 0;
      end
   end
end

% Multiply crossCorValue with distance, thus having introduced the distance criteria
distCriteria = crossCorValue .* distance;

clear distance;

% Find the coordinates of the absolute intensity maxima in the image
[maxDistCrit distCritIndex] = max (distCriteria);
[maxCorrelation correlationIndex] = max (maxDistCrit);
lastIndex = distCritIndex (correlationIndex);

% Add together everything that's necessary to receive the coordinates of the
% found spot in the frame (not within crossCorValue)
newCoord (1,1) = correlationIndex + xminNew - round (size (cellTemplate, 1) / 2);
newCoord (1,2) = lastIndex + yminNew - round (size (cellTemplate, 2) / 2);

% Again we have to make sure the new coordinates will fall within image boundaries
if newCoord (1,1) < 1
   newCoord (1,1) = 1;
end

if newCoord (1,1) > (newImage_w - 1)
   newCoord (1,1) = newImage_w - 1;
end

if newCoord (1,2) < 1
   newCoord (1,2) = 1;
end

if newCoord (1,2) > (newImage_h - 1)
   newCoord (1,2) = newImage_h - 1;
end

% Mark the new coordinates (using the fraction) as found by template
newCoord (1,1) = newCoord (1,1) + oldCellTemplateMarker;
newCoord (1,2) = newCoord (1,2) + oldCellTemplateMarker;

% Mark cells that were new cells as new cells propagated by templates
realCurCoord = floor (currentCoord + 0.1);
identifier = round (10 * (currentCoord - realCurCoord));

% If the newly found coordinates belong to a cell found by template matching, mark it appropriately
% (which is done by setting the fraction to the correct marker value)
if identifier == newCell | identifier == newCellTemplate
   newCoord = floor (newCoord) + newCellTemplateMarker;
end
