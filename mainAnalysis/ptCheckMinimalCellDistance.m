function newCoord = ptCheckMinimalCellDistance (inputCoord, newInputCoord, minDistCellCell)
% ptCheckMinimalCellDistance combines two lists of coordinates and ensures a
% minimal distance between every combination of two cells or in case the
% second list of coords is empty it ensures that there are no cells in the
% first set that are too close together
%
% SYNOPSIS       newCoord = ptCheckMinimalCellDistance (inputCoord, newInputCoord, minDistCellCell)   
%
% INPUT          inputCoord     : a set of nuclei coordinates
%                newInputCoord   : a set of possible nuclei coordinates.
%                This can be [] in case only the inputCoord has te be tested
%                minDistCellCell : minimal distance between two cells (in pixels)
%
% OUTPUT         newCoord : the combined coordinates, with minimal distance between them
%
% DEPENDENCIES   ptCheckMinimalCellDistance uses { nothing }
%                                  
%                ptCheckMinimalCellDistance is used by { ptTrackCells
%                                                        ptInitializeJob }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Apr 04          Cleaned up source
% Andre Kerstens        Jun 04          Changed variable names to something
%                                       more general than nucleicoord

% Ensure a minimal distance between nuclei coordinates
iCount = 1;
while iCount < length (inputCoord)
   distance = [];
   distance = min (sqrt ((inputCoord (iCount+1:end, 1) - inputCoord (iCount, 1)).^2 + ...
                         (inputCoord (iCount+1:end, 2) - inputCoord (iCount, 2)).^2));

   % Test the distance between them
   if distance < minDistCellCell
      % Throw away the coordinate, because the distance is to small
      inputCoord(iCount,:) = [];
      iCount = iCount - 1;
   end
   iCount = iCount + 1;
end
clear iCount;    
clear distance;

% Ensure minimal distance between the would be nuclei and the real ones if
% there are any (newInputCoord will be [] if no vector is provided)
if ~isempty (newInputCoord)
   jCount = 1;
   while jCount < length (newInputCoord)
      distance = [];
      distance = min (sqrt ((inputCoord (jCount+1:end, 1) - newInputCoord (jCount, 1)).^2 + ...
                            (inputCoord (jCount+1:end, 2) - newInputCoord (jCount, 2)).^2));

      % Test the distance between them
      if distance < minDistCellCell
         % Throw away the coordinate, because the distance is to small
         newInputCoord (jCount,:) = [];
         jCount = jCount - 1;
      end
      jCount = jCount + 1; 
   end
end

% If we found any new nuclei coordinates, cat them together with the nuclei ones
if ~isempty (newInputCoord)
   newCoord = cat (1, inputCoord, newInputCoord);
else
   newCoord = inputCoord;
end

% Just in case there are duplicate entries remove them now
newCoord = unique (newCoord, 'rows');


