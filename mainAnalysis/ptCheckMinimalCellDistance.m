function newCoord = ptCheckMinimalCellDistance (nucleiCoord, haloCoord, minDistCellCell)                                    
% ptCheckMinimalCellDistance combines two lists of coordinates and ensures a
% minimal distance between every combination of two cells
%
% SYNOPSIS       newCoord = ptCheckMinimalCellDistance(nucleiCoord,haloCoord,minDistCellCell)   
%
% INPUT          nucleiCoord     : a set of nuclei coordinates
%                haloCoord       : a set of halo coordinates
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

% Ensure a minimal distance between nuclei coordinates
iCount = 1;
while iCount < length (nucleiCoord)
   distance = [];
   distance = min (sqrt ((nucleiCoord (iCount+1:end, 1) - nucleiCoord (iCount, 1)).^2 + ...
                         (nucleiCoord (iCount+1:end, 2) - nucleiCoord (iCount, 2)).^2));

   % Test the distance between them
   if distance < minDistCellCell
      % Throw away the coordinate, because the distance is to small
      nucleiCoord(iCount,:) = [];
      iCount = iCount - 1;
   end
   iCount = iCount + 1;
end
clear iCount;    
clear distance;

% Ensure minimal distance between the halos
jCount = 1;
while jCount < length (haloCoord)
   distance = [];
   distance = min (sqrt ((haloCoord (jCount+1:end, 1) - haloCoord (jCount, 1)).^2 + ...
                         (haloCoord (jCount+1:end, 2) - haloCoord (jCount, 2)).^2));

   % Test the distance between them
   if distance < minDistCellCell
      % Throw away the coordinate, because the distance is to small
      haloCoord (jCount,:) = [];
      jCount = jCount - 1;
   end
  jCount = jCount + 1;
end
clear distance;
clear jCount;

% Ensure minimal distance between nuclei and halos
% If the distance is greater, add the halo coordinates to the mix. AK: why?????
haloCoordKept = [0,0];
if ~isempty (haloCoord)
   for hCount = 1 : size (haloCoord, 1)
      distance = [];
      distance = min (sqrt ((nucleiCoord (:, 1) - haloCoord (hCount, 1)).^2 + ...
                            (nucleiCoord (:, 2) - haloCoord (hCount, 2)).^2));

      % Note that here the minimal distance is larger than
      % between two cells found by the same routine, because... ahhhmmm... just to make sure
      if distance  >  (1.5 * minDistCellCell)
         haloCoordKept (end+1, 1) = haloCoord (hCount, 1);
         haloCoordKept (end, 2) = haloCoord (hCount, 2);
      end
   end

   % Throw away the [0,0] coordinates again
   haloCoordKept (1,:) = [];
     
   % If we found any new coordinates, cat them together with the nuclei ones
   if ~isempty (haloCoordKept)
      newCoord = cat (1, nucleiCoord, haloCoordKept);
   else
      newCoord = nucleiCoord;
   end
else
   newCoord = nucleiCoord;
end

newCoord = nucleiCoord;

clear haloCoord;    
clear namesnumbers;   
clear nucleiCoord;   
clear distance;


