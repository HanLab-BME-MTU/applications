function newCoord = ptCheckMinimalCellDistance (nucleiCoord, extraNucCoord, minDistCellCell)
% ptCheckMinimalCellDistance combines two lists of coordinates and ensures a
% minimal distance between every combination of two cells or in case the
% second list of coords is empty it ensures that there are no cells in the
% first set that are too close together
%
% SYNOPSIS       newCoord = ptCheckMinimalCellDistance(nucleiCoord,haloCoord,minDistCellCell)   
%
% INPUT          nucleiCoord     : a set of nuclei coordinates
%                extraNucCoord   : a set of possible nuclei coordinates.
%                This can be [] in case only the nucleiCoord has te be tested
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

% Ensure minimal distance between the would be nuclei and the real ones if
% there are any (extraNucCoord will be [] if no vector is provided)
if ~isempty (extraNucCoord)
   jCount = 1;
   while jCount < length (extraNucCoord)
      distance = [];
      distance = min (sqrt ((nucleiCoord (jCount+1:end, 1) - extraNucCoord (jCount, 1)).^2 + ...
                            (nucleiCoord (jCount+1:end, 2) - extraNucCoord (jCount, 2)).^2));

      % Test the distance between them
      if distance < minDistCellCell
         % Throw away the coordinate, because the distance is to small
         extraNucCoord (jCount,:) = [];
         jCount = jCount - 1;
      end
      jCount = jCount + 1; 
   end
end

% If we found any new nuclei coordinates, cat them together with the nuclei ones
if ~isempty (extraNucCoord)
   newCoord = cat (1, nucleiCoord, extraNucCoord);
else
   newCoord = nucleiCoord;
end

% Just in case there are duplicate entries remove them now
newCoord = unique (newCoord, 'rows');


