function [newM, lostCells] = ptCloseGapsInTracks (M, matchedLostCells, lostCells, frameNr, startFrame, increment, clusterDirectory)
% ptCloseGapsInTracks closes gaps in cell tracks based on lost cells that have been matched up with
% new cells
%
% SYNOPSIS       [newM, lostCells] = ptCloseGapsInTracks (M, matchedLostCells, lostCells, 
%                                                         frameNr, startFrame, increment, clusterDirectory)
%
% INPUT          M : M stack as returned by the tracker functions
%                        M = [y x y x]   [y x y x]   [y x y x]
%                               ...    ,    ...    ,    ...
%                             t1   t2     t2   t3     t3   t4
%                                1           2           3
%                matchedLostCells : matching coordinates in the format [y1 x1 y2 x2]
%                lostCells : a matrix of lost cell coordinates incl the framenr they 
%                            went lost: [y x frame]
%                frameNr : the frame where the new cells have been found
%                startFrame: the startFrame of the movie (needed for M position calculation)
%                increment: the rate of increment of the movie frames (needed for M position calculation)
%                clusterDirectory : the directory where the binary cluster mat files can be found
%
% OUTPUT         M : a new M stack where gaps have been closed using linear interpolation
%                lostCells : the updated matrix containing the cells that haven't been matched yet
%
% DEPENDENCIES   ptCloseGapsInTracks uses { nothing }
%
%                ptCloseGapsInTracks is used by { ptTrackCells }

% Get the lost cells that were matched up with new cells
matchedCells = matchedLostCells (find (matchedLostCells (:,3) ~= 0 & matchedLostCells (:,4) ~= 0), :);

% For each of these process the M matrix
for iCount = 1 : size (matchedCells,1)
   % Find the index of the matched lost cell in LostCells
   matchedCellIndex = find (lostCells (:,1) == matchedCells (iCount,1) & lostCells (:,2) == matchedCells (iCount,2));

   % This way we find out which frame that cell was lost
   frameLostCell = lostCells (matchedCellIndex, 3);
   
   % In case there are multiple entries (which means a lost cell exists
   % twice in the matrix) we only take the last one
   frameLostCell = frameLostCell (end,1);
   
   % Based on this info the position in M where the frame was lost can now
   % be calculated
   lostCellMEntry = ceil ((frameLostCell - startFrame) / increment);
   
   % Get the lost cell coordinates
   lostCell = matchedCells (iCount, 1:2);

   % Get the accompanying new cell coordinates as well
   newCell = matchedCells (iCount, 3:4);
   
   % Use all of this information to find the lost cell entry in M
   lostCellMInd = find (M (:,1,lostCellMEntry) == lostCell (1) & M (:,2,lostCellMEntry) == lostCell (2));

   % Same for the new cell entry in M
   newCellMInd = find (M (:,3,frameNr) == newCell (1) & M (:,4,frameNr) == newCell (2));

   % How many frames are in between the lost and new cell?
   numberOfFrames = frameNr - lostCellMEntry;

   % Using the lost and new cell coordinates and the number of frames in between, the
   % missing coordinates can be calculated
   y = []; x = [];
   for jCount = 1 : numberOfFrames
      y(jCount) = (floor(((newCell(1) - lostCell(1)) / numberOfFrames)) * jCount) + lostCell(1);
      x(jCount) = (floor(((newCell(2) - lostCell(2)) / numberOfFrames)) * jCount) + lostCell(2);
   
      % Using the cluster (binary) image of the cells (nuclei and halos combined) make sure
      % the calculated coordinate are actually in a cell area
      % First locate and load the correct binary cluster image
      cd (clusterDirectory);
      formatStr = sprintf ('%%.%dd', 3);
      imageNr = sprintf (formatStr, frameLostCell + jCount - 1);
      clusterFile = ['clusters' imageNr];
      load (clusterFile);

      % Check that the calculated coordinates fall into a cluster area (area with on-pixels)
      if clusterImage (x(end), y(end)) ~= 1
         % The track we are trying to close should not be closed since the coordinates
         % of the calculated cell do not fall into a cluster area
         break;         % Continue with the rest of the cells in lostCells
      end
      
      % Use these values to update M and this way close the track
      if (jCount == 1) && (numberOfFrames == 1)
         % Update the M entry (pos 3 and 4) where the cell was lost (lostCellMEntry)
         M (lostCellMInd,3,lostCellMEntry) = y(jCount);
         M (lostCellMInd,4,lostCellMEntry) = x(jCount);
         
         % Update pos 1 and 2 in the M entry where the new cell was found as well
         M (newCellMInd,1,frameNr) = y(jCount);
         M (newCellMInd,2,frameNr) = x(jCount);

      elseif (jCount == 1) && (numberOfFrames > 1)
         % Update the M entry (pos 3 and 4) where the cell was lost (lostCellMEntry)
         M (lostCellMInd,3,lostCellMEntry) = y(jCount);
         M (lostCellMInd,4,lostCellMEntry) = x(jCount);

         % Find the next zero entry so that it can be updated with the new values
         zeroInd = find (M (:,1,lostCellMEntry+jCount) == 0 & M (:,2,lostCellMEntry+jCount) == 0 & ...
                         M (:,3,lostCellMEntry+jCount) == 0 & M (:,4,lostCellMEntry+jCount) == 0);
         M (zeroInd(1),1,lostCellMEntry+jCount) = y(jCount);
         M (zeroInd(1),2,lostCellMEntry+jCount) = x(jCount);

      elseif (jCount > 1) && (jCount < numberOfFrames)
         % Find the entry (pos 1 and 2) that we updated for the previous jCount because we have
         % to update the second half of this as well (pos 3 and 4)
         entryInd = find (M (:,1,lostCellMEntry+jCount-1) == y(jCount-1) & M (:,2,lostCellMEntry+jCount-1) == x(jCount-1));
         M (entryInd,3,lostCellMEntry+jCount-1) = y(jCount);
         M (entryInd,4,lostCellMEntry+jCount-1) = x(jCount);

         % Find the next zero entry so that can be updated with the new values
         zeroInd = find (M (:,1,lostCellMEntry+jCount) == 0 & M (:,2,lostCellMEntry+jCount) == 0 & ...
                         M (:,3,lostCellMEntry+jCount) == 0 & M (:,4,lostCellMEntry+jCount) == 0);
         M (zeroInd(1),1,lostCellMEntry+jCount) = y(jCount);
         M (zeroInd(1),2,lostCellMEntry+jCount) = x(jCount);

      elseif (jCount > 1) && (jCount == numberOfFrames)    
         % Update the M entry before the entry where the new cell was found
         entryInd = find (M (:,1,lostCellMEntry+jCount-1) == y(jCount-1) & M (:,2,lostCellMEntry+jCount-1) == x(jCount-1));
         M (entryInd,3,frameNr-1) = y(jCount);
         M (entryInd,4,frameNr-1) = x(jCount);

         M (newCellMInd,1,frameNr) = y(jCount);
         M (newCellMInd,2,frameNr) = x(jCount);
      end
   end
   
   % Remove the processed entry from lostCells because we don't want to
   % process it again later on
   lostCells (matchedCellIndex,:) = [];
   
end    % for iCount = 1 : size (matchedCells,1)

% Return the new M matrix
newM = M;
