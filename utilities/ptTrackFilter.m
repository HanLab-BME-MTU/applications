function MPM = ptTrackFilter (MPM, plusFrames, minusFrames, maxRelinkDist, minTrackLength, savePath)
% ptTrackFilter filters and relinks the tracks given in MPM
%
% SYNOPSIS       MPM = ptTrackFilter (MPM, plusFrames, minusFrames, maxRelinkDist, minTrackLength, savePath)
%
% INPUT          MPM        : Magic Position Matrix 
%                               MPM = [ y  x  y  x  y  x ... ]
%                                        t1    t2    t3
%                plusFrames : how many frames after the currently analysed
%                             track stops must a candidate for gap closing survive
% 		  		 minusFrames : how many frames (at most!) before the currently analysed
%                              track stops can a candidate for gap closing be present
% 				 maxRelinkDist : max distance for relinking
% 				 minTrackLength : minimal length of tracks (or else they get erased)
% 				 savePath : where shall the new MPM be saved
%
%
% OUTPUT         MPM : the modified MPM matrix
%
% DEPENDENCIES   ptTrackFilter  uses {nothing}
%                                  
%                ptTrackFilter is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04

% Write a message to the screen
fprintf (1, 'Automatically filtering and relinking tracks using provided values...\n');

% If there are totally empty rows in MPM, erase them
% First find all unique rows unequal to zero
[notZeroEntryRows, notZeroEntryCols] = find (MPM);
notZeroEntryRows = unique (notZeroEntryRows);

% Get these entries from MPM and throw away the rest
MPM = MPM (notZeroEntryRows,:);

% First we find (for every cell/track) the frame in which it appears and the frame
% in which it disappears
for iCount = 1 : size (MPM,1)
    
   % Fetch one row of track coordinates from MPM
   track = MPM (iCount,:);
   
   % Find out where in this track coordinates can be found
   rowIndex = find (track);
   
   % Retrieve the start frame
   trackStart (iCount) = min (rowIndex);
   
   % Retrieve the ending frame
   trackEnd (iCount) = max (rowIndex);
end

% Initialize track counter
trackCount = 0;

% Loop over the number of tracks. Note: the older the track is, the lower 
% it's row index. So we always try to link tracks with higher
% row index (newer tracks) to ones with lower row index (older tracks) and 
% subsequently erase the newer track. In this way the likelyhood of assigning 
% a new number to a track already found is reduced
while trackCount < (size (MPM,1) - 0.3)
    
   % Start out with no links made 
   linked = 0;
    
   % Increase track counter
   trackCount = trackCount + 1;
    
   % Look for tracks which begin only a few frames before our current
   % track stops and try to link 
   neighbours = find ( (trackEnd (trackCount) - trackStart (trackCount + 1:end)) < minusFrames & ...
                      (trackEnd (trackCount) - trackStart (trackCount + 1:end)) > 0 & ...
                      (trackEnd (trackCount) - trackEnd (trackCount + 1:end)) < plusFrames & ...
                      (trackStart (trackCount + 1:end)) > 1.1 & ...
                      (trackEnd (trackCount) + 1.9) < size (MPM,2) );
    
   % In case we find any neighbours tracks that are in the specified range do:
   if ~isempty (neighbours)
       
      % Calculate the distance to all of these neighbours
      distance = sqrt((MPM (trackCount , trackEnd (trackCount)-1) - ...
                       MPM (neighbours(:) + trackCount , trackEnd(trackCount) + 1) ).^2 + ...
                      (MPM (trackCount , trackEnd (trackCount)) - ...
                       MPM (neighbours(:) + trackCount , trackEnd(trackCount) + 2) ).^2);
         
      % Take the smallest of these distances                        
      [minDist, minDistIndex] = min (distance);                        
      
      % If this distance is smaller than the specified max dist on the GUI
      if minDist < maxRelinkDist
                
         % Add the found track to the current track
         MPM (trackCount, (trackEnd (trackCount) + 1):end) = MPM (neighbours (minDistIndex) + trackCount, ...
                                                                 (trackEnd (trackCount) + 1):end);
                
         % Erase the redundant track
         MPM(neighbours(minDistIndex)+trackCount , :) = 0;
                
         % Make sure the current track gets processed again, with the new stop location. 
         % Maybe we can link it to yet another track 
         trackEnd(trackCount) = trackEnd(neighbours(minDistIndex)+trackCount);
                
         % Erase begin and stop location of our allocated track (now
         % part of the current track and no longer an individual track)
         trackEnd(neighbours(minDistIndex)+trackCount) = 0;
         trackStart(neighbours(minDistIndex)+trackCount) = 0;
                
         % Since we erased a track, decrease the track counter
         trackCount = trackCount - 1;
              
         % We linked two tracks together
         linked = 1; 
      end 
   end
    
   % If we still couldn't link any of the neighbours
   if ~linked
       
      % Look for tracks which begin two frames after our current
      % track stops and try to link 
      gaps = find( (trackEnd (trackCount) + 3.9) < size(MPM,1) & ...
                  ((trackEnd (trackCount) - trackStart (trackCount + 1:end)) == -2) );
                    
      % If we do find any of these do:
      if ~isempty (gaps)
          
         % Calculate the distance to these tracks
         distance = sqrt((MPM (trackCount , trackEnd (trackCount)-1) - ...
                          MPM (gaps(:) + trackCount , trackEnd(trackCount) + 3) ).^2 + ...
                         (MPM (trackCount , trackEnd (trackCount)) - ...
                          MPM (gaps(:) + trackCount , trackEnd(trackCount) + 4) ).^2);
                      
          % Take the smallest of these distances          
         [minDist, minDistIndex] = min (distance);
         
         % If this distance is smaller than the specified max dist on the GUI       
         if minDist < maxRelinkDist
                        
            % Add the found track to the current track
            MPM (trackCount, (trackEnd (trackCount) + 3):end) = MPM (gaps (minDistIndex) + trackCount, ...
                                                                    (trackEnd (trackCount) + 3):end);
                      
            % Fill in the gap
            MPM (trackCount, trackEnd (trackCount) + 1) = MPM (trackCount, trackEnd (trackCount) - 1) + ...
                                                          round ((MPM (trackCount, trackEnd (trackCount) + 3) - ...
                                                          MPM (trackCount, trackEnd (trackCount) - 1)) / 2);
            MPM (trackCount, trackEnd (trackCount) + 2) = MPM(trackCount , trackEnd(trackCount)) +...
                                                          round ((MPM (trackCount, trackEnd (trackCount) + 4) - ...
                                                          MPM (trackCount, trackEnd (trackCount))) / 2);
            % Erase the redundant track
            MPM (gaps (minDistIndex) + trackCount, :) = 0;
                        
            % Make sure the current track gets processed again, with the new stop location. 
            % Maybe we can link it to yet another track
            trackEnd (trackCount) = trackEnd (gaps (minDistIndex) + trackCount);
                        
            % Erase begin and end locations of our allocated track (now
            % part of the current track and no longer an individual track)
            trackEnd (gaps (minDistIndex) + trackCount) = 0;
            trackStart (gaps (minDistIndex) + trackCount) = 0;
                        
            % Since we erased a track, decrease the track counter
            trackCount = trackCount - 1;
            
         end   % if minDist < maxRelinkDist              
      end   % if ~isempty (gaps)     
   end   % if ~linked 
end   % while trackCount

% Just in case there are totally empty rows in MPM again, erase them
% First find all unique rows unequal to zero
[notZeroEntryRows, notZeroEntryCols] = find (MPM);
notZeroEntryRows = unique (notZeroEntryRows);

% Get these entries from MPM and throw away the rest
MPM = MPM (notZeroEntryRows,:);

% Next clear out tracks that are shorter than the min amount of frames
% specified on the GUI
[rows, cols] = find (MPM);

sortedTracks = sort (rows);
[uniqueEntries, uniqueIdx] = unique (sortedTracks);

% uniqueIdx returns the last occurence of the respective unique entry
% having sorted rows before, we can now count the number of occurences

if size (uniqueEntries,1) > size (uniqueEntries,2);
   uniqueIdx = [0;uniqueIdx];
else
   uniqueIdx = [0,uniqueIdx];
end 

% Get the number of tracks
numberOfOccurences = diff (uniqueIdx); 

% Since MPM has x and y coordinates for every entry we should test against
% (minTrackLength * 2 - 1)
minDistIndex = uniqueEntries (find (numberOfOccurences < minTrackLength*2-1));

% Clear these short tracks out of MPM
MPM(minDistIndex,:) = [];

% Store the modified MPM matrix in the data* directory
if ~isempty (savePath)
   cd (savePath);
   save ('MPM.mat', 'MPM');
end

% Tell the user that we've finished relinking
fprintf (1, 'Filtering and relinking finished!\n');
