function MPM = ptMinTrackLength (MPM, minTrackLength)
% ptMinTrackLength removes all tracks shorter than minTrackLength from MPM
%
% SYNOPSIS       MPM = ptMinTrackLength (MPM, minTrackLength)
%
% INPUT          MPM        : Magic Position Matrix 
%                               MPM = [ y  x  y  x  y  x ... ]
%                                        t1    t2    t3
% 				 minTrackLength : minimal length of tracks
%
% OUTPUT         MPM : a modified MPM matrix
%
% DEPENDENCIES   ptMinTrackLength  uses {nothing}
%                                  
%                ptMinTrackLength is used by { ptTrackCells }
%

% Write a message to the screen
fprintf (1, '\n     ptMinTrackLength: making sure that tracks have minimum length of %d...\n', minTrackLength);

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
MPM (minDistIndex,:) = [];

% Tell the user that we've finished relinking
fprintf (1, '     ptMinTrackLength: finished!\n');
