function [movieNr, movieLength] = ptMinMovieLength (allValidFrames)
% ptMinMovieLength returns the length of the shortest movie in the list of
% movies identified by validFrames
%
% SYNOPSIS       [movieNr, movieLength] = ptMinMovieLength (allValidFrames) 
%
% INPUT          allValidFrames : cell containing a number of validFrames arrays
%                
% OUTPUT         movieNr : the number of the shortest movie in the list
%                movieLength : the length of the shortest movie
%
% DEPENDENCIES   ptMinMovieLength uses { nothing }
%                                  
%                ptMinMovieLength is used by {  }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Nov 04          Initial version

% Returns the length of the shortest MPM in the list of MPM's

% Initialize vars
prevLength = 10000;  % Should be larger than the biggest possible movie
movieNr = 0;
movieLength = 0;
      
% Go throught the list of MPMs
for iCount = 1 : length(allValidFrames)
   
   % Test for length and keep if shorter
   %curLength = size(allvalidFrames{iCount},2);
   curLength = allValidFrames{iCount}(1,end) - allValidFrames{iCount}(1,1) + 1;
   if curLength < prevLength
      prevLength = curLength;
  
      movieNr = iCount;
      movieLength = curLength;
   end
end
