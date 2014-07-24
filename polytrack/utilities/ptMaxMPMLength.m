function [mpmNr, mpmLength] = ptMaxMPMLength (allMPM)
% ptMaxMPMLength returns the length of the longest MPM in the list of MPM's
%
% SYNOPSIS       [mpmNr, mpmLength] = ptMaxMPMLength (allMPM) 
%
% INPUT          allMPM : cell containing a number of MPM matrices
%                
% OUTPUT         mpmNr : the number of the longest MPM in the list
%                mpmLength : the length of the longest MPM
%
% DEPENDENCIES   ptMaxMPMLength  uses { nothing }
%                                  
%                ptMaxMPMLength is used by {  }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

% Returns the length of the longest MPM in the list of MPM's

% Initialize vars
prevLength = 0;
mpmNr = 0;
mpmLength = 0;
      
% Go throught the list of MPMs
for iCount = 1 : length(allMPM)
   
   % Test for length and keep if longer
   curLength = size(allMPM{iCount},2);
   if curLength > prevLength
      prevLength = curLength;
  
      mpmNr = iCount;
      mpmLength = curLength;
   end
end
