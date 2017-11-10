function [mpmNr, mpmLength] = ptMinMPMLength (allMPM)
% ptMaxMPMLength returns the length of the shortest MPM in the list of MPM's
%
% SYNOPSIS       [mpmNr, mpmLength] = ptMinMPMLength (allMPM) 
%
% INPUT          allMPM : cell containing a number of MPM matrices
%                
% OUTPUT         mpmNr : the number of the shortest MPM in the list
%                mpmLength : the length of the shortest MPM
%
% DEPENDENCIES   ptMinMPMLength  uses { nothing }
%                                  
%                ptMinMPMLength is used by {  }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial version

% Returns the length of the shortest MPM in the list of MPM's

% Initialize vars
prevLength = 10000;  % Should be larger than the biggest possible MPM
mpmNr = 0;
mpmLength = 0;
      
% Go throught the list of MPMs
for iCount = 1 : length(allMPM)
   
   % Test for length and keep if shorter
   curLength = size(allMPM{iCount},2);
   if curLength < prevLength
      prevLength = curLength;
  
      mpmNr = iCount;
      mpmLength = curLength;
   end
end
