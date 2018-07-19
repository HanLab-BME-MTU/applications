function [MPM, lengthMPM] = ptAddMPM (allMPM)
% ptAddMPM adds the MPM matrix newMPM to currentMPM and fills up with zeros if necessary 
%
% SYNOPSIS       [MPM, lengthMPM] = ptAddMPM (currentMPM, newMPM)
%
% INPUT          allMPM : all the MPM's in a cell
%                
% OUTPUT         MPM : the merged input MPM's
%                lengthMPM : the length of the longest MPM (column wise)
%
% DEPENDENCIES   ptAddMPM   uses {nothing}
%                                  
%                ptAddMPM  is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial release

for iCount = 1 : length (allMPM)
   if (iCount+1) <= length (allMPM)
       
      % Assign MPM's
      currentMPM = allMPM{iCount};
      newMPM = allMPM{iCount+1};
       
      % See which one is longest
      sizeCurrent = size(allMPM{iCount},2);
      sizeNew = size(allMPM{iCount+1},2);
      if sizeCurrent > sizeNew
         % Create a zeros matrix that we can cat to the shortest one so
         % that they become equal in size
         zeroMPM = zeros(size(newMPM,1), (sizeCurrent - sizeNew));
         
         catMPM = [newMPM zeroMPM];
         newMPM = catMPM;
         
         % Keep the size of the longest MPM
         lengthMPM = sizeCurrent;
      elseif sizeNew > sizeCurrent
         % Create a zeros matrix that we can cat to the shortest one so
         % that they become equal in size
         zeroMPM = zeros(size(currentMPM,1), (sizeNew - sizeCurrent));
         
         catMPM = [currentMPM zeroMPM];
         currentMPM = catMPM;
         
         % Keep the size of the longest MPM
         lengthMPM = sizeNew;
      else
         % It doesn't matter which length we take because they are equal
         lengthMPM = sizeNew;
      end
      
      % Cat the new and current MPMs together
      MPM = [currentMPM ; newMPM];
   end
end
