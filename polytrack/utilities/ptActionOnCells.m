function ptActionOnCells
% ptActionOnCells is the callback function for when the user clicks on a
% cell in the slider window
%
% SYNOPSIS       ptActionOnCells
%
% INPUT          none (it gets values from the figure object created in ptManualPostProcessJob)
%
% OUTPUT         none (it updates the handles object (including MPM) directly)
%
% DEPENDENCIES   ptActionOnCells uses {nothing}
%                                  
%                ptActionOnCells is used by { ptManualPostProcessJob,
%                                             ptShowSlidingFrames }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Feb 05          Initial release

fprintf(1,'Pressed cell\n');