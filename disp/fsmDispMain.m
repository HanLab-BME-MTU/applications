function [fsmParam,status,SCORE]=fsmDispMain(fsmParam,SCORE)
% fsmDispMain is the main function of the fsmDisp module
%
% SYNOPSIS   [fsmParam,status]=fsmDispMain(fsmParam,speckleArray)
%
% INPUT      fsmParam     :   general parameter structure
%            SCORE        :   scores rearranged into a matrix with the form
%                             [t y x s]n     t : time point (frame)
%                                            y : y coordinate of the event position
%                                            x : x coordinate of the event position
%                                            s : score
%                                            n : total number of events
%
% OUTPUT     fsmParam     :   modified (when needed) parameter structure
%            status       :   it decides whether the program should continue after this module 
%                             or should stop because of errors;
%                             status is set to 0 (error) in the beginning of a module and will
%                             be set to 1 at the end if the module completed successfully.
%                             status = 1 - if the module completed successfully
%                             status = 0 - if the module did not complete successfully
%
% DEPENDENCES   fsmKinMain uses {}
%               fsmKinMain is used by { fsmMain }
%
% Aaron Ponti, October 9th, 2002

% Set initial module status
status=0;

% Check input parameter
if nargin~=2
    error('Two input parameters expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CREATE IMAGES WITH OVERLAID COLOR-CODED SCORES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsmDispScoreMaps(SCORE,fsmParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING MODULE STATUS TO 1 AND RETURNING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the status to 1 to mean that the module successfully finished
status=1;
