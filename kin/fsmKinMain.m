function [fsmParam,status,speckleArray,SCORE]=fsmKinMain(fsmParam,speckleArray)
% fsmKinMain is the main function of the fsmKin module
%
% SYNOPSIS   [fsmParam,status]=fsmKinMain(fsmParam,speckleArray)
%
% INPUT      fsmParam     :   general parameter structure
%            speckleArray :   structure containing all speckle information from a movie 
%
% OUTPUT     fsmParam     :   modified (when needed) parameter structure
%            status       :   it decides whether the program should continue after this module 
%                             or should stop because of errors;
%                             status is set to 0 (error) in the beginning of a module and will
%                             be set to 1 at the end if the module completed successfully.
%                             status = 1 - if the module completed successfully
%                             status = 0 - if the module did not complete successfully
%            speckleArray :   speckleArray augmented with score/activity information
%            SCORE        :   scores rearranged into a matrix with the form
%                             [t y x s]n     t : time point (frame)
%                                            y : y coordinate of the event position
%                                            x : x coordinate of the event position
%                                            s : score
%                                            n : total number of events
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
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

userPath=fsmParam.main.path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE TO WORK PATH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(userPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CLASSIFY BIRTH AND DEATH EVENS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

speckleArray=fsmKinEventClassifier(speckleArray,fsmParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SAVE speckleArray STRUCTURE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoH=fsmGuiInfo; ch=get(infoH,'Children');
if length(userPath)>30
    strUserPath=[userPath(1:27),'...'];
else
    strUserPath=userPath;
end
textString=['Saving augmented speckleArray to ',strUserPath,'\speckleArray.mat'];
set(ch(3),'String','KINETIC ANALYSIS MODULE');
set(ch(2),'String',textString);
set(ch(1),'String','Prease wait...');
save speckleArray.mat speckleArray;
set(ch(1),'String','Done!');
pause(1);
close(infoH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CALCULATE SCORES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SCORE=fsmKinSaveScoreLists(speckleArray,fsmParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING MODULE STATUS TO 1 AND RETURNING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the status to 1 to mean that the module successfully finished
status=1;
