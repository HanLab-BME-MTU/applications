function fsmParam=fsmMarkModuleUpToDate(fsmParam)
% fsmMarkModuleUpToDate marks modules in fsmParam as up-to-date (fsmParam.module.uptodate=1) or not (fsmParam.module.uptodate=0)
%
% SYNOPSIS      fsmParam=fsmMarkModuleUpToDate(fsmParam)
%
% INPUT         fsmParam    : parameter structure created and used by fsmGuiMain (speckTackle)
%
% OUTPUT        fsmParam    : updated fsmParam structure. The fields fsmParam.{module}.uptodate, where
%                             {module} is one of {'prep','track',build','kin','disp'} are set to 1 if 
%                             the corresponding module is up-to-date or 0 otherwise.
%
% DEPENDENCES   fsmMarkModuleUpToDate uses { }
%               fsmMarkModuleUpToDate is used by { }
%
% Aaron Ponti, August 27th, 2004


if nargin~=1
    error('One input parameter - fsmParam - expected.');
end

% Extract information on which modules have been run
status=[fsmParam.prep.enable  ...
        fsmParam.track.enable ...
        fsmParam.build.enable ...
        fsmParam.kin.enable   ...
        fsmParam.disp.enable];

% Set flags according to which module is up-to-date
indx=find(status);
status(1:min(indx))=1; % All modules BEFORE the first run are up-to-date

% Update fsmParam
fsmParam.prep.uptodate  = status(1);
fsmParam.track.uptodate = status(2);
fsmParam.build.uptodate = status(3);
fsmParam.kin.uptodate   = status(4);
fsmParam.disp.uptodate  = status(5);



        
    