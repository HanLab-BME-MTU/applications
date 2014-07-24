function cord = makiCoord2Cord(initCoord,forCutoff,timeRange)
%MAKICOORD2CORD converts initCoord to a cord-array for the detector
%
% SYNOPSIS: cord = makiCoord2Cord(initCoord)
%
% INPUT initCoord: initCoord-cell from makiInitCoord
%       forCutoff: only use the coordinates for establishing cutoff
%       timeRange: list of time points for which to fill cord (so that you
%           can avoid MMF for certain frames)
%
% OUTPUT cord: cord-structure used for detectSpots-subfunctions
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 02-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTimepoints = length(initCoord);
if nargin < 2 || isempty(forCutoff)
    forCutoff = false;
end
if ~forCutoff
    fname = 'allCoordPix';
else
    fname = 'data4MMF';
end

% write sp, COM. Should anything more be necessary, check spotfind.m for
% how the information is filled in
cord(1:nTimepoints) = ...
    struct('sp',[],'COM',[]);
for t=1:nTimepoints
    if ~isempty(initCoord(t).nSpots)
        for i=1:size(initCoord(t).(fname),1)
            cord(t).sp(i).cord = ...
                initCoord(t).(fname)(i,[2,1,3]);
        end
        cord(t).COM = ...
            mean(initCoord(t).(fname)(:,[2,1,3]),1);
    end
end