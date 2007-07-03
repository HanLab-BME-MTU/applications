function cord = makiCoord2Cord(initCoord)
%MAKICOORD2CORD converts initCoord to a cord-array for the detector
%
% SYNOPSIS: cord = makiCoord2Cord(initCoord)
%
% INPUT initCoord: initCoord-cell from makiInitCoord
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

cord(1:nTimepoints) = ...
    struct('sp',[],'COM',[]);
for t=1:nTimepoints
    cord(t).sp.cord = job(iJob).dataStruct.initCoord{t}(:,1:3);
end