
function [anglesFinal meanAngleVector,stdAngleVector] =plusTipPlotTrackAngles2(xCoord,yCoord)
% plusTipPlotTrackAngles read the project tracks, get the angles and plot them
%
% Input:
% 
%   xCoord/yCoord: coordinates in structure such that xCoord/YCoord
%   (trackID (row) frame Number(col)) 
% anglesFinal: vector of orientations in degrees using x/y coordinate system (4 quadrant)
% of all frame to frame displacements (in future probably should have some 
% displacement criteria? 
% 
% meanAngleVector: average orientation of each subtrack (length of # of
% subtracks) (averaged over frame to frame displacements)
% 
% stdAngleVector: std of frame-2-frame displacements of each subtrack (length of 
% # of subtracks)
%
% Sebastien Besson, August 2011
% Maria Bagonis, October 2011
%
% Copyright (C) 2011 LCCB 
%
% This file is part of plusTipTracker.
% 
% plusTipTracker is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% plusTipTracker is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with plusTipTracker.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Input check
%ip =inputParser;
%ip.addRequired('xCoord',@isdouble)
%ip.addRequired('yCoord',@isdouble)
%ip.addRequired('saveDir',@ischar);
%ip.parse(xCoord,yCoord);

% Load tracking result
%trackDir = [projData.anDir filesep 'track'];
%try
 %   files = dir([trackDir filesep '*.mat']);
  %  s=load([trackDir filesep files(1).name],'tracksFinal');
%catch ME
 %   error(['--plusTipPostTracking:' ME.message]);
%end
% ask sebastian about better ways to write this
[nTracks nFrames] = size(xCoord(:,1));
meanAngleVector = zeros(nTracks,1); 
stdAngleVector = zeros(nTracks,1); 

for iTrack = 1:nTracks
x = xCoord(iTrack,:);
y = yCoord(iTrack,:); 
deltaX = x(2:end)-x(1:end-1); 
deltaY = -(y(2:end)-y(1:end-1)); 

avgDeltaX = nanmean(deltaX); 
avgDeltaY = nanmean(deltaY); 

meanAngleVector(iTrack)= rad2deg(atan2(avgDeltaY, avgDeltaX)); 

% use negative y here to correct for the y coordinate system of images
% so it is easier for the user to interpret 
angles = atan2(-(y(2:end)-y(1:end-1)),x(2:end)-x(1:end-1)); 
% convert to degrees 
anglesDeg = (180/pi)*(angles(:)) ;
% convert to four quadrant degrees
for i = 1:length(anglesDeg)
    if anglesDeg(i) <0 
        anglesDeg(i) = 360+anglesDeg(i);
    end 

end


% save the mean orientation for each subtrack in a vector 
%meanAngleVector(iTrack,1) =nanmean(anglesDeg(:)) ;
%stdAngleVector(iTrack,1) = nanstd(anglesDeg(:)); 

% angleForPlotting Rose are in Radians (right now the Rose plot 
% plots all frame to frame orientations NOT the average subtrack orientations) 

if iTrack == 1 
    anglesFinal = angles ;
else 
    anglesFinal = [anglesFinal; angles];
end 
end 

% Get all track angles and concatenate them
%%angles = arrayfun(@getAngles,s.tracksFinal,'UniformOutput',false);
%angles =horzcat(angles{:});
%angles=angles(~isnan(angles));
anglesFinal= anglesFinal(~isnan(anglesFinal)); 


end

%function angles = getAngles(track)
%x=track.tracksCoordAmpCG(1:8:end);
%y=track.tracksCoordAmpCG(2:8:end);
%angles=atan2(y(2:end)-y(1:end-1),x(2:end)-x(1:end-1));
%end



