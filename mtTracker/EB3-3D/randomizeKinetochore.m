function [randTracks]=randomizeKinetochore(tracks,fixedDistance)
% This function randomize the position of the kinTracks. 
% INPUT
% - kinTracks: trajectory in isotropized pixel space
% - fixed distance (pixel): number of pixel used to change the track coordinate in a random direction. 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('tracks',@(x) isa(x,'Tracks'));
ip.parse(tracks);
p=ip.Results;

randTracks=tracks.copy();
for kIdx=1:length(tracks)
   track=randTracks(kIdx);
   track.x=track.x + sign((rand(1)-0.5))*fixedDistance;
   track.y=track.y + sign((rand(1)-0.5))*fixedDistance;
   track.z=track.z + sign((rand(1)-0.5))*fixedDistance;
end

