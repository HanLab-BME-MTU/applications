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
   randDir=rand(1,3)-0.5;
   randDir=fixedDistance*randDir/norm(randDir);
   track.x=track.x + randDir(1);
   track.y=track.y + randDir(2);
   track.z=track.z + randDir(3);
end

