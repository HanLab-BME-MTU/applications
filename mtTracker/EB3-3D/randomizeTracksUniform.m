function [randTracks]=randomizeTracksUniform(tracks,maxRandomDist)
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
   randDir=rand(3,1)-0.5;
   randPoint=maxRandomDist*randDir;
   track.x=track.x + randPoint(1);
   track.y=track.y + randPoint(2);
   track.z=track.z + randPoint(3);
end

