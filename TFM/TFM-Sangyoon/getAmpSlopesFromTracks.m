function [ ampSlopes ] = getAmpSlopesFromTracks( tracksNA )
%[ ampSlopes ] = getAmpSlopesFromTracks( tracksNA ) calculates slopes in
%amplitude in tracksNA via regression.
%   input:
%               tracksNA: tracksNA
%   output:
%               ampSlopes: Nx1 array containing slopes.
ampSlopes = zeros(numel(tracksNA),1);
for ii=1:numel(tracksNA)
    curTrack = tracksNA(ii);
    [~,ampSlopes(ii)] = regression((1:curTrack.lifeTime+1),curTrack.amp(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
end
end

