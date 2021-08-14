function tracksNA = removeOutOfRangeAmpFromTrack(tracksNA)
%function tracksNA = removeOutOfRangeAmpFromTrack(tracksNA) removes amp
%that is outside startingFrameExtraExtra and endingFrameExtraExtra and save
%to tracksNA again.
%   Sangyoon Han, August 2021
tic
nTracks = numel(tracksNA);
for ii=1:nTracks
    sF=tracksNA(ii).startingFrameExtraExtra;
    eF=tracksNA(ii).endingFrameExtraExtra;
    nFrames = tracksNA(ii).iFrame(end);
    % Making out of range NaNs
    if sF>1
        tracksNA(ii).amp(1:sF-1)=NaN;
    end
    if eF<nFrames
        tracksNA(ii).amp(eF+1:nFrames)=NaN;
    end
end
toc
end

