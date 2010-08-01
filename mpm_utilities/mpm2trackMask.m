function trackMask = mpm2trackMask(MPM,keepCompleteTracksOnly)

% Get the track mask
[nrows ncols] = size(MPM);
nFrames = ncols / 2;
trackID = (1:nrows)';
trackMask = MPM(:,1:2:end) ~= 0;
    
if keepCompleteTracksOnly
    % remove any track that begins at 1st frame
    startAtFirstFrame = trackMask(:,1);
    trackMask(:,1) = false;
    for iFrame = 2:nFrames-1
        idxDead = trackID(~trackMask(:,iFrame));
        trackMask(:,iFrame) = trackMask(:,iFrame) & ~startAtFirstFrame;
        startAtFirstFrame(idxDead) = false;
    end
    % remove any track that ends at last frame
    endAtLastFrame = trackMask(:,end);
    trackMask(:,end) = false;
    for iFrame = nFrames-1:-1:1
        idxDead = trackID(~trackMask(:,iFrame));
        trackMask(:,iFrame) = trackMask(:,iFrame) & ~endAtLastFrame;
        endAtLastFrame(idxDead) = false;
    end
end
