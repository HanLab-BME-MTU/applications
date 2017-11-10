function tracksNA = recalculateLifeTimeTracks(tracksNA)
numTracks=numel(tracksNA);
for k=1:numTracks
    try
        sF=tracksNA(k).startingFrameExtra;
        eF=tracksNA(k).endingFrameExtra;
        if isempty(sF)
            sF=tracksNA(k).startingFrame;
        end
        if isempty(eF)
            eF=tracksNA(k).endingFrame;
        end
    catch
        sF=tracksNA(k).startingFrame;
        eF=tracksNA(k).endingFrame;
        tracksNA(k).lifeTime = eF-sF;    
    end
    tracksNA(k).lifeTime = eF-sF;    
end