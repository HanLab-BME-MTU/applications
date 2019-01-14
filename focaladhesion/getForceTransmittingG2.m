function idGroup2f = getForceTransmittingG2(idGroup2,tracksNAG2,tInterval)
% input:
%           idGroup2        logical indices for G2
%           tracksNAG2      specific for only G2, equivalent as yracksNA(idGroup2)
% Sangyoon Han

ampSlopeG2 = zeros(sum(idGroup2),1); %arrayfun(@(x) x.ampSlope,tracksNA(idGroup2));
ii=0;
for pp=1:numel(tracksNAG2) %find(idGroup2)'
    ii=ii+1;
    tracksNAG2(pp).lifeTime = tracksNAG2(pp).endingFrameExtra - tracksNAG2(pp).startingFrameExtra;
    curTrack = tracksNAG2(pp);
    
    sF=curTrack.startingFrameExtraExtra; eF=curTrack.endingFrameExtraExtra;
    [~,ampSlopeG2(ii)] = regression((1:eF-sF+1),curTrack.ampTotal(sF:eF));
end
% initForceG2 = arrayfun(@(x) x.forceMag(x.startingFrame),tracksNAG2);
lifeTimeG2 = arrayfun(@(x) x.lifeTime,tracksNAG2);
ampEndingG2 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNAG2);
ampStartingG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),tracksNAG2);
idxIncreasingAmpG2 = ampSlopeG2>0 & ampEndingG2>ampStartingG2;
% idxLowInitForceG2= initForceG2<500;
% maturing NA should have at least 3 min of lifetime
thresLT_G2_frames = 2*60/tInterval;
idxLongLifeTimeG2=lifeTimeG2>thresLT_G2_frames;
idGroup2valid = idxIncreasingAmpG2 & idxLongLifeTimeG2; %& idxLowInitForceG2 
idGroup2f=false(size(idGroup2)); idGroup2index=find(idGroup2);
idGroup2f(idGroup2index(idGroup2valid))=true;

end
