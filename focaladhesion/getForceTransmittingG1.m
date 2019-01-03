function idGroup1f = getForceTransmittingG1(idGroup1,tracksNAG1)
% input:
%           idGroup1        logical indices for G1
%           tracksNAG2      specific for only G2, equivalent as yracksNA(idGroup2)
% Sangyoon Han

% 1. Based on amplitude 
curAmpSlopeGroup = arrayfun(@(x) x.ampSlope, tracksNAG1);
curEarlyAmpSlopeGroup = arrayfun(@(x) x.earlyAmpSlope, tracksNAG1);
% We decided to regard amplitude with flat slope as noise.
%             figure, plot(curAmpSlopeGroup,curEarlyAmpSlopeGroup,'*')
indFlatAmp = curAmpSlopeGroup<=0 & curEarlyAmpSlopeGroup<=0;
% 2. Based on forceMag
curVeryEarlyAmpSlopeGroup = NaN(sum(idGroup1),1);
curForceSlopeGroup = NaN(sum(idGroup1),1);
curForceEarlySlopeGroup = NaN(sum(idGroup1),1);
% curAmpOverallSlopeGroup = NaN(sum(idGroup1),1);
% curAmpOverallSlopeUpErrGroup = NaN(sum(idGroup1),1);
% curAmpOverallSlopeDownErrGroup = NaN(sum(idGroup1),1);
curForceLateSlopeGroup = NaN(sum(idGroup1),1);

periodFrames = 30;
ii=0;
for pp = 1:numel(tracksNAG1) %find(idGroup1)'
    ii=ii+1;
    curTrack = tracksNAG1(pp);
    [~,curForceSlopeGroup(ii)] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
    curEndFrame = min(curTrack.startingFrameExtra+periodFrames-1,curTrack.endingFrame);
    curEarlyPeriod = curEndFrame - curTrack.startingFrameExtra+1;
    [~,curForceEarlySlopeGroup(ii)] = regression((1:curEarlyPeriod),curTrack.forceMag(curTrack.startingFrameExtra:curEndFrame));
    sF=curTrack.startingFrameExtraExtra; eF=curTrack.endingFrameExtraExtra;
    curEndFrame = min(sF+periodFrames-1,curTrack.endingFrame);
    curEarlyPeriod = curEndFrame - sF+1;
    [~,curVeryEarlyAmpSlopeGroup(ii)] = regression((1:curEarlyPeriod),curTrack.ampTotal(sF:curEndFrame));
    
%     tVar=(1:eFshifted)';
%     [curCurAmpOverallSlopeGroup,curCurAmpOverallSlopeErrGroup] = regress(curTrack.ampTotal(sF:eF)',[ones(size(tVar)) tVar]);
%     curAmpOverallSlopeGroup(ii)= curCurAmpOverallSlopeGroup(2);
%     curAmpOverallSlopeUpErrGroup(ii) = curCurAmpOverallSlopeErrGroup(2,2);
%     curAmpOverallSlopeDownErrGroup(ii) = curCurAmpOverallSlopeErrGroup(2,1);
    
    % There is amplitude that has rising phase but just stops in the middle
    % of course. (It means the the particle loses it's point-ness while
    % increasing the intensity. This does not help in characterizing true
    % G1 behavior
    % I'll filter these out with lateAmpSlope
    curStartFrame = max(eF-periodFrames,curTrack.startingFrame);
    curLatePeriod = curTrack.endingFrameExtraExtra - curStartFrame + 1;
    [~,curForceLateSlopeGroup(ii)] = regression((1:curLatePeriod),curTrack.ampTotal(curStartFrame:eF));
end
indFlatForce = curForceSlopeGroup<=0 & curForceEarlySlopeGroup<=0;
validTrackID = (~indFlatForce & ~indFlatAmp) & curVeryEarlyAmpSlopeGroup>0 ...
    & curForceLateSlopeGroup<0;
idGroup1f=false(size(idGroup1)); idGroup1index=find(idGroup1);
idGroup1f(idGroup1index(validTrackID))=true;

end

