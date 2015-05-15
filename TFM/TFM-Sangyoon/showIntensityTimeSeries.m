function []=showIntensityTimeSeries(tracksNA,selTrackID)
    h=figure;
    curEndFrame=tracksNA(selTrackID).endingFrame;
    curStartFrame = tracksNA(selTrackID).startingFrame;
%     curTimeRange = 1:(curEndFrame-curStartFrame+1);
    curFrameRange= curStartFrame:curEndFrame;
    subplot(1,2,1), plot(curFrameRange,tracksNA(selTrackID).ampTotal(curFrameRange))
    subplot(1,2,2), plot(curFrameRange,tracksNA(selTrackID).forceMag(curFrameRange),'r')
    set(h,'Position',[900,300,400,150]),title(['ID:' num2str(selTrackID) ', CC-score:' num2str(tracksNA(selTrackID).CCscore)])
end

