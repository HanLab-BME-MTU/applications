function [tracksNA, lateForceSlopes] = calculateTrackLateForceSlope(tracksNA,tInterval)
% tracksNA = calculateTrackLateForceSlope(tracksNA,tInterval) calculates
% late force slope during the last minute or half lifetime.
% Sangyoon Han, December, 2015
tIntervalMin = tInterval/60;
prePeriodFrame = ceil(10/tInterval); %pre-10 sec
lateForceSlopes = zeros(numel(tracksNA),1);
for k=1:numel(tracksNA)
    eF=min(tracksNA(k).endingFrameExtra+prePeriodFrame,tracksNA(k).endingFrameExtraExtra);
    curLT = tracksNA(k).lifeTime;
    halfLT = ceil(curLT/2);
    latePeriod = min(halfLT,floor(60/tInterval)); % frames per a minute or half life time

    earlyFrame = max(tracksNA(k).startingFrameExtraExtra,eF-latePeriod-prePeriodFrame);
    earlyFrameFromLast = eF-earlyFrame;
    
    [~,curFS] = regression(tIntervalMin*(1:earlyFrameFromLast),tracksNA(k).forceMag(earlyFrame:eF));
%         figure, plot(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame))
%         figure, plotregression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame))
    tracksNA(k).lateForceSlope = curFS; % in Pa/min
    lateForceSlopes(k)=curFS;
end
