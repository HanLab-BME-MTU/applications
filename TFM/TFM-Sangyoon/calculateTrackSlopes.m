function tracksNA = calculateTrackSlopes(tracksNA,tInterval)
% tracksG1 = calculateTrackSlopes(tracksG1) calculates slopes of force and
% amplitude from first minute or half lifetime.
% Sangyoon Han, December, 2015

tIntervalMin = tInterval/60;
prePeriodFrame = ceil(10/tInterval); %pre-10 sec
periodFrames = 30;
for k=1:numel(tracksNA)
    sF=max(tracksNA(k).startingFrameExtra-prePeriodFrame,tracksNA(k).startingFrameExtraExtra);
    
%     eF=tracksNA(k).endingFrameExtra;
    curLT = tracksNA(k).lifeTime;
    halfLT = ceil(curLT/2);
    earlyPeriod = min(halfLT,floor(20/tInterval)); % frames per 20 sec or half life time
    oneMinPeriod = min(halfLT,floor(60/tInterval)); % frames per a minute or half life time

    lastFrame = min(tracksNA(k).endingFrameExtraExtra,sF+earlyPeriod+prePeriodFrame-1);
    lastFrameFromOne = lastFrame - sF+1;

    lastFrameOneMin = min(tracksNA(k).endingFrameExtraExtra,sF+oneMinPeriod+prePeriodFrame-1);
    lastFrameFromOneOneMin = lastFrameOneMin - sF+1;
    
    try
        [~,curM] = regression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).amp(sF:lastFrame));
        tracksNA(k).earlyAmpSlope = curM; % in a.u./min
    catch
        tracksNA(k).earlyAmpSlope = NaN;
    end
    if isfield(tracksNA,'ampTotal2')
        [~,curM2] = regression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).ampTotal2(sF:lastFrame));
        tracksNA(k).earlyAmpSlope2 = curM2; % in a.u./min
        [~,curM2] = regression(tIntervalMin*(1:lastFrameFromOneOneMin),tracksNA(k).ampTotal2(sF:lastFrameOneMin));
        tracksNA(k).ampSlope2 = curM2; % in a.u./min
    end
    if isfield(tracksNA,'forceMag')  
        [~,curForceM] = regression(tIntervalMin*(1:lastFrameFromOneOneMin),tracksNA(k).forceMag(sF:lastFrameOneMin));
    %         figure, plot(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame))
    %         figure, plotregression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame))
        tracksNA(k).forceSlope = curForceM; % in Pa/min

        curEndFrame = min(sF+periodFrames-1,tracksNA(k).endingFrame);
        curEarlyPeriod = curEndFrame - sF+1;
        [~,curForceEarlySlopeGroup] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).forceMag(sF:curEndFrame));

        tracksNA(k).earlyForceSlope = curForceEarlySlopeGroup; % in Pa/min
    end
end
