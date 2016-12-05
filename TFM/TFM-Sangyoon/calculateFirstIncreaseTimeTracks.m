function tracksNA = calculateFirstIncreaseTimeTracks(tracksNA,splineParamInit,preDetecFactor,tInterval)
%  tracksNA = calculateFirstIncreaseTimeTracks(tracksNA,splineParamInit) 
    if nargin<4
        tInterval=1; %for vinculin 
    end
    if nargin<3
        preDetecFactor=1/5; %for vinculin 
    end
    if nargin<2
        splineParamInit=1;
    end
    useSmoothing=false;
    if splineParamInit<1
        useSmoothing=true;
    end
%     for ii=curIndices
    for ii=1:numel(tracksNA)
        curTrack = tracksNA(ii);
        sFEE = curTrack.startingFrameExtraExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
        % See how many frames you have before the startingFrameExtra
        numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
        numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
        numPreSigStart = min(20,numFramesBefore);
        sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
        sF10before = max(sFEE,curTrack.startingFrameExtra-5*numPreFrames);
        d = tracksNA(ii).ampTotal;
        nTime = length(d);
        if useSmoothing
%             d = tracksNA(ii).ampTotal;
            tRange = tracksNA(ii).iFrame;
            d(d==0)=NaN;
            warning('off','SPLINES:CHCKXYWP:NaNs')
            try
                sd_spline= csaps(tRange,d,splineParamInit);
            catch
                d = tracksNA(ii).amp;
                d(sFEE:curTrack.endingFrameExtraExtra) = tracksNA(ii).ampTotal(sFEE:curTrack.endingFrameExtraExtra);
                sd_spline= csaps(tRange,d,splineParamInit);
            end
            sd=ppval(sd_spline,tRange);
            sd(isnan(d))=NaN;

%             tRange = 1:nTime;
%             numNan = find(isnan(d),1,'last');
%             if isempty(numNan)
%                 numNan=0;
%             end
%             tRange(isnan(d)) = [];
%             d(isnan(d)) = [];
%             sd_spline= csaps(tRange,d,splineParamInit);
%             sd=ppval(sd_spline,tRange);
%             d = [NaN(1,numNan) d];
%     %         tRange = [NaN(1,numNan) tRange];
%             sd = [NaN(1,numNan) sd];
            bkgMaxInt = nanmax(sd(sF10before:sF5before));
            firstIncreseTimeInt = find(sd>bkgMaxInt & 1:length(sd)>sF5before,1);
        else
            bkgMaxInt = max(curTrack.ampTotal(sF10before:sF5before));
            firstIncreseTimeInt = find(curTrack.ampTotal>bkgMaxInt & 1:nTime>sF5before,1);
        end
%         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstIncreseTimeInt)
            if useSmoothing
                curForce=d;
                curForce(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra) = tracksNA(ii).forceMag(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra);
%                 curForce(isnan(curForce)) = [];
                sCurForce_spline= csaps(tRange,curForce,splineParamInit);
                sCurForce_sd=ppval(sCurForce_spline,tRange);
                sCurForce_sd(isnan(curForce))=NaN;
%                 sCurForce = [NaN(1,numNan) sCurForce];
                bkgMaxForce = nanmax(sCurForce_sd(sF10before:sF5before));
                firstIncreseTimeForce = find(sCurForce_sd>bkgMaxForce & 1:length(sCurForce_sd)>sF5before,1);
            else
                bkgMaxForce = max(curTrack.forceMag(sF10before:sF5before));
                firstIncreseTimeForce = find(curTrack.forceMag>bkgMaxForce & 1:nTime>sF5before,1);
            end
            if isempty(firstIncreseTimeForce) || firstIncreseTimeForce>curTrack.endingFrameExtraExtra
                tracksNA(ii).forceTransmitting = false;
                tracksNA(ii).firstIncreseTimeInt = [];
                tracksNA(ii).firstIncreseTimeForce = [];
                tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            else
                tracksNA(ii).forceTransmitting = true;
                tracksNA(ii).firstIncreseTimeInt = firstIncreseTimeInt*tInterval; % in sec
                tracksNA(ii).firstIncreseTimeForce = firstIncreseTimeForce*tInterval;
                tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeInt*tInterval - firstIncreseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            end
        else
            tracksNA(ii).forceTransmitting = false;
            tracksNA(ii).firstIncreseTimeInt = [];
            tracksNA(ii).firstIncreseTimeForce = [];
            tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
            tracksNA(ii).bkgMaxInt = bkgMaxInt;
            tracksNA(ii).bkgMaxForce = [];
        end
    end
