function [peakTimeIntAgainstForceAll,peakTimeIntAll,peakTimeForceAll,tracksNA] ...
    = calculatePeakTimeLagFromTracks(tracksNA,splineParam,tInterval,varargin)
% [peakTimeIntAgainstForceAll,peakTimeIntAll,peakTimeForceAll] ...
%     = calculatePeakTimeLagFromTracks(tracksNA,splineParam,tInterval,varargin)
%  calculates the time lag of the main intensity peak (ampTotal) against
%  the slave source's peak.
% input:
%       splineParamInit: smoothing parameter (0-1). Use 1 if you don't want to
%       smooth the signal.
%       'slaveSource': either 'forceMag', 'ampTotal2' or 'ampTotal3'
% output:
%       
% Big change: I gave up updating tracksNA becuase it increase too much 
% the file size
ip =inputParser;
ip.addRequired('tracksNA',@isstruct)
ip.addOptional('splineParam',0.1,@isscalar)
ip.addOptional('tInterval',1,@(x)isscalar(x))
ip.addParamValue('slaveSource','forceMag',@(x)ismember(x,{'forceMag','ampTotal2','ampTotal3'})); % collect NA tracks that ever close to cell edge
ip.parse(tracksNA,splineParam,tInterval,varargin{:});
slaveSource=ip.Results.slaveSource;

peakTimeIntAgainstForceAll=NaN(numel(tracksNA),1);
peakTimeIntAll=NaN(numel(tracksNA),1);
peakTimeForceAll=NaN(numel(tracksNA),1);
for ii=1:numel(tracksNA)
    sF = max(1,tracksNA(ii).startingFrameExtra-50);
    tRange = sF:tracksNA(ii).endingFrameExtraExtra;
    d = tracksNA(ii).ampTotal(tRange);
    warning('off','SPLINES:CHCKXYWP:NaNs')
    d(d==0)=NaN;
    try
        sd_spline= csaps(tRange,d,splineParam);
    catch
        d = tracksNA(ii).amp;
        d(tRange) = ...
            tracksNA(ii).ampTotal(tRange);
        sd_spline= csaps(tRange,d,splineParamInit);
    end
    sd=ppval(sd_spline,tRange);
    sd(isnan(d))=NaN;
    [~,curFrameMaxAmp]=nanmax(sd);
    if tracksNA(ii).lifeTime>31
        movingWindowSize=31;
    else
        movingWindowSize=5;
    end

    curFrameLocMaxes=locmax1d(sd,movingWindowSize);

    if ismember(curFrameMaxAmp,curFrameLocMaxes) && curFrameMaxAmp>tRange(1) && ...
            curFrameMaxAmp<tRange(end) && isfield(tracksNA(ii),slaveSource)
        peakTimeIntAll(ii) = curFrameMaxAmp;
        tracksNA(ii).intenPeakness = true;
        tracksNA(ii).intenPeakFrame = curFrameMaxAmp;

%         curForce=d;
        curForce = getfield(tracksNA(ii),{1},slaveSource,{tRange});
        if all(isnan(curForce))
            disp(['There is something wrong with ' slaveSource '. Please run the associated process again and run the step 11.'])
            peakTimeForceAll(ii) = NaN;
            tracksNA(ii).forcePeakMag = NaN; 
            tracksNA(ii).forcePeakFrame = NaN; 
            tracksNA(ii).forcePeakMagRel = NaN; % this is a relative difference
            peakTimeIntAgainstForceAll(ii) = NaN; %
            continue;
        end
        sCurForce_spline= csaps(tRange,curForce,splineParam);

        sCurForce_sd=ppval(sCurForce_spline,tRange);
        sCurForce_sd(isnan(curForce))=NaN;

        curForceLocMaxes=locmax1d(sCurForce_sd,movingWindowSize);
        % delete the first and last frame from loc maxes
        curForceLocMaxes = setdiff(curForceLocMaxes, [tracksNA(ii).startingFrameExtra tracksNA(ii).endingFrameExtra]);
        if ~isempty(curForceLocMaxes)
%                 [forceMacMax,indMaxForceAmongLMs] = max(sCurForce(curForceLocMaxes));
            forceMagLMCand=sCurForce_sd(curForceLocMaxes)';
            weightForceMag = (forceMagLMCand - min(sCurForce_sd))/(max(forceMagLMCand) - min(sCurForce_sd));
            [forceMacMax,indMaxForceAmongLMs] = min((abs(curForceLocMaxes-curFrameMaxAmp)).^0.5/length(tracksNA(ii).lifeTime)./weightForceMag);
%                 [forceMacMax,indMaxForceAmongLMs] = min(abs(curForceLocMaxes-curFrameMaxAmp)./weightForceMag);
            tracksNA(ii).forcePeakness = true;
            peakTimeForceAll(ii) = curForceLocMaxes(indMaxForceAmongLMs);
            tracksNA(ii).forcePeakMag = forceMacMax; 
            tracksNA(ii).forcePeakMagRel = forceMacMax - min(sCurForce_sd); % this is a relative difference
            tracksNA(ii).forcePeakFrame = peakTimeForceAll(ii); % this is a relative difference
            peakTimeIntAgainstForceAll(ii) = -peakTimeForceAll(ii)*tInterval + peakTimeIntAll(ii)*tInterval; % - means intensity comes first; + means force peak comes first
        else
            peakTimeForceAll(ii) = NaN;
            tracksNA(ii).forcePeakMag = NaN; 
            tracksNA(ii).forcePeakFrame = NaN; 
            tracksNA(ii).forcePeakMagRel = NaN; % this is a relative difference
            peakTimeIntAgainstForceAll(ii) = NaN; %
        end
    end
end

