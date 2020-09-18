function [firstIncreseTimeIntAgainstSlaveAll,forceTransmittingAll...
    ,firstIncreaseTimeIntAll,firstIncreaseTimeSlaveAll,bkgMaxIntAll,bkgMaxSlaveAll, ...
    tracksNA, ampIncreasingAll] ...
    = calculateFirstIncreaseTimeTracks(tracksNA,numAveragingWind,preDetecPeriod,tInterval,varargin)
% function [firstIncreseTimeIntAgainstSlaveAll,forceTransmittingAll...
%     ,firstIncreaseTimeIntAll,firstIncreaseTimeSlaveAll,bkgMaxIntAll,bkgMaxSlaveAll, ...
%     tracksNA, ampIncreasingAll] ...
%     = calculateFirstIncreaseTimeTracks(tracksNA,numAveragingWind,preDetecPeriod,tInterval,varargin)
%  calculates the time lag of the main intensity (ampTotal) against the slave source.
% input:
%       splineParamInit: smoothing parameter (0-1). Use 1 if you don't want to
%       smooth the signal.
%       preDetecPeriod: how much ahead of the inital start point you want
%       to detect as a pre-signal maximum (default: 60 seconds). It is used as:
%       sF10before = effectiveSF - round(preDetecPeriod/tInterval)
%       Thus, if you put larger number, you are selecting the maximum from
%       much earlier time points.
%       'slaveSource': either 'forceMag', 'ampTotal2' or 'ampTotal3'
% output:
%       firstIncreseTimeIntAgainstSlaveAll,firstIncreseTimeIntAll,
%       firstIncreseTimeSlaveAll are all in the unit of time, not frame. To
%       get the frame, you should divide them with tInterval.
%       
% Big change: I gave up updating tracksNA becuase it increase too much 
% the file size
% Sangyoon Han, developed 2016, revised August 2018.
ip =inputParser;
ip.addRequired('tracksNA',@isstruct)
ip.addOptional('numAveragingWind',5,@isscalar)
ip.addOptional('preDetecPeriod',60,@(x)isscalar(x)) % in second
ip.addOptional('tInterval',1,@(x)isscalar(x))
ip.addOptional('plotEachTrack',false,@(x)islogical(x)||isempty(x))
ip.addParamValue('slaveSource','forceMag',@(x)ismember(x,{'forceMag','amp2','amp3','ampTotal2','ampTotal3'})); % collect NA tracks that ever close to cell edge
ip.parse(tracksNA,numAveragingWind,preDetecPeriod,tInterval,varargin{:}); %splineParamInit,preDetecParam
slaveSource=ip.Results.slaveSource;
useSmoothing=false;
if numAveragingWind>1
    useSmoothing=true;
end
firstIncreseTimeIntAgainstSlaveAll=NaN(numel(tracksNA),1);
forceTransmittingAll=false(numel(tracksNA),1);
firstIncreaseTimeIntAll=NaN(numel(tracksNA),1);
firstIncreaseTimeSlaveAll=NaN(numel(tracksNA),1);
bkgMaxIntAll=NaN(numel(tracksNA),1);
bkgMaxSlaveAll=NaN(numel(tracksNA),1);
differentInitialMargin=50;
ampIncreasingAll=false(numel(tracksNA),1);

for ii=1:numel(tracksNA)
    curTrack = tracksNA(ii);
%     curEarlyAmpSlope = curTrack.earlyAmpSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
    curTrack.lifeTime = curTrack.endingFrameExtra-curTrack.startingFrameExtra;
%     curAmpSlope = curTrack.ampSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
%     curForceSlope = curTrack.earlyAmpSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
    % See how many frames you have before the startingFrameExtra
    % (20181003) I decided to discard this not-clear method to determine
    % sF5before and sF10before because it a bit hard to explain and now 
%     stepFrame =5;
%     effectiveSF = curTrack.startingFrameExtra - stepFrame; sF5before=1; sF10before=1;
%     sFEE = max(1,curTrack.startingFrameExtra-differentInitialMargin); %curTrack.startingFrameExtraExtra;
%     while sF5before==sF10before
%         effectiveSF = effectiveSF + stepFrame;
%         numFramesBefore = effectiveSF - sFEE;
%         numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
%         numPreSigStart = min([20,numFramesBefore 3*numPreFrames]);
%         sF5before = max(effectiveSF-numPreSigStart,effectiveSF-numPreFrames);
%         sF10before = max(sFEE,effectiveSF-3*numPreFrames);
%     end
    effectiveSF = curTrack.startingFrameExtra; 
    sFEE = max(1,curTrack.startingFrameExtra-differentInitialMargin); %curTrack.startingFrameExtraExtra;
    sF5before = max(sFEE,effectiveSF -1); % So the variabl name should be sF1before
    sF10before = max(sFEE,effectiveSF - round(preDetecPeriod/tInterval)); % So the variabl name should be sF1before
    
%     while sF5before==sF10before
%         numFramesBefore = effectiveSF - sFEE;
%         numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
%         numPreSigStart = min([20,numFramesBefore 3*numPreFrames]);
%         sF5before = max(effectiveSF-numPreSigStart,effectiveSF-numPreFrames);
%         sF10before = max(sFEE,effectiveSF-3*numPreFrames);
%     end
    pp=0;
    lastPeriod = min(curTrack.lifeTime+1, 100);
    firstPeriod = 20; pStep=20;
    nSlopes = floor((lastPeriod-firstPeriod)/pStep);
    curEarlyAmpSlope = NaN(nSlopes,1);
    for  initialTime=firstPeriod:pStep:lastPeriod
        pp=pp+1;
        ealryFrames = min(initialTime, curTrack.endingFrameExtra-sF10before+1);
        [~,curEarlyAmpSlope(pp)] = regression((1:ealryFrames),curTrack.ampTotal(sF10before:sF10before+ealryFrames-1));
    end
%     [~,curForceSlope] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));

%     sFEE = curTrack.startingFrameExtraExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
    if any(curEarlyAmpSlope>0) %&& (numFramesBefore>1) % && curForceSlope>0 % intensity should increase. We are not 
        % interested in decreasing intensity which will 
        d = tracksNA(ii).amp;
        nTime = length(d);
        if useSmoothing
%             d = tracksNA(ii).ampTotal;
%             tRange = tracksNA(ii).iFrame;
            d(d==0)=NaN;
%             warning('off','SPLINES:CHCKXYWP:NaNs')
%             try
%                 sd_spline= csaps(tRange,d,splineParamInit);
%             catch
%                 d = tracksNA(ii).amp;
%                 d(sFEE:curTrack.endingFrameExtraExtra) = tracksNA(ii).ampTotal(sFEE:curTrack.endingFrameExtraExtra);
%                 sd_spline= csaps(tRange,d,splineParamInit);
%             end
%             sd=ppval(sd_spline,tRange);
%             sd(isnan(d))=NaN;
            % Trying the change point analysis with mean of 3 moving window
            windowSize = 5; 
            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
%             sd = filter(b,a,d);
            sd = conv(d,b,'same'); %filter(b,a,d);
%             zeroPhaseFilter = designfilt('lowpassiir','FilterOrder',12, ...
%                 'HalfPowerFrequency',0.15,'DesignMethod','butter');            
%             sd = filtfilt(zeroPhaseFilter,d); %filter(b,a,d);
            % End points correction
            sF = find(~isnan(d),1);
            eF = find(~isnan(d),1,'last');
            halfW = floor(windowSize/2);
            for jj=1:halfW
                % first point
                sd(sF+jj-1) = sd(sF+jj-1) * windowSize / (jj+halfW);
                sd(eF-jj+1) = sd(eF-jj+1) * windowSize / (jj+halfW);
            end
%             % median filtering - based outlier filtering
%             filteredSD = medfilt1(sd,13,'omitnan');
%             [outlierIdx] = detectOutliers(sd,9);
            
            tracksNA(ii).amp = sd;

            bkgMaxInt = nanmax(sd(sF10before:sF5before));
            firstIncreaseTimeInt = find(sd>bkgMaxInt & 1:length(sd)>sF5before,1);
        else
            bkgMaxInt = max(d(sF10before:sF5before));
            firstIncreaseTimeInt = find(d>bkgMaxInt & 1:nTime>sF5before,1);
        end
%         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstIncreaseTimeInt) && isfield(curTrack,slaveSource)
            ampIncreasingAll(ii)=true;
            if useSmoothing
                curSlave=d;
                curSlave(curTrack.startingFrameExtraExtra:curTrack.endingFrameExtraExtra) = ...
                 getfield(curTrack,{1},slaveSource,{curTrack.startingFrameExtraExtra:curTrack.endingFrameExtraExtra});
                %                 curForce(isnan(curForce)) = [];
                pp=0;
                curEarlyForceSlope = NaN(5,1);
                for  initialTime=10:10:50
                    pp=pp+1;
                    ealryFrames = min(initialTime, curTrack.endingFrameExtra-sF10before+1);
                    [~,curEarlyForceSlope(pp)] = regression((1:ealryFrames),curSlave(sF10before:sF10before+ealryFrames-1));
                end
                
                if any(curEarlyForceSlope>0) && sum(~isnan(curSlave))>1
%                     sCurForce_spline= csaps(tRange,curSlave,splineParamInit);
%                     sCurForce_sd=ppval(sCurForce_spline,tRange);
                    sCurForce_sd = conv(curSlave,b,'same'); %filter(b,a,d);
%                     sCurForce_sd = filter(b,a,curSlave);
                    for jj=1:halfW
                        % first point
                        sCurForce_sd(sF+jj-1) = sCurForce_sd(sF+jj-1) * windowSize / (jj+halfW);
%                         sCurForce_sd(eF-jj+1) = sCurForce_sd(eF-jj+1) * windowSize / jj;
                    end
                    
                    sCurForce_sd(isnan(curSlave))=NaN;
                    tracksNA(ii).(slaveSource) = sCurForce_sd;
    %                 sCurForce = [NaN(1,numNan) sCurForce];
                    bkgMaxForce = max(0,nanmax(sCurForce_sd(sF10before:sF5before))); % 10 is the tolr value in L1-reg
                    firstIncreaseTimeForce = find(sCurForce_sd>bkgMaxForce & 1:length(sCurForce_sd)>sF5before,1);
                else
                    bkgMaxForce = NaN;
                    firstIncreaseTimeForce = NaN;
                end
            else
                bkgMaxForce = max(0,max(curTrack.(slaveSource)(sF10before:sF5before)));
                try
                    firstIncreaseTimeForce = find(curTrack.(slaveSource)>bkgMaxForce & 1:nTime>sF5before,1);
                catch
                    nTime = length(curTrack.(slaveSource));
                    firstIncreaseTimeForce = find(curTrack.(slaveSource)>bkgMaxForce & 1:nTime>sF5before,1);
                end
            end
            if isempty(firstIncreaseTimeForce) || firstIncreaseTimeForce>curTrack.endingFrameExtraExtra
                bkgMaxIntAll(ii) = bkgMaxInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
            else
                forceTransmittingAll(ii) = true;
                firstIncreaseTimeIntAll(ii) = firstIncreaseTimeInt*tInterval; % in sec
                firstIncreaseTimeSlaveAll(ii) = firstIncreaseTimeForce*tInterval;
                bkgMaxIntAll(ii) = bkgMaxInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
                firstIncreseTimeIntAgainstSlaveAll(ii)=firstIncreaseTimeInt*tInterval - firstIncreaseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
            end
        else
            bkgMaxIntAll(ii) = bkgMaxInt;
        end
    else
        disp(''); %disp('Amplitude not increasing')
    end
    if nargout>6
        tracksNA(ii).forceTransmitting=forceTransmittingAll(ii);
        tracksNA(ii).firstIncreaseTimeInt=firstIncreaseTimeIntAll(ii); % in sec
        tracksNA(ii).firstIncreaseTimeForce = firstIncreaseTimeSlaveAll(ii);
        tracksNA(ii).bkgMaxInt(ii) = bkgMaxIntAll(ii);
        tracksNA(ii).bkgMaxSlave(ii) = bkgMaxSlaveAll(ii);
        tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeIntAgainstSlaveAll(ii);
    end
end
