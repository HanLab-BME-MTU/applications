function [firstIncreseTimeIntAgainstSlaveAll,forceTransmittingAll...
    ,firstIncreseTimeIntAll,firstIncreseTimeSlaveAll,bkgMaxIntAll,bkgMaxSlaveAll] ...
    = calculateFirstIncreaseTimeTracks(tracksNA,splineParamInit,preDetecFactor,tInterval,varargin)
% [tracksNA,firstIncreseTimeIntAgainstForceAll] =
%  calculateFirstIncreaseTimeTracks(tracksNA,splineParamInit,preDetecFactor,tInterval)
%  calculates the time lag of the main intensity (ampTotal) against the slave source.
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
ip.addOptional('splineParamInit',0.99,@isscalar)
ip.addOptional('preDetecFactor',1/5,@(x)isscalar(x))
ip.addOptional('tInterval',1,@(x)isscalar(x))
ip.addOptional('plotEachTrack',false,@(x)islogical(x)||isempty(x))
ip.addParamValue('slaveSource','forceMag',@(x)ismember(x,{'forceMag','ampTotal2','ampTotal3'})); % collect NA tracks that ever close to cell edge
ip.parse(tracksNA,splineParamInit,preDetecFactor,tInterval,varargin{:});
slaveSource=ip.Results.slaveSource;
    useSmoothing=false;
    if splineParamInit<1
        useSmoothing=true;
    end
    firstIncreseTimeIntAgainstSlaveAll=NaN(numel(tracksNA),1);
    forceTransmittingAll=false(numel(tracksNA),1);
    firstIncreseTimeIntAll=NaN(numel(tracksNA),1);
    firstIncreseTimeSlaveAll=NaN(numel(tracksNA),1);
    bkgMaxIntAll=NaN(numel(tracksNA),1);
    bkgMaxSlaveAll=NaN(numel(tracksNA),1);
    for ii=1:numel(tracksNA)
        curTrack = tracksNA(ii);
        sFEE = curTrack.startingFrameExtraExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
        % See how many frames you have before the startingFrameExtra
        numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
        numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
        numPreSigStart = min([20,numFramesBefore 3*numPreFrames]);
        sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
        sF10before = max(sFEE,curTrack.startingFrameExtra-3*numPreFrames);
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

            bkgMaxInt = nanmax(sd(sF10before:sF5before));
            firstIncreaseTimeInt = find(sd>bkgMaxInt & 1:length(sd)>sF5before,1);
        else
            bkgMaxInt = max(curTrack.ampTotal(sF10before:sF5before));
            firstIncreaseTimeInt = find(curTrack.ampTotal>bkgMaxInt & 1:nTime>sF5before,1);
        end
%         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstIncreaseTimeInt)
            if useSmoothing
                curSlave=d;
                curSlave(curTrack.startingFrameExtraExtra:curTrack.endingFrameExtraExtra) = ...
                 getfield(curTrack,{1},slaveSource,{curTrack.startingFrameExtraExtra:curTrack.endingFrameExtraExtra});
                %                 curForce(isnan(curForce)) = [];
                if sum(~isnan(curSlave))>1
                    sCurForce_spline= csaps(tRange,curSlave,splineParamInit);
                    sCurForce_sd=ppval(sCurForce_spline,tRange);
                    sCurForce_sd(isnan(curSlave))=NaN;
    %                 sCurForce = [NaN(1,numNan) sCurForce];
                    bkgMaxForce = max(10,nanmax(sCurForce_sd(sF10before:sF5before))); % 10 is the tolr value in L1-reg
                    firstIncreaseTimeForce = find(sCurForce_sd>bkgMaxForce & 1:length(sCurForce_sd)>sF5before,1);
                else
                    bkgMaxForce = NaN;
                    firstIncreaseTimeForce = NaN;
                end
            else
                bkgMaxForce = max(10,max(curTrack.forceMag(sF10before:sF5before)));
                firstIncreaseTimeForce = find(getfield(curTrack,slaveSource)>bkgMaxForce & 1:nTime>sF5before,1);
            end
            if isempty(firstIncreaseTimeForce) || firstIncreaseTimeForce>curTrack.endingFrameExtraExtra
                bkgMaxIntAll(ii) = bkgMaxInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
           else
                forceTransmittingAll(ii) = true;
                firstIncreseTimeIntAll(ii) = firstIncreaseTimeInt*tInterval; % in sec
                firstIncreseTimeSlaveAll(ii) = firstIncreaseTimeForce*tInterval;
                bkgMaxIntAll(ii) = bkgMaxInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
                firstIncreseTimeIntAgainstSlaveAll(ii)=firstIncreaseTimeInt*tInterval - firstIncreaseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
            end
        else
            bkgMaxIntAll(ii) = bkgMaxInt;
        end
    end
