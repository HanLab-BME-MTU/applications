function [endingTimeIntAgainstForceAll,endingTimeIntAll,endingTimeSlaveAll] ...
    = calculateEndingTimeLagFromTracks(tracksNA,splineParam,preDetecFactor,tInterval,varargin)
% [endingTimeIntAgainstForceAll,endingTimeIntAll,endingTimeSlaveAll] ...
%     = calculateEndingTimeLagFromTracks(tracksNA,splineParam,tInterval,varargin)% calculates the time lag of the main intensity ending (ampTotal) against 
% ending of the slave source.
% input:
%       splineParam: smoothing parameter (0-1). Use 1 if you don't want to
%       smooth the signal.
%       'slaveSource': either 'forceMag', 'ampTotal2' or 'ampTotal3'
% output:
%       
% Sangyoon Han Oct 2017
ip =inputParser;
ip.addRequired('tracksNA',@isstruct)
ip.addOptional('splineParam',0.99,@isscalar)
ip.addOptional('preDetecFactor',1/5,@(x)isscalar(x))
ip.addOptional('tInterval',1,@(x)isscalar(x))
ip.addParamValue('slaveSource','forceMag',@(x)ismember(x,{'forceMag','ampTotal2','ampTotal3'})); % collect NA tracks that ever close to cell edge
ip.parse(tracksNA,splineParam,preDetecFactor,tInterval,varargin{:});
slaveSource=ip.Results.slaveSource;

useSmoothing=false;
if splineParam<1
    useSmoothing=true;
end
endingTimeIntAgainstForceAll=NaN(numel(tracksNA),1);
disassemblingAll=false(numel(tracksNA),1);
endingTimeIntAll=NaN(numel(tracksNA),1);
endingTimeSlaveAll=NaN(numel(tracksNA),1);
bkgMaxIntAll=NaN(numel(tracksNA),1);
bkgMaxSlaveAll=NaN(numel(tracksNA),1);
%% Calculation
for ii=1:numel(tracksNA)
    curTrack = tracksNA(ii);
    sFEE = curTrack.startingFrameExtraExtra;
    eFEE = curTrack.endingFrameExtraExtra;
    eFE = curTrack.endingFrameExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
    % See how many frames you have after the endingFrameExtra
    numFramesAfter = eFEE - curTrack.endingFrameExtra;
    numPostFrames = max(1,floor(preDetecFactor*numFramesAfter));
    eF10before = max(min(eFE,eFE-3*numPostFrames),1);
    eF5after = eFE;
    d = tracksNA(ii).ampTotal;
    nTime = length(d);
    if eFE == eFEE % In this case this adhesion is not completely disassembled
        disassemblingAll(ii)=false;
    else
        if useSmoothing
    %             d = tracksNA(ii).ampTotal;
            tRange = tracksNA(ii).iFrame;
            d(d==0)=NaN;
            warning('off','SPLINES:CHCKXYWP:NaNs')
            try
                sd_spline= csaps(tRange,d,splineParam);
            catch
                d = tracksNA(ii).amp;
                d(sFEE:eFEE) = tracksNA(ii).ampTotal(sFEE:eFEE);
                sd_spline= csaps(tRange,d,splineParam);
            end
            sd=ppval(sd_spline,tRange);
            sd(isnan(d))=NaN;

            endingSigMinInt = nanmin(sd(eF10before:eF5after));
            firstDecreseTimeInt = find(sd<endingSigMinInt & 1:length(sd)>eF10before,1);
        else
            endingSigMinInt = max(curTrack.ampTotal(eF10before:eF5after));
            firstDecreseTimeInt = find(curTrack.ampTotal<endingSigMinInt & 1:nTime>eF10before,1);
        end
    %         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstDecreseTimeInt)
            if useSmoothing
                curSlave=d;
                curSlave(sFEE:eFEE) = ...
                 getfield(curTrack,{1},slaveSource,{sFEE:eFEE});
                %                 curForce(isnan(curForce)) = [];
                sCurForce_spline= csaps(tRange,curSlave,splineParam);
                sCurForce_sd=ppval(sCurForce_spline,tRange);
                sCurForce_sd(isnan(curSlave))=NaN;
    %                 sCurForce = [NaN(1,numNan) sCurForce];
                bkgMaxForce = nanmin(sCurForce_sd(eF10before:eF5after));
                firstDecreaseTimeForce = find(sCurForce_sd<bkgMaxForce & 1:length(sCurForce_sd)>eF10before,1);
            else
                bkgMaxForce = max(curTrack.forceMag(eF10before:eF5after));
                firstDecreaseTimeForce = find(getfield(curTrack,slaveSource)<bkgMaxForce & 1:nTime>eF10before,1);
            end
            if isempty(firstDecreaseTimeForce) || firstDecreaseTimeForce>eFEE
                bkgMaxIntAll(ii) = endingSigMinInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
           else
                disassemblingAll(ii) = true;
                endingTimeIntAll(ii) = firstDecreseTimeInt*tInterval; % in sec
                endingTimeSlaveAll(ii) = firstDecreaseTimeForce*tInterval;
                bkgMaxIntAll(ii) = endingSigMinInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
                endingTimeIntAgainstForceAll(ii)=firstDecreseTimeInt*tInterval - firstDecreaseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
            end
        else
            bkgMaxIntAll(ii) = endingSigMinInt;
        end
    end
end


