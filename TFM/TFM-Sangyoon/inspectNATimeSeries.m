function [BccMatrix,avgBcc,tStartingFrame,indexValidBccInTracks,firstTimeAboveZeroBccAll,...
    firstTimeAboveHalfBccAll,firstTimeAboveOneBccAll] = inspectNATimeSeries(tracksNA,tInterval,varargin)
% [tracksNA,firstIncreseTimeIntAgainstForceAll] =
%  inspectNATimeSeries(tracksNA,splineParamInit,preDetecFactor,tInterval)
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
ip.addOptional('tInterval',1,@(x)isscalar(x))
ip.addParameter('splineParamInit',0.99,@isscalar)
ip.addParameter('preDetecFactor',1/5,@(x)isscalar(x))
ip.addParameter('slaveSource','forceMag',@(x)ismember(x,{'forceMag','ampTotal2','ampTotal3'})); % collect NA tracks that ever close to cell edge
ip.parse(tracksNA,tInterval,varargin{:});
splineParamInit=ip.Results.splineParamInit;
preDetecFactor=ip.Results.preDetecFactor;
slaveSource=ip.Results.slaveSource;

if isempty(tracksNA)
    BccMatrix = [];
    avgBcc = [];
    tStartingFrame = NaN;
    indexValidBccInTracks = [];
    firstTimeAboveZeroBccAll = [];
    firstTimeAboveHalfBccAll = [];
    firstTimeAboveOneBccAll = [];
    return
end
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
firstTimeAboveZeroBccAll=NaN(numel(tracksNA),1);
firstTimeAboveHalfBccAll=NaN(numel(tracksNA),1);
firstTimeAboveOneBccAll=NaN(numel(tracksNA),1);

iii=0;
tFluc = 10; % this is actually in frame
tStartingFrame = tFluc/2 + 1; %30 + 
for ii=1:numel(tracksNA)
    curTrack = tracksNA(ii);
    curEarlyAmpSlope = curTrack.earlyAmpSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
%     curAmpSlope = curTrack.ampSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
%     curForceSlope = curTrack.earlyAmpSlope; if isnan(curEarlyAmpSlope); curEarlyAmpSlope=-1000; end
    try
        [~,curAmpSlope] = regression((1:curTrack.lifeTime+1),curTrack.amp(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
    catch
        % There was a case where  startingFrameExtra and endingFrameExtra
        % are incorrectly recorded. 
        curTrack.startingFrameExtra = find(~isnan(curTrack.amp),1);
        curTrack.endingFrameExtra = find(~isnan(curTrack.amp),1,'last');
        curTrack.lifeTime = curTrack.endingFrameExtra -curTrack.startingFrameExtra;
        [~,curAmpSlope] = regression((1:curTrack.lifeTime+1),curTrack.amp(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
    end
%     [~,curForceSlope] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
%     disp(['curAmpSlope = ' num2str(curAmpSlope) '  curEarlyAmpSlope = ' num2str(curEarlyAmpSlope)])

    sFEE = curTrack.startingFrameExtraExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
    % See how many frames you have before the startingFrameExtra
    numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
    numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
    numPreSigStart = min([20,numFramesBefore 3*numPreFrames]);
    sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
    sF10before = max(sFEE,curTrack.startingFrameExtra-3*numPreFrames);

    if (curAmpSlope>0 || curEarlyAmpSlope>0) %&& (numFramesBefore>1) % && curForceSlope>0 % intensity should increase. We are not 
        % interested in decreasing intensity which will 
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
                 curTrack.(slaveSource)(curTrack.startingFrameExtraExtra:curTrack.endingFrameExtraExtra);
                %                 curForce(isnan(curForce)) = [];
                if sum(~isnan(curSlave))>1
                    sCurForce_spline= csaps(tRange,curSlave,splineParamInit);
                    sCurForce_sd=ppval(sCurForce_spline,tRange);
                    sCurForce_sd(isnan(curSlave))=NaN;
    %                 sCurForce = [NaN(1,numNan) sCurForce];
                    bkgMaxForce = max(0,nanmax(sCurForce_sd(sF10before:sF5before))); % 10 is the tolr value in L1-reg
                    firstIncreaseTimeForce = find(sCurForce_sd>bkgMaxForce & 1:length(sCurForce_sd)>sF5before,1);
                else
                    bkgMaxForce = NaN;
                    firstIncreaseTimeForce = NaN;
                end
            else
                bkgMaxForce = max(0,max(curTrack.forceMag(sF10before:sF5before)));
                firstIncreaseTimeForce = find(curTrack.(slaveSource)>bkgMaxForce & 1:nTime>sF5before,1);
            end
%             if isempty(firstIncreaseTimeForce) || firstIncreaseTimeForce>curTrack.endingFrameExtraExtra
%                 bkgMaxIntAll(ii) = bkgMaxInt;
%                 bkgMaxSlaveAll(ii) = bkgMaxForce;
%                 
% %                 h=figure;  subplot(2,1,1), plot(sd,'o-'), hold on, plot(tRange(firstIncreaseTimeInt),sd(firstIncreaseTimeInt),'ro'); 
% %                 subplot(2,1,2), plot(curSlave,'o-') 
% %                 uiwait(); 
%             else
                iii=iii+1;
                forceTransmittingAll(ii) = true;
                firstIncreseTimeIntAll(ii) = firstIncreaseTimeInt*tInterval; % in sec
%                 firstIncreseTimeSlaveAll(ii) = firstIncreaseTimeForce*tInterval;
                bkgMaxIntAll(ii) = bkgMaxInt;
                bkgMaxSlaveAll(ii) = bkgMaxForce;
%                 firstIncreseTimeIntAgainstSlaveAll(ii)=firstIncreaseTimeInt*tInterval - firstIncreaseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
                
                % Bcc (co-fluctuation) calculation: Here the rule is that
                % real starting frame is from 36 which is the maximum frame
                % difference between startingFrameExtra and
                % startingFrameExtraExtra.
%                 lastFrameCC = curTrack.endingFrameExtraExtra-tFluc; clear Bcc
%                 tShift = tStartingFrame +curTrack.startingFrameExtraExtra-curTrack.startingFrameExtra;
%                 for jj=curTrack.startingFrameExtraExtra:lastFrameCC
%                     sd_segment = sd(jj:jj+tFluc); avgSD = mean(sd_segment);
%                     force_seg = sCurForce_sd(jj:jj+tFluc); avgForceSeg = mean(force_seg);
%                     sigma2cc=1/tFluc*sum((sd_segment-avgSD).*(force_seg-avgForceSeg));
%                     Bcc(jj-curTrack.startingFrameExtraExtra+tShift) = sigma2cc/sqrt(avgSD*avgForceSeg);
%                 end
%                 if tShift>1
%                     Bcc(1:tShift-1)=NaN; 
%                 end
                Bcc = crossVariance(sd,sCurForce_sd,tFluc);
                
                BccAll{iii}=Bcc;
                indexValidBccInTracks(iii) = ii;
                if ~isempty(find(Bcc > 0, 1))
                    firstTimeAboveZeroBccAll(ii) = find(Bcc > 0, 1);
                    if ~isempty(find(Bcc > 0.5, 1))
                        firstTimeAboveHalfBccAll(ii) = find(Bcc > 0.5, 1);
                        if ~isempty(find(Bcc > 1.0, 1))
                            firstTimeAboveOneBccAll(ii) = find(Bcc > 1.0, 1);
                        end
                    end
                end
%             end
        else
            bkgMaxIntAll(ii) = bkgMaxInt;
%             h=figure; subplot(2,2,1), plot(sd,'o-')
%             subplot(2,2,3), plot(curSlave,'o-')
%             subplot(2,2,2), plot(curTrack.xCoord,'o-'); hold on; 
%             subplot(2,2,4), plot(curTrack.yCoord,'o-')
%             uiwait(); 
        end
    else
%         h=figure; subplot(2,2,1), plot(curTrack.ampTotal,'o-')
%         subplot(2,2,3), plot(curTrack.(slaveSource),'o-')
%         subplot(2,2,2), plot(curTrack.xCoord,'o-'); hold on; 
%         subplot(2,2,4), plot(curTrack.yCoord,'o-')
%         uiwait();
    end
end
if exist('BccAll','var')
    numConditions = numel(BccAll);
    [lengthLongest]=max(cellfun(@(x) length(x),BccAll));
    BccMatrix = NaN(numConditions,lengthLongest);
    for k=1:numConditions
        BccMatrix(k,1:length(BccAll{k})) = BccAll{k};
    end
%     disp(['The BccMatrix is aligned with ' num2str(iRepTal2) ' from first.'])

    avgBcc = nanmean(BccMatrix,1);
else
    BccMatrix = [];
    avgBcc = [];
    tStartingFrame = NaN;
    indexValidBccInTracks = [];
    firstTimeAboveZeroBccAll = [];
    firstTimeAboveHalfBccAll = [];
    firstTimeAboveOneBccAll = [];
end
