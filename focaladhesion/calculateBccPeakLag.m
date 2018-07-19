function [ tLagBccPosG2,tLagBccNegG2,BccInLaterG2,BccInPriorG2 ] = calculateBccPeakLag( ts1, ts2 , tRef, tInterval)
%function [ tLagBccPosG2,tLagBccNegG2,BccInLaterG2,BccInPriorG2 ] =
%calculateBccPeakLag( ts1, ts2 , SFE) calculates time delay in cross
%variance time series between ts1 and ts2 against frame of reference
%(tRef).
%   Sangyoon Han
if length(ts2)<length(ts1)
    tempForce=ts2;
    ts2=ts1;
    ts2(1:length(tempForce))=tempForce;
elseif length(ts1)<length(ts2)
    tempForce=ts1;
    ts1=ts2;
    ts1(1:length(tempForce))=tempForce;
end
curBcc = crossVariance(ts1,ts2,11);
%Will find the peak against curTrack.startingFrameExtra to the
%front and back
maxBcc=nanmax(curBcc);
indCands=locmax1d(curBcc,5);
% Exclude the locmaxes at the ends
firstInd = find(~isnan(curBcc),1);
lastInd = find(~isnan(curBcc),1,'last');
indCands(ismember(indCands,[firstInd,lastInd]))=[];
% if the loc max is close to curTrack.startingFrameExtra and the
% actual value is above 0.8*maxBcc, we'll admit the time as
% BccShift time

% Sort based on proximity to curTrack.startingFrameExtra
%         SFE = curTrack.startingFrameExtra;
frameDiff = indCands-tRef;
tLagBccNegG2 = NaN;
BccInPriorG2 = NaN;
tLagBccPosG2 = NaN;
BccInLaterG2 = NaN;

if any(frameDiff<0)
    indNegDiff= indCands(frameDiff<0);
    for curIndCand = indNegDiff'
        if curBcc(curIndCand)>0.3*maxBcc
            tLagBccNegG2 = frameDiff(indCands==curIndCand)*tInterval;
            BccInPriorG2 = curBcc(curIndCand);
            break
        end
    end
end
indPosDiff= indCands(frameDiff>=0);

for curIndCand = indPosDiff'
    if curBcc(curIndCand)>0.3*maxBcc
        tLagBccPosG2 = frameDiff(indCands==curIndCand)*tInterval;
        BccInLaterG2 = curBcc(curIndCand);
        break
    end
end

end

