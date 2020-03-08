function [assemRate,bestModel,bestSummary,yEntire] = getAssemRateFromTracks(curTrack,tIntervalMin, whichComp, minLifeTime,plotFit)
%function assemRate = getAssemRateFromTracks(curTrack) is a wrapper
%function for tracksNA to run getAssemRate 
% input:
%       curTrack    tracksNA struct variable which contains amp, ampTotal
%                   or ampTotal2 or ampTotal3
%       tIntervalMin 
% output:
%       assemRate:  assembly rate calculated as first rate constant of
%                   logarithm of normalized amplitude (from amp or
%                   ampTotal)
% Sangyoon Han 2020/3/4
if nargin<3
    plotFit=false;
    minLifeTime=15;
    whichComp='amp';
end
if nargin<4
    plotFit=false;
    minLifeTime=15;
end
if nargin<5
    plotFit=false;
end
%% Setting
tRange = curTrack.startingFrameExtra:curTrack.endingFrameExtra;
% curAmpTotal =  curTrack.ampTotal(tRange);
try
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
catch
    tRange = find(~isnan(curTrack.amp));
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
%     curAmpTotal =  curTrack.amp(tRange);
end

%% GET IT
if length(tRange)<minLifeTime+1
    assemRate=NaN;
    bestModel=[]; bestSummary=[]; yEntire=NaN;
    return
end
[assemRate,bestModel,bestSummary] = getAssemRate(tIntervalMin*(tRange),curAmpTotal,minLifeTime);
yEntire = curAmpTotal/curAmpTotal(1);


%% To look at the fit from entire TS
if plotFit
    subplot(2,1,1)
    plot(tIntervalMin*(tRange),yEntire, 'o-.')
    hold on
    plot(bestModel.Variables.x1,exp(bestModel.Fitted), 'r-')
    subplot(2,1,2)
    plot(tIntervalMin*(tRange),log(yEntire), 'o-.')
    hold on
    plot(bestModel.Variables.x1,bestModel.Fitted, 'r-')
end
end
