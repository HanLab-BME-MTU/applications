function [assemRate,bestModel,bestSummary,yEntire,tRangeSelected] = ...
    getAssemRateFromTracks(curTrack,tIntervalMin, whichComp, minLifeTime,plotFit,tRangeSelected)
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
    whichComp='ampTotal';
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
if nargin<6
    tRangeSelected=tRange;
elseif isempty(tRangeSelected)
    assemRate = NaN;
    bestModel=[]; bestSummary=[]; yEntire=NaN;
    return
end


try
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
    if tRange(end) > length(curTrack.amp) %Just in case amp2 is well-fit in the tRange than amp, we stick with 'catch' section.
        tRange = find(~isnan(curTrack.amp));
        curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
    end
catch
    tRange = find(~isnan(curTrack.amp));
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
%     curAmpTotal =  curTrack.amp(tRange);
end

% There is a case where the first element is negative. In that case, we
% elevate the entire sereis above zero.
if curAmpTotal(1)<0
    curAmpTotal = curAmpTotal - curAmpTotal(1) + 0.1*std(curAmpTotal);
end
%% GET IT
if length(tRange)<minLifeTime+1
    assemRate=NaN;
    bestModel=[]; bestSummary=[]; yEntire=NaN;
    return
end
if strcmp(whichComp,'ampTotal') || strcmp(whichComp,'amp')
    [assemRate,bestModel,bestSummary,tRangeSelected] = getAssemRate(tIntervalMin*(tRange),curAmpTotal,minLifeTime);
elseif strcmp(whichComp,'ampTotal2') || strcmp(whichComp,'amp2')
    [assemRate,bestModel,bestSummary] = getAssemRateDirect(tIntervalMin*(tRange),curAmpTotal,tRangeSelected);
    tRangeSelected=[];
end
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
