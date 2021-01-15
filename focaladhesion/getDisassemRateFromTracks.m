function [disassemRate,bestModel,bestSummary,yEntire,tRangeSelected] = ...
    getDisassemRateFromTracks(curTrack,tIntervalMin, whichComp, minLifeTime,plotFit,tRangeSelected)
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
    minLifeTime=25;
    whichComp='ampTotal';
end
if nargin<4
    plotFit=false;
    minLifeTime=25;
end
if nargin<5
    plotFit=false;
end
%% Set up
tRange = curTrack.startingFrameExtra:curTrack.endingFrameExtra;
if nargin<6
    tRangeSelected=tRange;
end

% curAmpTotal =  curTrack.ampTotal(tRange);
try
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
    if tRange(end) > length(curTrack.amp) %Just in case amp2 is well-fit in the tRange than amp, we stick with 'catch' section.
        tRange = find(~isnan(curTrack.amp));
        curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
    end
catch
    tRange = find(~isnan(curTrack.amp));
    curAmpTotal =  getfield(curTrack,{1},whichComp,{tRange});
end
%% GET IT
if length(tRange)<minLifeTime+1
    disassemRate=NaN;
    bestModel=[]; bestSummary=[]; yEntire=NaN;
    return
end

% Now I am revising this. For amp2, the time points should come from the
% amp1 directly.
if strcmp(whichComp,'ampTotal') || strcmp(whichComp,'amp')
    [disassemRate,bestModel,bestSummary,tRangeSelected] = getDisassemRate(tIntervalMin*(tRange),curAmpTotal,minLifeTime);
    yEntire = bestSummary.startingTS./curAmpTotal;
elseif strcmp(whichComp,'ampTotal2') || strcmp(whichComp,'amp2')
    [disassemRate,bestModel,bestSummary] = getDisassemRateDirect(tIntervalMin*(tRange),curAmpTotal,tRangeSelected);
    yEntire = bestSummary.startingTS./curAmpTotal;
    tRangeSelected=[];
end

%% Plot
if plotFit
    % To look at the fit from entire TS
    hold off
    subplot(3,1,1)
    plot(tIntervalMin*(tRange),yEntire, 'o-.')
    hold on
    plot(bestModel.Variables.x1,exp(bestModel.Fitted), 'r-')
    subplot(3,1,2)
    hold off
    plot(tIntervalMin*(tRange),log(yEntire), 'o-.')
    hold on
    plot(bestModel.Variables.x1,bestModel.Fitted, 'r-')
    subplot(3,1,3)
    hold off
    plot(bestModel.Variables.x1,bestModel.Variables.y, 'o-.')
    hold on
    plot(bestModel.Variables.x1,bestModel.Fitted, 'r-')
end

end