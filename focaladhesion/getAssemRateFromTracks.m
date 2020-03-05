function [assemRate,bestModel,bestSummary,yEntire] = getAssemRateFromTracks(curTrack,tIntervalMin, plotFit)
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
end

tRange = curTrack.startingFrameExtra:curTrack.endingFrameExtra;
% curAmpTotal =  curTrack.ampTotal(tRange);
curAmpTotal =  curTrack.amp(tRange);

%% GET IT
[assemRate,bestModel,bestSummary] = getAssemRate(tIntervalMin*(tRange),curAmpTotal);
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
