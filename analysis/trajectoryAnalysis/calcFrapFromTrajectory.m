function [tOneHalf, recoveryCurve] = calcFrapFromTrajectory(trajectoryData,verbose)
%CALCFRAPFROMTRAJECTROY attempts to calculate FRAP statistics from a single MT trajectory
%
% SYNOPSIS [tOneHalf, recoveryCurve] = calcFrapFromTrajectory(trajectoryData,verbose)
%
% INPUT    trajectoryData(1:n) : output of calculateTrajectoryFromIdlist
%          verbose             : (opt) [{0}/1] plot of recovery
%          
% OUTPUT   tOneHalf      : (nx1) t1/2 of recovery. Inf if 50% recovery is never
%                             reached
%          recoveryCurve : (nx1) struct with fields 
%                     .recovery : individual recovery measurements
%                     .time     : corresponding time
%
% could easily be modified to accept nanLists from simulations
%
% c: 04/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=====================
%  TEST INPUT
%=====================

% defaults
% verbose: 0

if isempty(trajectoryData) | ~isfield(trajectoryData,'distance') | ~isfield(trajectoryData,'time') | ~isfield(trajectoryData,'timePoints')
    error('CALCFRAPFROMTRAJECTROY needs a non-empty trajectoryData with the correct fields. See CALCULATETRAJECTORYFROMIDLIST for details')
end

if nargin < 2 | isempty(verbose)
    verbose = 0;
end

%-------------------


%===============================
% LOOP & FRAP
%===============================

numData = length(trajectoryData);
tOneHalf = repmat(NaN,numData,1);

for iData = 1:numData
    
    % INIT DATA
    
    %-----make nan-separated list of distance, add 'corrected' distance at
    %-----end. Throw away leading deleted frames
    
    % throw away number of leading NaNs 
    minTimePoint = trajectoryData(iData).timePoints(1)-1;
    timePoints = trajectoryData(iData).timePoints - minTimePoint;
    
    trajectoryLength = timePoints(end);
    % init distList
    distanceList = repmat(NaN,[2*trajectoryLength,1]);
    
    % first distances
    distanceList(timePoints) = trajectoryData(iData).distance(:,1);
    
    % calc correction: difference of mean of first and last five data
    % points
    diffDist = robustMean(trajectoryData(iData).distance(1:5,1)) -...
        robustMean(trajectoryData(iData).distance(end-4:end,1));
    
    % second distances
    distanceList(timePoints + trajectoryLength) = ...
        trajectoryData(iData).distance(:,1) + diffDist;
    
    %---------
    
    % find mean time by finding where we had just one timepoint
    % difference
    timeAndTp = [trajectoryData(iData).time,trajectoryData(iData).timePoints];
    deltaTaT = diff(timeAndTp,1,1);
    timeInterval = mean(deltaTaT(find(deltaTaT(:,3)==1),1));
    
    %-----calculate initial length: extrapolate forward (use KJ's algorithm
    %-----in the future).
    % L0: shortest distance of the mt over lifetime, init as initial length
    distZero = repeatEntries(trajectoryData(iData).distance(:,1),[diff(timePoints);1]);
    % distInit : for the calculation we want to know the total initial
    % tubulin fluorescence, so that we can compare the number of tubulin
    % subunits into the MT
    distInitSum = sum(distZero);
    
    currentTime = [1:trajectoryLength]';
    
    % RECOVERY
    
    for t = 1:trajectoryLength-1
        
        % find good entries in distanceList
        currentDistance = distanceList(currentTime);
        
        goodTimePoints  = find(isfinite(currentDistance));
        goodDistance    = currentDistance(goodTimePoints);
        
        % update L0
        distZero(goodTimePoints) = min(distZero(goodTimePoints),goodDistance);
        
        % calculate FRAP with good points: (new fluorescent tubulin)/(preFrap fluorescent tubulin) 
        % sum(Lt-L0)/sum(Linit)
        recoveryCurve(iData).recovery(t) = sum(goodDistance - distZero(goodTimePoints))/distInitSum;
        recoveryCurve(iData).time(t) = (t-1)*timeInterval;
        
        % next current time
        currentTime = currentTime + 1;
        
    end %for t = 1:trajectoryLength-1
    
    if verbose
        figure,plot(recoveryCurve(iData).time,recoveryCurve(iData).recovery)
    end
    
    rec05 = find(recoveryCurve(iData).recovery > 0.5);
    if ~isempty(rec05)
        tOneHalf(iData) = recoveryCurve(iData).time(rec05(1));
    end
    
end %for iData = 1:numData

