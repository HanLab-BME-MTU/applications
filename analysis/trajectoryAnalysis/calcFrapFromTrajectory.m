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

if nargin == 0 || isempty(trajectoryData)
    run = loadRunsFromFile(1);
    trajectoryData = run(1).data;
elseif ~isfield(trajectoryData,'distance') || ~isfield(trajectoryData,'time') || ~isfield(trajectoryData,'timePoints')
    error('CALCFRAPFROMTRAJECTROY needs a non-empty trajectoryData with the correct fields. See CALCULATETRAJECTORYFROMIDLIST for details')
end

if nargin < 2 || isempty(verbose)
    verbose = 0;
end

%-------------------


%===============================
% LOOP & FRAP
%===============================

% init vars
numData = length(trajectoryData);
tOneHalf = repmat(NaN,numData,3);
recoveryCurve(1:numData,1) = struct('time',[],'recovery',[]);
options = optimset('Display','off');

% find the first timePoint
firstTimePoint = catStruct(1,'trajectoryData.timePoints(1)');

% find first and last five distance measurements (5-by-numData)
firstFiveDist = catStruct(2,'trajectoryData.distance(1:5,1)');
lastFiveDist  = catStruct(2,'trajectoryData.distance(end-4:end,1)');

% convert trajectories to nan-separated lists
nanTrajectories = convertTrajectoryData(trajectoryData,0);

for iData = 1:numData
    
    % INIT DATA
    
    % create full distanceList. Throw away leading nans, then catenate the
    % trajectory with a copy of itself which has been corrected so that
    % there will not be a severe jump in distance at the end. 
    
    % read distances
    distanceList = nanTrajectories(iData).distance(firstTimePoint(iData):end,1);
    
    % get timePoints where we have not a NaN
    timePoints = find(~isnan(distanceList));
    trajectoryLength = timePoints(end);
    
    % calculate difference of distances (afterwards, we will add this to
    % the second trajectory. Ideally, d(1)+diffDist=d(end)
    diffDist = robustMean(lastFiveDist(:,iData)-firstFiveDist(:,iData));
    
    % catenate distanceList
    distanceList = [distanceList;distanceList+diffDist];
    
    
        % correct for drift - subtract the drift, but conserve the mean
    % B = A*X
    % we need to go back to not-nans, so let's read make a goodIdx. It is
    % still pretty fast and less error-prone than playing with timePoints
    
    goodIdx = find(~isnan(distanceList));
    goodLength = length(goodIdx);
    distListLength = length(distanceList);
    
    A = ones(distListLength,2); %allocate A
    A(:,2) = [1:distListLength]'; % fill in time
    B = distanceList; % fill in distance
    %diagonal elements of the covariance matrix
    sigmaDist = repmat(trajectoryData(iData).distance(:,2),[2,1]);
    V = diag([sigmaDist].^2); 

    %calculate linear fit. a0+a1*x=y
    [X] = myLscov(A(goodIdx,:),B(goodIdx),V);
    
    % calculate mean
    meanDist = weightedStats(B(goodIdx),sigmaDist,'s');
    
    % correct distanceList
    distanceList = distanceList - A*X + meanDist;

    
    
    % read timeInterval
    timeInterval = nanTrajectories(iData).timeInterval;
    
    %-----calculate initial length: extrapolate forward (use KJ's algorithm
    %-----in the future).
    % L0: shortest distance of the mt over lifetime, init as initial length
    
    % distZero with interpolation
    % distZero = repeatEntries(distanceList(goodIdx(1:goodLength/2)),[diff(timePoints);1]);
    % goodZero = ones(distListLength/2,1);
    
    distZero = distanceList(1:distListLength/2);
    goodZero = ~isnan(distZero);
    
    
    
    
    
    % distInit : for the calculation we want to know the total initial
    % tubulin fluorescence, so that we can compare the number of tubulin
    % subunits in the MT
    distInitSum = nansum(distZero);
    
    timeZero = [1:trajectoryLength]'-1;
    
    
    % RECOVERY
    
    % init
    recoveryCurve(iData).time = timeZero * timeInterval;
    recoveryCurve(iData).recovery = zeros(trajectoryLength,1);
    
    for t = timeZero'+1
        
        % find good entries in distanceList
        currentDistance = distanceList(timeZero+t);
        
        goodTimePoints  = find(isfinite(currentDistance) & goodZero);
        goodDistance    = currentDistance(goodTimePoints);
        
        % update L0
        distZero(goodTimePoints) = min(distZero(goodTimePoints),goodDistance);
        
        % calculate FRAP with good points: (new fluorescent tubulin)/(preFrap fluorescent tubulin) 
        % sum(Lt-L0)/sum(Linit)
        recoveryCurve(iData).recovery(t) = sum(goodDistance - distZero(goodTimePoints))/distInitSum;        
        
    end %for t = 1:trajectoryLength-1
    
    if verbose
        figure,plot(recoveryCurve(iData).time,recoveryCurve(iData).recovery)
    end
    
    % calculate recovery: y=a(1-exp(-x/tau))
    
    % prepare estimates
    u0(1) = recoveryCurve(iData).recovery(end); 
    biggerIdx = find(recoveryCurve(iData).recovery > u0(1)*exp(-1));
    if ~isempty(biggerIdx)
        u0(2) = biggerIdx(1) * timeInterval;


        tOneHalf(iData,:) = [lsqnonlin(inline('y-u(1)*(1-exp(-x/u(2)))','u','x','y'),...
            u0,[],[],options,recoveryCurve(iData).time,recoveryCurve(iData).recovery),trajectoryLength];
    else
        tOneHalf(iData,3) = trajectoryLength;
    end
    
end %for iData = 1:numData


% calculate FRAP for averaged recovery
[maxTp, maxTpIdx] = max(tOneHalf(:,3));
allCurves = repmat(NaN,maxTp,numData);
for iData = 1:numData
    allCurves(1:tOneHalf(iData,3),iData) = recoveryCurve(iData).recovery;
end

globalRecovery = zeros(maxTp,1);
for t=1:maxTp
    goodIdx = find(~isnan(allCurves(t,:)));
    globalRecovery(t) = weightedStats(allCurves(t,goodIdx),tOneHalf(goodIdx,3));
end

% fit
u0(1) = globalRecovery(end);
biggerIdx = find(globalRecovery > u0(1)*exp(-1));
if ~isempty(biggerIdx)
    u0(2) = biggerIdx(1) * timeInterval;


    globalT = [lsqnonlin(inline('y-u(1)*(1-exp(-x/u(2)))','u','x','y'),...
        u0,[],[],options,recoveryCurve(maxTpIdx(1)).time,globalRecovery)];
end

recoveryCurve(1).global = struct('curve',globalRecovery,'parms',globalT);

