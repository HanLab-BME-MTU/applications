function [H,pValue,errFlag] = rankTest(traj,significance)
%RANKTEST tests the hypothesis that a time series is IID by checking whether it has any linear trend.
%
%SYNOPSIS [H,pValue,errFlag] = rankTest(traj,significance)
%
%INPUT  traj        : Observations of time series to be tested. Either an 
%                     array of structures traj(1:nTraj).observations, or 
%                     a column representing one single trajectory. Missing 
%                     points should be indicated with NaN.
%       significance: Significance level of hypothesis test. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       pValue  : Probability of obtaining a test statistic >= observed
%                 value assuming that null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.37. It can be used to detect linear
%       trends in a time series.
%
%       MUST BE FIXED FOR MULTIPLE TRAJECTORIES WITH MISSING DATA POINTS!
%
%Khuloud Jaqaman, August 2004

%%% FIX STD %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = [];
pValue = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin < 1
    disp('--rankTest: You should at least input time series to be analyzed!');
    errFlag = 1;
    return
end

%check trajectory and turn into struct if necessary
if ~isstruct(traj)
    tmp = traj;
    clear traj
    traj.observations = tmp;
    clear tmp
elseif ~isfield(traj,'observations')
    disp('--rankTest: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%check that trajectories are column vectors and get their overall length
numTraj = length(traj);
for i=1:numTraj
    nCol = size(traj(i).observations,2);
    if nCol > 1
        disp('--rankTest: Each trajectory should be a column vector!');
        errFlag = 1;
    else
        trajLength(i) = length(traj(i).observations);
        availLength(i) = length(find(~isnan(traj(i).observations)));
    end
end

%assign default value of significance, if needed
if nargin < 2
    significance = 0.05;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate number of pairs where x(t1) > x(t2) for t1>t2 and get total 
%number of pairs tested
numIncPairs = 0;
numPairsTested = 0

for j = 1:numTraj
    for i=1:trajLength(j)-1

        %calculate x(t1) - x(t2)
        trajDiff = traj(j).observations(i+1:end) - traj(j).observations(i);
        
        %get number of pairs with x(t1) > x(t2)
        numIncPairs = numIncPairs + length(find(trajDiff>0));
        
        %get number of pairs tested
        numPairsTested = numPairsTested + length(find(~isnan(trajDiff)));
    
    end
end

%calculate mean and standard deviation of distribution
meanIncPairs = numPairsTested/2;
stdIncPairs = sqrt(sum(availLength)*(sum(availLength)-1)*...
    (2*sum(availLength)+5)/72);

%calculate the test statistic
testStatistic = (numIncPairs-meanIncPairs)/stdIncPairs;

%get the p-value of the test statistic assuming a standard normal distribution
pValue = 1 - normcdf(abs(testStatistic),0,1);

if pValue < significance/2 %if p-value is smaller than probability of type I error
    H = 1; %reject hypothesis that series is IID => there is a linear trend
else %if p-value is larger than probability of type I error
    H = 0; %cannot reject hypothesis
end


%%%%% ~~ the end ~~ %%%%%

