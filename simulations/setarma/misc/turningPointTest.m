function [H,errFlag] = turningPointTest(traj,significance)
%TURNINGPOINTTEST tests the hypothesis that a time series is IID by looking at the number of turning points in it.

%SYNOPSIS [H,errFlag] = turningPointTest(traj,significance)
%
%INPUT  traj        : Trajectory to be tested. Missing observations should be
%                     indicated with NaN.
%       significance: Significance level of hypothesis test. Default: 0.95.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.36.

%Khuloud Jaqaman, August 2004

%initialize output
H = [];
errFlag = 0;

%check input data
if nargin < 1
    disp('--turningPointTest: You should at least input time series to be analyzed!');
    errFlag = 1;
    return
end
trajLength = length(find(~isnan(traj)));

%assign default value of significance, if needed
if nargin < 2
    significance = 0.95;
end

%check which points are turning points
turnPointTest = (traj(2:end-1)-traj(1:end-2)).*(traj(2:end-1)-traj(3:end));

%get number of turning points
numTurnPoints = length(find(turnPointTest>0));

%calculate mean and standard deviation of distribution of number of turning points
meanTurnPoints = 2*(trajLength-2)/3;
stdTurnPoints = sqrt((16*trajLength-29)/90);

%calculate the test statistic
testStatistic = (numTurnPoints-meanTurnPoints)/stdTurnPoints;

%get the p-value of the test statistic assuming a standard normal distribution
pValue = 1 - normcdf(abs(testStatistic),0,1);

if pValue < (1-significance)/2 %if p-value is smaller than limit
    H = 1; %reject hypothesis that series is IID
else %if p-value is larger than limit
    H = 0; %cannot reject hypothesis
end
