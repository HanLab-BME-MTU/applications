function [H,errFlag] = rankTest(traj,significance)
%RANKTEST tests the hypothesis that a time series is IID by checking whether it has any linear trend.

%SYNOPSIS [H,errFlag] = rankTest(traj,significance)
%
%INPUT  traj        : Trajectory to be tested. Missing observations should be
%                     indicated with NaN.
%       significance: Significance level of hypothesis test. Default: 0.95.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.37. It can be used to detect linear
%       trends in a time series. 

%Khuloud Jaqaman, August 2004

%initialize output
H = [];
errFlag = 0;

%check input data
if nargin < 1
    disp('--rankTest: You should at least input time series to be analyzed!');
    errFlag = 1;
    return
end
availLength = length(find(~isnan(traj)));
trajLength = length(traj);

%assign default value of significance, if needed
if nargin < 2
    significance = 0.95;
end

%calculate number of points where x(t) > x(t-1)
numIncPairs = 0;
for i=1:trajLength-1
    numIncPairs = numIncPairs + length(find((traj(i+1:end)-traj(i))>0));
end

%calculate mean and standard deviation of distribution
meanIncPairs = availLength*(availLength-1)/4;
stdIncPairs = sqrt(availLength*(availLength-1)*(2*availLength+5)/72);

%calculate the test statistic
testStatistic = (numIncPairs-meanIncPairs)/stdIncPairs;

%get the p-value of the test statistic assuming a standard normal distribution
pValue = 1 - normcdf(abs(testStatistic),0,1);

if pValue < (1-significance)/2 %if p-value is smaller than limit
    H = 1; %reject hypothesis that series is IID => there is a linear trend
else %if p-value is larger than limit
    H = 0; %cannot reject hypothesis
end
