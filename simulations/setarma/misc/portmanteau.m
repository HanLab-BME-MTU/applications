function [H,errFlag] = portmanteau(traj,maxLag,significance)
%PORTMANTEAU tests the hypothesis that a time series is IID by looking at its autocorrelation function.

%SYNOPSIS [H,errFlag] = portmanteau(traj,maxLag,significance)
%
%INPUT  traj        : Trajectories to be tested. An array of structures
%                     with the field "observations". Missing data points
%                     in a trajectory should be indicated with NaN.
%       maxLag      : Maximum lag used to calculate autocorrelation
%                     function. Default: 40.
%       significance: Significance level of hypothesis test. Default: 0.95.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.36. Use with caution when there are
%       missing observations.

%Khuloud Jaqaman, August 2004

%initialize output
H = [];
errFlag = 0;

%check input data
if nargin < 1
    disp('--portmanteau: You should at least input time series to be analyzed!');
    errFlag = 1;
    return
end

%check that trajectories are column vectors and get their overall length
trajLength = 0;
for i=1:length(traj)
    nCol = size(traj(i).observations,2);
    if nCol > 1
        disp('--portmanteau: Each trajectory should be a column vector!');
        errFlag = 1;
    else
        trajLength = trajLength + length(find(~isnan(traj(i).observations)));
    end
end

%assign default value of maxLag, if needed
if nargin < 2
    maxLag = 40;
end

%assign default value of significance, if needed
if nargin < 3
    significance = 0.95;
end

%get autocorrelation function of time series
[gamma,errFlag] = autoCorr(traj,maxLag);
if errFlag
    disp('--portmanteau: "autoCorr" did not function properly!');
    return
end
gamma = gamma(2:end); %get rid of autocorrelation at 0 lag

%calculate the Ljung and Box test statistic
QLB = trajLength*(trajLength+2)*sum(gamma.^2./[trajLength-1:-1:trajLength-maxLag]');

%get the p-value of the test statistic assuming a chi2 distribution
pValue = 1 - chi2cdf(QLB,maxLag);

if pValue < 1-significance %if p-value is smaller than limit
    H = 1; %reject hypothesis that series is IID
else %if p-value is larger than limit
    H = 0; %cannot reject hypothesis
end

