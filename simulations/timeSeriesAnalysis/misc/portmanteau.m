function [H,pValue,errFlag] = portmanteau(traj,maxLag,significance)
%PORTMANTEAU tests the hypothesis that a time series is IID by looking at its autocorrelation function.
%
%SYNOPSIS [H,pValue,errFlag] = portmanteau(traj,maxLag,significance)
%
%INPUT  traj        : Observations of time series to be tested. Either an 
%                     array of structures traj(1:nTraj).observations, or 
%                     a column representing one single trajectory. Missing 
%                     points should be indicated with NaN.
%       maxLag      : Maximum lag used to calculate autocorrelation
%                     function. Default: 40.
%       significance: Significance level of hypothesis test. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       pValue  : Probability of obtaining a test statistic >= observed
%                 value assuming that null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.36. In the original test, 1 time series 
%       which does not have any missing observations is used. Here the test
%       can be applied to several time series that are considered to be
%       realizations of the same process and which could have missing data
%       points. Thus I substitute the total number of points used in the test
%       for trajectory length in the original test.

%Khuloud Jaqaman, August 2004

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
    disp('--portmanteau: You should at least input time series to be analyzed!');
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
    disp('--portmaneau: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%check that trajectories are column vectors and get the total number of 
%available points
numAvail = 0; %initialize number of points to zero
for i=1:length(traj)
    nCol = size(traj(i).observations,2);
    if nCol > 1 %check that each trajectory is a colum vector
        disp('--portmanteau: Each trajectory should be a column vector!');
        errFlag = 1;
    else
        numAvail = numAvail + length(find(~isnan(traj(i).observations)));
    end
end

%assign default value of maxLag, if needed
if nargin < 2
    maxLag = 40;
end

%assign default value of significance, if needed
if nargin < 3
    significance = 0.05;
end

%exit if there are problems with input data
if errFlag
    disp('--portmanteau: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get autocorrelation function of time series
[gamma,errFlag] = autoCorr(traj,maxLag);
if errFlag
    disp('--portmanteau: "autoCorr" did not function properly!');
    return
end
gamma = gamma(2:end,1); %get rid of autocorrelation at 0 lag

%calculate the Ljung and Box test statistic
QLB = numAvail*(numAvail+2)*sum(gamma.^2./(numAvail-1:-1:numAvail-maxLag)');

%get the p-value of the test statistic assuming a chi2 distribution
pValue = 1 - chi2cdf(QLB,maxLag);

if pValue < significance %if p-value is smaller than probability of type I error
    H = 1; %reject null hypothesis that series is IID
else %if p-value is larger than probability of type I error
    H = 0; %cannot reject null hypothesis
end


%%%%% ~~ the end ~~ %%%%%

