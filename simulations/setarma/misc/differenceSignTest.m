function [H,errFlag] = differenceSignTest(traj,significance)
%DIFFERENCESIGNTEST tests the hypothesis that a time series is IID by checking whether it has any linear trend.

%SYNOPSIS [H,errFlag] = differenceSignTest(traj,significance)
%
%INPUT  traj        : Trajectory to be tested. An array of structures
%                     with the field "observations". Missing data points
%                     in a trajectory should be indicated with NaN.
%       significance: Significance level of hypothesis test. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK This test is taken from Brockwell and Davis, "Introduction to Time
%       Series and Forecasting", p.37. It can be used to detect linear
%       trends in a time series. However, it cannot detect any cyclic
%       element. 
%       The original test is applied to a single time series which does not 
%       have any missing observations.
%       Here the test can be applied to several time series that are 
%       considered to be realizations of the same process and which could 
%       have missing data points. Thus I generalize the mean to be 1/2 all
%       points tested. As for the variance, in the book it is given by
%       (n+1)/12 = (n-1)/12+1/6 = (all points tested)/12+1/6. For
%       large n, 1/6 can be ignored. Thus I approximate the variance by
%       (all points tested)/12.
%
%Khuloud Jaqaman, August 2004

%initialize output
H = [];
errFlag = 0;

%check input data
if nargin < 1
    disp('--differenceSignTest: You should at least input time series to be analyzed!');
    errFlag = 1;
    return
end

%check that trajectories are column vectors
numTraj = length(traj);
for i=1:numTraj
    nCol = size(traj(i).observations,2);
    if nCol > 1
        disp('--turningPointTest: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end

%assign default value of significance, if needed
if nargin < 2
    significance = 0.05;
end

%calculate number of points where x(t) > x(t-1) and total number of
%points tested
numIncPoints = 0;
numPointsTested = 0;

for i = 1:numTraj

    %calculate x(t) - x(t-1)
    trajDiff = traj(i).observations(2:end) - traj(i).observations(1:end-1);
    
    %get number of points where x(t) > x(t-1) 
    numIncPoints = numIncPoints + length(find(trajDiff>0));

    %obtain number of points tested
    numPointsTested = numPointsTested + length(find(~isnan(trajDiff)));

end

%calculate mean and standard deviation of distribution
meanIncPoints = numPointsTested/2;
stdIncPoints = sqrt(numPointsTested/12);

%calculate the test statistic
testStatistic = (numIncPoints-meanIncPoints)/stdIncPoints;

%get the p-value of the test statistic assuming a standard normal distribution
pValue = 1 - normcdf(abs(testStatistic),0,1);

if pValue < significance/2 %if p-value is smaller than probability of type I error
    H = 1; %reject hypothesis that series is IID => there is linear trend
else %if p-value is larger than probability of type I error
    H = 0; %cannot reject hypothesis
end
