function [H,pValue,errFlag] = turningPointTest(traj,significance)
%TURNINGPOINTTEST tests the hypothesis that a time series is IID by looking at the number of turning points in it.
%
%SYNOPSIS [H,pValue,errFlag] = turningPointTest(traj,significance)
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
%       Series and Forecasting", p.36. 
%       The original test is applied to a single time series which does not 
%       have any missing observations.
%       Here the test can be applied to several time series that are 
%       considered to be realizations of the same process and which could 
%       have missing data points. Thus I generalize the mean to be 2/3 all
%       points tested. As for the variance, in the book it is given by
%       (16n-29)/90 = 8(n-2)/45+1/30 = 8(all points tested)/45+1/30. For
%       large n, 1/30 can be ignored. Thus I approximate the variance by
%       8(all points tested)/45.

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
    disp('--turningPointTest: You should at least input time series to be analyzed!');
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
    disp('--turningPointTest: Please input the trajectories in fields ''observations''')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of turning points and total number of points tested
numTurnPoints = 0;
numPointsTested = 0;

for i=1:numTraj %for each trajectory

    %check which points are turning points
    turnPointTest = (traj(i).observations(2:end-1)-traj(i).observations(1:end-2))...
        .*(traj(i).observations(2:end-1)-traj(i).observations(3:end));

    %calculate number of turning points
    numTurnPoints = numTurnPoints + length(find(turnPointTest>0));
    
    %obtain number of points tested
    numPointsTested = numPointsTested + length(find(~isnan(turnPointTest)));

end

%calculate mean and standard deviation of distribution of number of turning points
meanTurnPoints = 2*numPointsTested/3;
stdTurnPoints = sqrt(8*numPointsTested/45);

%calculate the test statistic
testStatistic = (numTurnPoints-meanTurnPoints)/stdTurnPoints;

%get the p-value of the test statistic assuming a standard normal distribution
pValue = 1 - normcdf(abs(testStatistic),0,1);

if pValue < significance/2 %if p-value is smaller than probability of type I error
    H = 1; %reject null hypothesis that series is IID
else %if p-value is larger than probability of type I error
    H = 0; %cannot reject null hypothesis
end


%%%%% ~~ the end ~~ %%%%%

