function [newTraj,trend,errFlag] = movAvSmooth(traj,window)
%movAvSmooth detects trends in a series by calculating a 2-sided moving average with a cut-off
%
%SYNPOSIS [newTraj,trend,errFlag] = movAvSmooth(traj,window)
%
%INPUT  traj  : Column vector of series to be de-trended.
%       window: Number of points to be used on each side for averaging.
%               Optional. Default: 2.
%
%OUTPUT newTraj: Column vector of de-trended series.
%       trend  : Estimated trend.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, October 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newTraj = [];
trend = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin < 1 || nargin > 2
    disp('--movAvSmooth: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
[trajLength,nCol] = size(traj);
if nCol ~= 1
    disp('--movAvSmooth: Trajectory should be a column vector!');
    errFlag = 1;
end

if nargin < 2 %if "window" was not input
    window = 2; %assign default
else
    if window <= 0 %check that it has accetable value
        disp('--movAvSmooth: Number of points used for averaging on each side should be >= 1!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--movAvSmooth: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trend calculation and subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%append trajectory at both ends to allow trend calculation for points
%less than "window" points away from either end
traj = [traj(1)*ones(window,1); traj; traj(end)*ones(window,1)];

%compute moving average
trend = NaN*ones(trajLength,1);
for i=1:trajLength
    trend(i) = nanmean(traj(i:i+2*window));
end

%subtract trend from trajectory
traj = traj(window+1:end-window);
newTraj = traj - trend;


%%%%% ~~ the end ~~ %%%%%

