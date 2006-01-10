function [newTraj,trend,errFlag] = expoSmooth(traj,alpha)
%EXPOSMOOTH detects trends in a series by calculating a 1-sided moving average with exponential weights
%
%SYNPOSIS [newTraj,trend,errFlag] = expoSmooth(traj,alpha)
%
%INPUT  traj  : Column vector of series to be de-trended.
%       alpha : Factor by which new time point contributes to moving average.
%               Must be between 0 and 1, exclusive. 
%               Optional. Default: 0.25.
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
    disp('--expoSmooth: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
[nRow,nCol] = size(traj);
if nRow ~= 1 && nCol ~= 1
    disp('--expoSmooth: Trajectory should be a (column or row) vector!');
    errFlag = 1;
else
    trajLength = length(traj);
end

if nargin < 2 %if "alpha" was not input
    alpha = 0.25; %assign default
else
    if alpha <= 0 || alpha >= 1 %check that it has accetable value
        disp('--expoSmooth: "alpha" should be in the interval (0,1)!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--expoSmooth: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trend calculation and subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create vector
trend = NaN*ones(trajLength,1);

%find first time point that is available (i.e. not NaN)
iStart = 1;
while isnan(traj(iStart))
    iStart = iStart + 1;
end

%compute trend
trend(iStart) = traj(iStart);
for i=iStart+1:trajLength
    if ~isnan(traj(i))
        trend(i) = alpha*traj(i) + (1-alpha)*trend(i-1);
    else
        trend(i) = trend(i-1);
    end
end

%subtract trend from trajectory
if nCol == 1
    newTraj = traj - trend;
else
    newTraj = traj - trend';
end


%%%%% ~~ the end ~~ %%%%%

