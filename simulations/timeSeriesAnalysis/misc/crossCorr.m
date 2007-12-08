function [gamma,errFlag] = crossCorr(traj1,traj2,maxLag)
%CROSSCORR calculates the cross-correlation (at lag 0) between 2 series which might have missing observations
%
%SYNOPSIS [gamma,errFlag] = crossCorr(traj1,traj2,maxLag)
%
%INPUT  traj1, traj2: Observations of the 2 time series to be tested.
%                     Each is either an array of structures 
%                     traj(1:nTraj).observations, or a 2D array (value +
%                     std) representing a single trajectory. Missing 
%                     points should be indicated with NaN.
%       maxLag      : Maximu lag to calculate cross-correlation.
%
%OUTPUT gamma       : 2*maxLag+1-by-2 array of the normalized
%                     cross-correlation and its std.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS Input trajectories could have a nonzero mean. The algorithm accounts
%        for that.
%
%Khuloud Jaqaman, August 2007

%% Output

errFlag = 0;
gamma = [];

%% Input

%check if correct number of arguments was used when function was called
if nargin < 3
    disp('--crossCorr: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if maxLag <= 0
    disp('--crossCorr: Variable "maxLag" should be a positive integer!');
    errFlag = 1;
end

%check first trajectory and turn into struct if necessary
if ~isstruct(traj1)
    tmp = traj1;
    clear traj1
    traj1.observations = tmp;
    clear tmp
elseif ~isfield(traj1,'observations')
    disp('--crossCorr: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%check second trajectory and turn into struct if necessary
if ~isstruct(traj2)
    tmp = traj2;
    clear traj2
    traj2.observations = tmp;
    clear tmp
elseif ~isfield(traj2,'observations')
    disp('--crossCorr: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%get number of trajectories, length of each trajectory, and make sure that
%trajectories are column vectors - first
numTraj1 = length(traj1);
trajLength1 = zeros(numTraj1,1);
for i=1:numTraj1
    [trajLength1(i),nCol] = size(traj1(i).observations); %length of each trajectory
    if nCol > 2
        disp('--crossCorr: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end

%get number of trajectories, length of each trajectory, and make sure that
%trajectories are column vectors - second
numTraj2 = length(traj2);
trajLength2 = zeros(numTraj2,1);
for i=1:numTraj2
    [trajLength2(i),nCol] = size(traj2(i).observations); %length of each trajectory
    if nCol > 2
        disp('--crossCorr: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end

%make sure that both trajectories contain the same number of data points
if max(abs([numTraj1; trajLength1] - [numTraj2; trajLength2])) ~= 0
    disp('--crossCorr: Both trajectories must have the same number of data points!');
    errFlag = 1;
end

%ad-hoc criterion to ensure that there are enough data points
if sum(trajLength1((trajLength1>maxLag))-maxLag) < 3*maxLag
    disp('--crossCorr: Trajectories not long enough! Increase trajectory length or reduce maxLag');
    errFlag = 1;
end

if errFlag
    disp('--crossCorr: Please fix input data!');
    return
end

%% Crosscorrelation calculation

%shift each trajectory to make its mean zero
for iTraj = 1 : numTraj1
    observations = traj1(iTraj).observations(:,1);
    trajMean = nanmean(observations);
    observations = observations - trajMean;
    traj1(iTraj).observations(:,1) = observations;
end
for iTraj = 1 : numTraj2
    observations = traj2(iTraj).observations(:,1);
    trajMean = nanmean(observations);
    observations = observations - trajMean;
    traj2(iTraj).observations(:,1) = observations;
end

%calculate the collective standard deviation of each trajectory
observations = vertcat(traj1.observations);
observations = observations(:,1);
traj1Std = sqrt(nanmean((observations).^2));
observations = vertcat(traj2.observations);
observations = observations(:,1);
traj2Std = sqrt(nanmean((observations).^2));

%calculate crosscorrelation function for lags -maxLag through maxLag
gamma = zeros(2*maxLag+1,2);

for lag = -maxLag : -1 %for all lags < 0

    vecMult = [];
    for j = 1 : numTraj1 %(=numTraj2) go over all trajectories
        
        if trajLength1(j) > abs(lag) %which are longer than lag

            %find pairs of points to be used in calculation
            vec1 = traj1(j).observations(abs(lag)+1:end,1);
            vec2 = traj2(j).observations(1:end-abs(lag),1);

            %multiply vectors
            vecMult = [vecMult; vec1.*vec2];

        end
        
    end

    %calculate autocorrelation function (omit pairs with missing
    %observations)
    vecMult = vecMult(~isnan(vecMult))/traj1Std/traj2Std;
    gamma(lag+maxLag+1,:) = [mean(vecMult) std(vecMult)/sqrt(length(vecMult))];

end

for lag = 0 : maxLag %for all lags >= 0

    vecMult = [];
    for j = 1 : numTraj1 %(=numTraj2) go over all trajectories
        
        if trajLength1(j) > lag %which are longer than lag

            %find pairs of points to be used in calculation
            vec1 = traj1(j).observations(1:end-lag,1);
            vec2 = traj2(j).observations(lag+1:end,1);

            %multiply vectors
            vecMult = [vecMult; vec1.*vec2];

        end
        
    end

    %calculate autocorrelation function (omit pairs with missing
    %observations)
    vecMult = vecMult(~isnan(vecMult))/traj1Std/traj2Std;
    gamma(lag+maxLag+1,:) = [mean(vecMult) std(vecMult)/sqrt(length(vecMult))];

end

%% ~~ the end ~~
