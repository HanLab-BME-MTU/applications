function [corrVal,errFlag] = crossCorr(traj1,traj2)
%CROSSCORR calculates the cross-correlation (at lag 0) between 2 series which might have missing observations
%
%SYNOPSIS [corrVal,errFlag] = crossCorr(traj1,traj2)
%
%INPUT  traj1, traj2: Observations of the 2 time series to be tested.
%                     Each is either an array of structures 
%                     traj(1:nTraj).observations, or a 2D array (value +
%                     std) representing a single trajectory. Missing 
%                     points should be indicated with NaN.
%
%OUTPUT corrVal     : Normalized cross-correlation (between -1 and 1).
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS Input trajectories could have a nonzero mean. The algorithm accounts
%        for that.
%
%Khuloud Jaqaman, June 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = 0;
gamma = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin < 2 || nargin > 3
    disp('--autoCorr: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if maxLag <= 0
    disp('--autoCorr: Variable "maxLag" should be a positive integer!');
    errFlag = 1;
end

%check trajectory and turn into struct if necessary
if ~isstruct(traj)
    tmp = traj;
    clear traj
    traj.observations = tmp;
    clear tmp
elseif ~isfield(traj,'observations')
    disp('--autoCorr: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%check input "correct"
if nargin > 2 && ~isempty(correct)
    % use user input
else %use default
    correct = -1;
end

%get number of trajectories, length of each trajectory, and make sure that
%trajectories are column vectors
numTraj = length(traj);
trajLength = zeros(numTraj,1);
for i=1:numTraj
    [trajLength(i),nCol] = size(traj(i).observations); %length of each trajectory
    if nCol > 2
        disp('--autoCorr: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end

%ad-hoc criterion to ensure that there are enough data points
if sum(trajLength((trajLength>maxLag))-maxLag) < 3*maxLag
    disp('--autoCorr: Trajectories not long enough! Increase trajectories or reduce maxLag');
    errFlag = 1;
end

if errFlag
    if nargout == 2
        disp('--autoCorr: Please fix input data!');
        return
    else
        error('--autoCorr: Please fix input data!')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trend correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajMean = zeros(numTraj,1);

for i=1:numTraj

    %correct trajectories if necessary and remember the mean
    switch correct
        case -1
            %remember mean
            trajMean(i) = nanmean(traj(i).observations(:,1));
        case 0
            %remove trend with linear fit
            traj(i).observations = removeLinearTrend(traj(i).observations(:,1));
            trajMean(i) = 0;
        otherwise %remove trend with differencing at specified lag
            tmp1 = traj(i).observations(correct+1:end,1) - ...
                traj(i).observations(1:end-correct,1);
            tmp2 = sqrt(traj(i).observations(correct+1:end,2).^2 + ...
                traj(i).observations(1:end-correct,2).^2);
            traj(i).observations = [tmp1 tmp2];
            trajMean(i) = nanmean(traj(i).observations(:,1));
            %modify trajectory length
            trajLength(i) = trajLength(i) - correct;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Autocorrelation function calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate unnormalized autocorrelation function for lags 0 through maxLag
gamma = zeros(maxLag+1,1);
for lag = 0:maxLag %for each lag

    vecMult = [];
    for j = 1:numTraj %go over all trajectories
        if trajLength(j) > lag %which are longer than lag

            %find pairs of points to be used in calculation
            vec1 = traj(j).observations(1:end-lag,1);
            vec2 = traj(j).observations(lag+1:end,1);

            %multiply mean-subtracted vectors
            vecMult = [vecMult; (vec1-trajMean(j)).*(vec2-trajMean(j))];

        end
    end

    %calculate autocorrelation function (omit pairs with missing observations)
    gamma(lag+1) = nanmean(vecMult);
    %     gamma(lag+1) = nansum(vecMult)/(length(vecMult)-1);

end

%normalize autocorrelation function
gamma = gamma/gamma(1);


%%%%% ~~ the end ~~ %%%%%
