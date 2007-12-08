function [gamma,errFlag] = autoCorr(traj,maxLag,correct)
%AUTOCORR calculates the unbiased autocorrelation function of a time series with missing observations
%
%SYNOPSIS [gamma,errFlag] = autoCorr(traj,maxLag,correct)
%
%INPUT  traj   : Observations of time series to be tested. Either an
%                array of structures traj(1:nTraj).observations, or a 2D
%                array (value + std) representing a single trajectory.
%                Missing points should be indicated with NaN.
%       maxLag : Maximum lag at which autocorrelation function is calculated.
%       correct: (optional) Switch for correction of trends:
%                  -1 :if no correction.
%                   0 :if correcting via linear fit to data.
%                   n :if correcting by taking the 1st difference at lag n
%                     (reduces traj length!).
%                   Default: -1.
%
%OUTPUT gamma  : Unbiased autocorrelation function of series,
%                where gamma(i) is the autocorrelation at lag i-1.
%                2nd column: std of gamma.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%REMARKS Input trajectories could have a nonzero mean. The algorithm accounts
%        for that.
%
%Khuloud Jaqaman, April 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = 0;
gamma = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin < 2
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

for i=1:numTraj

    %correct trajectories if necessary
    switch correct
        case -1
            %don't do anything
        case 0
            %remove trend with linear fit
            traj(i).observations = removeLinearTrend(traj(i).observations(:,1));
        otherwise %remove trend with differencing at specified lag
            tmp1 = traj(i).observations(correct+1:end,1) - ...
                traj(i).observations(1:end-correct,1);
            tmp2 = sqrt(traj(i).observations(correct+1:end,2).^2 + ...
                traj(i).observations(1:end-correct,2).^2);
            traj(i).observations = [tmp1 tmp2];
            %modify trajectory length
            trajLength(i) = trajLength(i) - correct;
    end

    %get trajectory mean
    trajMean(i) = nanmean(traj(i).observations(:,1));

end

% %shift trajectories so that the overall compound trajectory has a mean = 0
% meanAll = nanmean(vertcat(traj.observations));
% meanAll = meanAll(1);
% for i=1:numTraj
%     trajMean(i) = meanAll;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Autocorrelation function calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate unnormalized autocorrelation function for lags 0 through maxLag
gamma = zeros(maxLag+1,2);
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
    vecMult = vecMult(~isnan(vecMult));
    gamma(lag+1,:) = [mean(vecMult) std(vecMult)/sqrt(length(vecMult))];
    %     gamma(lag+1) = nansum(vecMult)/(length(vecMult)-1);

end

%normalize autocorrelation function
gamma = gamma/gamma(1);


%%%%% ~~ the end ~~ %%%%%
