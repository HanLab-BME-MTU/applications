function [gamma,errFlag] = autoCorr(traj,maxLag,correct)
%AUTOCORR calculates the unbiased autocorrelation function of a time series with missing observations
%
%SYNOPSIS [gamma,errFlag] = autoCorr(traj,maxLag)
%
%INPUT  traj   : Observations of time series whose autocorrelation function
%                is to be calculated. 
%                Either an array of structures traj(1:nTraj).observations, or 
%                one single trajectory. 
%                Missing points should be indicated with NaN.
%       maxLag : Maximum lag at which autocorrelation function is calculated.
%       correct: (optional) Switch for correction of trends.
%                  {-1} if no correction
%                   0   if correct via linear fit to data
%                   n   if correct via nth difference (reduces traj length!)
%
%OUTPUT gamma  : Unbiased autocorrelation function of series, 
%                where gamma(i) is the autocorrelation at lag i-1.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%REMARKS Input trajectories could have a nonzero mean. The algorithm accounts 
%        for that.
%
%Khuloud Jaqaman, April 2004

%initialize output and set defaults
errFlag = 0;
gamma = [];
defCorrect = -1;

%=============
% TEST INPUT
%=============

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

% check trajectory and turn into struct if necessary
if ~isstruct(traj)
    tmp = traj(:);
    clear traj
    traj.observations = tmp;
    clear tmp
elseif ~isfield(traj,'observations')
    disp('--autoCorr: Please input the trajectories in fields ''observations''')
    errFlag = 1;
end

% check input correct
if nargin > 2 && ~isempty(correct)
    % use user input
else
    correct = defCorrect;
end

numTraj = length(traj); %number of trajectories
trajMean = zeros(numTraj,1);
trajLength = zeros(numTraj,1);

for i=1:numTraj    
    
    % correct the trajectories if necessary and remember the mean
    switch correct
        case -1
            % remember mean
            trajMean(i) = nanmean(traj(i).observations);
        case 0
            % remove trend with linear fit
            traj(i).observations = removeLinearTrend(traj(i).observations);
            trajMean(i) = 0;
        otherwise
            traj(i).observations = ...
                traj(i).observations(correct+1:end) - traj(i).observations(1:end-correct);
            trajMean(i) = nanmean(traj(i).observations);
    end
    
    [trajLength(i),nCol] = size(traj(i).observations); %length of each trajectory
    if nCol > 1
        disp('--autoCorr: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end

% ad-hoc criterion to ensure that we have enough data points
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

%==============



%==============
% CALCULATE
%==============

%calculate unnormalized autocorrelation function for lags 0 through maxLag
gamma = zeros(maxLag+1,1);
for lag = 0:maxLag %for each lag
    
    vecMult = [];
    for j = 1:numTraj %go over all trajectories
        if trajLength(j) > lag %which are longer than lag

            %find pairs of points to be used in calculation
            vec1 = traj(j).observations(1:end-lag);
            vec2 = traj(j).observations(lag+1:end);

            %multiply mean-subtracted vectors
            vecMult = [vecMult; (vec1-trajMean(j)).*(vec2-trajMean(j))];

        end
    end

    %remove all pairs with missing observations
    vecMult = vecMult(~isnan(vecMult));

    %calculate autocorrelation function
    gamma(lag+1) = mean(vecMult);

end

%normalize autocorrelation function
gamma = gamma/gamma(1);
