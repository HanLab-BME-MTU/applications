function [gamma,errFlag] = autoCorr(traj,maxLag)
%AUTOCORR calculates the unbiased autocorrelation function of a time series with missing observations
%
%SYNOPSIS [gamma,errFlag] = autoCorr(traj,maxLag)
%
%INPUT  traj   : Observations of time series whose autocorrelation function
%                is to be calculated. An array of structures with the field
%                "observations". Missing points should be indicated with NaN.
%       maxLag : Maximum lag at which autocorrelation function is calculated.
%
%OUTPUT gamma  : Unbiased autocorrelation function of series, 
%                where gamma(i) is the autocorrelation at lag i-1.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

%initialize output
errFlag = 0;
gamma = [];

%check if correct number of arguments were used when function was called
if nargin ~= nargin('autoCorr')
    disp('--autoCorr: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if maxLag <= 0
    disp('--autoCorr: Variable "maxLag" should be a positive integer!');
    errFlag = 1;
end
numTraj = length(traj); %number of trajectories
for i=1:numTraj
    [trajLength(i),nCol] = size(traj(i).observations); %length of each trajectory
    if nCol > 1
        disp('--autoCorr: Each trajectory should be a column vector!');
        errFlag = 1;
    end
end
if sum(trajLength(find(trajLength>maxLag))-maxLag) < 3*maxLag
    disp('--autoCorr: Trajectories not long enough! Increase trajectories or reduce maxLag');
    errFlag = 1;
end
if errFlag
    disp('--autoCorr: Please fix input data!');
    return
end

%get mean of each trajectory
for i=1:numTraj
    trajMean(i) = mean(traj(i).observations(find(~isnan(traj(i).observations))));
end

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
    vecMult = vecMult(find(~isnan(vecMult)));

    %calculate autocorrelation function
    gamma(lag+1) = mean(vecMult);

end

%normalize autocorrelation function
gamma = gamma/gamma(1);
