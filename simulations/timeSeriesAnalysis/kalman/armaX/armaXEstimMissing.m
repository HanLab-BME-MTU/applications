function [trajFilled,errFlag] = armaEstimMissing(traj,arParam,maParam,...
    wnVariance)
%ARMAESTIMMISSING estimated missing observations using Kalman smoothing
%
%SYNOPSIS [trajFilled,errFlag] = armaEstimMissing(traj,arParam,maParam,...
%    wnVariance)
%
%INPUT  traj      : Observations of time series to be fitted. Either an
%                   array of structures traj(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Missing points should be indicated with NaN.
%       arParam   : Autoregressive coefficients (row vector).
%       maParam   : Moving average coefficients (row vector).
%       wnVariance: White noise variance.
%
%OUTPUT trajFilled: Trajectories with estimated values of missing
%                   observations. Structure with fields observations.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%MATLAB VERSION (originally written on): 7.0.4.352 (R14) Service Pack 2 
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajFilled = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('armaEstimMissing')
    disp('--armaEstimMissing: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check trajectory and turn it into struct if necessary
if ~isstruct(traj)
    tmp = traj;
    clear traj
    traj.observations = tmp;
    clear tmp
elseif ~isfield(traj,'observations')
    disp('--armaEstimMissing: Please input the trajectories in fields ''observations''!')
    errFlag = 1;
    return
end

%get number of trajectories
numTraj = length(traj);

%get length of each trajectory, number of missing observations in it,
%and add column for observational error if not input
trajLength = zeros(numTraj,1);
numMissing = zeros(numTraj,1);
for i=1:numTraj
    trajT = traj(i).observations;
    [trajLength(i),nCol] = size(trajT);
    numMissing(i) = length(find(isnan(trajT(:,1))));
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            trajT = [trajT zeros(trajLength,1)];
        else
            disp('--armaEstimKalman: "traj.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    traj(i).observations = trajT;
end

%get arOrder
[nRow,arOrder] = size(arParam);
if ~isempty(arParam)
    if nRow ~= 1
        disp('--armaEstimMissing: "arParam" should be a row vector!');
        errFlag = 1;
    end
end

%get maOrder
[nRow,maOrder] = size(maParam);
if ~isempty(maParam)
    if nRow ~= 1
        disp('--armaEstimMissing: "maParam" should be a row vector!');
        errFlag = 1;
    end
end

%check white noise variance
if wnVariance < 0
    disp('--armaEstimMissing: "wnVariance" should be nonnegative!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--armaEstimMissing: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smoothing algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go over all supplied trajectories
for i=1:numTraj
    
    %get forward prediction of state vector and corresponding covariance matrix
    %at missing observations
    [dummy1,dummy2,dummy3,stateVecFor,stateCovFor,errFlag] = ...
        armaKalmanInnov(traj(i).observations,arParam,maParam);

    %get backward prediction of state vector and corresponding covariance
    %matrix at missing observations
    [stateVecBack,stateCovBack,errFlag] = armaKalmanInnov(traj(i).observations,...
        arParam,maParam);

    %calculate the best estimate of missing observations using all observations
    for j=1:numMissing(i)
        tmp = stateCovFor(:,:,j)*stateCovBack(:,:,j);
        sizeT = size(tmp,1);
        identityMat = eye(sizeT);
        gainSmooth = tmp*inv(identityMat + tmp);
        stateCovSmooth = (identityMat - gainSmooth)*stateCovFor(:,:,j);
        stateVecSmooth = (identityMat - gainSmooth)*stateVecFor(:,j) + ...
            stateCovSmooth*stateVecBack(:,j);
    end
    
end %(for i=1:numTraj)


%%%%% ~~ the end ~~ %%%%%
