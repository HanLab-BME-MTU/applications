function [varCovMat,arParam,maParam,errFlag] = armaLeastSquares(...
    trajectories,wnVector,arOrder,maOrder,wnVariance)
%ARMALEASTSQUARES estimates the ARMA coefficients and their variance-covariance matrix of a trajectory whose residuals are known.
%
%SYNOPSIS [varCovMat,arParam,maParam,errFlag] = armaLeastSquares(...
%    trajectories,wnVector,arOrder,maOrder,wnVariance)
%
%INPUT  trajectories: Observations of time series to be fitted. Either an 
%                     array of structures traj(1:nTraj).observations, or a
%                     2D array representing one single trajectory. 
%           .observations: 2D array of measurements and their uncertainties.
%                     Missing points should be indicated with NaN.
%       wnVector    : Structure containing 1 field:
%           .observations: Estimated white noise series in a trajectory.
%       arOrder     : Order of AR part of process.
%       maOrder     : Order of MA part of process.
%       wnVariance  : Estimated variance of white noise in process.
%                     Optional. Default: Zero.
%
%OUTPUT varCovMat : Variance-Covariance matrix of estimated coefficients:
%           .cofactorMat  : Cofactor matrix.
%           .posterioriVar: A posteriori estimate of residuals' variance.
%       arParam   : Estimated AR coefficients.
%       maParam   : Estimated MA coefficients.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varCovMat = [];
arParam = [];
maParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--armaLeastSquares: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%check trajectory and turn it into struct if necessary
if ~isstruct(trajectories)
    tmp = trajectories;
    clear trajectories
    trajectories.observations = tmp;
    clear tmp
elseif ~isfield(trajectories,'observations')
    disp('--armaLeastSquares: Please input the trajectories in fields "observations"')
    errFlag = 1;
    return
end

if nargin < 5 || isempty(wnVariance)
    wnVariance = 0;
else
    if wnVariance < 0
        disp('--armaLeastSquares: White Noise Variance should be nonnegative!');
        errFlag = 1;
    end
end

for i=1:length(trajectories);
    [trajLength,nCol] = size(trajectories(i).observations);
    switch nCol
        case 1 %if no error is supplied
            if wnVariance == 0
                trajectories(i).observations = [trajectories(i).observations ...
                    ones(trajLength,1)]; %assume that there is no observational error
            else
                trajectories(i).observations = [trajectories(i).observations ...
                    sqrt(wnVariance)*ones(trajLength,1)]; %assume that there is no observational error
            end
        case 2 %if there is observational error
            trajectories(i).observations(:,2) = ...
                sqrt(trajectories(i).observations(:,2).^2+wnVariance);
        otherwise
            disp('--armaLeastSquares: "trajectories.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaLeastSquares: Please fix input data!');
    return
end

%get maximum order and sum of orders
maxOrder = max(arOrder,maOrder);
sumOrder = arOrder + maOrder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fitLength = 0;
prevPoints = [];
observations = [];
epsilon = [];

%go over all trajectories
for i = 1:length(trajectories)
    
    %obtain available points in this trajectory
    available = find(~isnan(trajectories(i).observations(:,1)));

    %find data points to be used in fitting
    fitSet = available(find((available(maxOrder+1:end)-available(1:end-maxOrder))...
        ==maxOrder) + maxOrder);
    fitLength1 = length(fitSet);
    fitLength = fitLength + fitLength1;
    
    %construct weighted matrix of previous points and errors multiplying 
    %AR and MA coefficients, respectively
    %[size: fitLength1 by sumOrder]
    prevPoints1 = zeros(fitLength1,sumOrder);
    for j = 1:arOrder
        prevPoints1(:,j) = trajectories(i).observations(fitSet-j,1)...
            ./trajectories(i).observations(fitSet,2);
    end    
    for j = arOrder+1:sumOrder
        prevPoints1(:,j) = wnVector(i).observations(fitSet-j+arOrder)...
            ./trajectories(i).observations(fitSet,2);
    end
    %put points from all trajectories together in 1 matrix
    prevPoints = [prevPoints; prevPoints1];
    
    %form vector of weighted observations
    observations = [observations; trajectories(i).observations(fitSet,1)...
            ./trajectories(i).observations(fitSet,2)];
    
    %form vector of weighted errors
    epsilon = [epsilon; wnVector(i).observations(fitSet)...
            ./trajectories(i).observations(fitSet,2)];
    
end %(for i = 1:length(trajectories))

%estimate ARMA coefficients
maParam = (prevPoints\observations)';
arParam = maParam(1:arOrder);
maParam = maParam(arOrder+1:sumOrder);

%calculate variance-covariance matrix
varCovMat.cofactorMat = inv(prevPoints'*prevPoints);
varCovMat.posterioriVar = epsilon'*epsilon/(fitLength-sumOrder);


%%%%% ~~ the end ~~ %%%%%
