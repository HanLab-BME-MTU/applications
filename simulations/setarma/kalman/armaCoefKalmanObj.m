function neg2LnLikelihood = armaCoefKalmanObj(param,arOrder,trajectories)
%ARMACOEFKALMANOBJ calculates -2ln(likelihood) of the fit of an ARMA model to time series which could have missing data points
%
%SYNOPSIS neg2LnLikelihood = armaCoefKalmanObj(param,arOrder,trajectories)
%
%INPUT  param       : Set of partial ARMA coefficients.
%       arOrder     : Order of autoregressive part of process.
%       trajectories: Structure array containing trajectories to be modeled:
%           .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%           .available: Indices of existing observations.
%
%OUTPUT neg2LnLikelihood: Value of -2ln(likelihood).
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper.
%
%
%Khuloud Jaqaman, July 2004

%initialize output
neg2LnLikelihood = [];

%check if correct number of arguments were used when function was called
if nargin ~= nargin('armaCoefKalmanObj')
    disp('--armaCoefKalmanObj: Incorrect number of input arguments!');
    return
end

%assign parameters
arParamP = param(1:arOrder);
maParamP = param(arOrder+1:end);

%get AR and MA coefficients
if ~isempty(arParamP)
    [arParam,errFlag] = levinsonDurbinAR(arParamP);
else
    arParam = [];
end
if ~isempty(maParamP)
    [maParam,errFlag] = levinsonDurbinMA(maParamP);
else
    maParam = [];
end

%go over all trajectories to get likelihood
sum1 = 0;
sum2 = 0;
numAvail = 0;
for i = 1:length(trajectories)
    
    %calculate the number of available points in this trajectory
    available = trajectories(i).available;
%     available = available(find(available>10));
    numAvail = numAvail + length(available);
    
    %get the innovations, their variances and process white noise
    %using Kalman prediction and filtering
    [innovation,innovationVar,wnVector,errFlag] = ...
        armaKalmanInnov(trajectories(i).observations,arParam,maParam);
    if errFlag
        return
    end
        
    %1st sum in Eq. 3.15
    sum1 = sum1 + sum(log(innovationVar(available)));
    %2nd sum in Eq. 3.15
    sum2 = sum2 + sum(innovation(available).^2./innovationVar(available));
    
end

%construct -2ln(likelihood)
neg2LnLikelihood = sum1 + numAvail*log(sum2);

% % %TOMLAB stuff
% % % arOrder = prob.user.arOrder;
% % % trajectories = prob.user.trajectories;

