function neg2LnLikelihood = armaCoefKalmanObj(param,arOrder,maOrder,traj,available)
%ARMACOEFKALMANOBJ calculates -2ln(likelihood) of the fit of an ARMA model to time series which could have missing data points
%
%SYNOPSIS neg2LnLikelihood = armaCoefKalmanObj(param,arOrder,maOrder,traj,available)
%
%INPUT  param    : Set of ARMA coefficients and observational error variance.
%       arOrder  : Order of autoregressive part of process.
%       maOrder  : Order of moving average part of process.
%       traj     : Observed trajectory (with measurement uncertainties).
%                  Missing points should be indicated with NaN.
%       available: Indices of available observations.
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
arParam = param(1:arOrder);
maParam = param(arOrder+1:end-1);
obsVariance = param(end);

%get the innovations and their variances using Kalman prediction and filtering
[innovation,innovationVar,errFlag] = armaKalmanInnov(traj,arOrder,maOrder,...
    arParam,maParam,obsVariance);

%construct -2ln(likelihood)
neg2LnLikelihood = sum(ln(innovationVar(available))) + length(available)*ln(sum(innovation(available).^2./innovationVar(available)));
