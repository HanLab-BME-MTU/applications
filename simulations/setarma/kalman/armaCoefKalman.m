function [arParam,maParam,wnVariance,obsVariance,aic,aicc,errFlag] = ...
    armaCoefKalman(traj,arOrder,maOrder,arParam0,maParam0,obsVariance0)
%ARMACOEFKALMAN fits an ARMA(p,q) model to, and estimates the observation error variance of, a time series which could have missing data points using Kalman prediction and filtering.
%
%SYNOPSIS [arParam,maParam,wnVariance,obsVariance,aic,aicc,errFlag] = ...
%    armaCoefKalman(traj,arOrder,maOrder,arParam0,maParam0,obsVariance0)
%
%INPUT  traj        : Trajectory to be modeled (with measurement uncertainties).
%                     Missing points should be indicated with NaN.
%       arOrder     : Order of AR part of proposed ARMA model.
%       maOrder     : Order of MA part of proposed ARMA model.
%       arParam0    : Initial guess of autoregressive coefficients (row vector).
%       maParam0    : Initial guess of moving average coefficients (row vector).
%       obsVariance0: Initial guess of observational error variance.
%
%OUTPUT arParam     : Estimated AR parameters.
%       maParam     : Estimated MA parameters.
%       wnVariance  : Estimated variance of white noise.
%       obsVariance : Estimated variance of observational error.
%       aic         : Akaike's Information Criterion.
%       aicc        : Akaike's Information Criterion with bias Correction.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper.
%
%Khuloud Jaqaman, July 2004

%initialize output
arParam = [];
maParam = [];
wnVariance = [];
obsVariance = [];
errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < nargin('armaCoefKalman')
    disp('--armaCoefKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
[trajLength,nCol] = size(traj);
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--armaCoefKalman: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if arOrder < 0
    disp('--armaCoefKalman: AR order cannot be negative!');
    errFlag = 1;
end
if maOrder < 0
    disp('--armaCoefKalman: MA order cannot be negative!');
    errFlag = 1;
end
if arOrder == 0 && maOrder == 0
    disp('--armaCoefKalman: Either Ar order or MA order should be different from zero!');
    errFlag = 1;
end
if ~isempty(arParam0)
    [nRow,nCol] = size(arParam0);
    if nRow ~= 1
        disp('--armaCoefKalman: "arParam0" should be a row vector!');
        errFlag = 1;
    end
    if nCol ~= arOrder
        disp('--armaCoefKalman: "arParam0" should have "arOrder" entries!');
        errFlag = 1;
    end
    r = abs(roots([-arParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Initial AR part should be causal (i.e. all roots of AR polynomial should be greater than 0)!');
        errFlag = 1;
    end
end
if ~isempty(maParam0)
    [nRow,nCol] = size(maParam0);
    if nRow ~= 1
        disp('--armaCoefKalman: "maParam0" should be a row vector!');
        errFlag = 1;
    end
    if nCol ~= maOrder
        disp('--armaCoefKalman: "maParam0" should have "maOrder" entries!');
        errFlag = 1;
    end
    r = abs(roots([maParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Initial MA part should be invertible (i.e. all roots of MA polynomial should be greater than 0)!');
        errFlag = 1;
    end
end
if obsVariance0 < 0
    disp('--armaCoefKalman: Observational error variance cannot be negative!');
    errFlag = 1;
end
if errFlag
    disp('--armaCoefKalman: Please fix input data!');
    return
end

%initial values of parameters to be estimated
param0 = [arParam0 maParam0 obsVariance0];

%obtain number of available observations and their indices
available = find(~isnan(traj(:,1)));
numAvail = length(available);

%define optimization options.
options = optimset('Display','iter','maxIter',1e6,'maxFunEvals',1e6');

%minimize -2ln(likelihood) to get ARMA coefficients and variance of observational error
[params,minFunc,exitFlag,output] = fmincon(@armaCoefKalmanObj,param0,[],[],[],[],...
    [-inf*ones(1,arOrder+maOrder) 0],inf*ones(1,arOrder+maOrder+1),...
    @armaCoefKalmanConst,options,arOrder,maOrder,traj,available);

%assign parameters obtained through minimization
arParam = params(1:arOrder);
maParam = params(arOrder+1:end-1);
obsVariance = params(end);

%check for causality and invertibility of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--armaCoefKalman: Warning: Predicted model not causal!');
end
r = abs(roots([maParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--armaCoefKalman: Warning: Predicted model not invertible!');
end

%get the innovations and their variances using Kalman prediction and filtering
[innovation,innovationVar,errFlag] = armaKalmanInnov(traj,arOrder,maOrder,...
    arParam,maParam,obsVariance);

%calculate white noise variance
wnVariance = mean(innovation(available).^2./innovationVar(available));

%calculate -2ln(likelihood) 
neg2LnLikelihood = sum(log(innovationVar(available))) ...
    + numAvail*log(sum(innovation(available).^2./innovationVar(available)));

%get number of parameters estimated: arOrder AR coefficients, maOrder MA
%coefficients, white noise variance and observational error variance
numParam = arOrder + maOrder + 2;

%evaluate Akaike's Information Criterion
aic = neg2LnLikelihood + 2*numParam;

%evaluate the bias-corrected Akaike's Information Criterion
aicc = neg2LnLikelihood + 2*numParam*numAvail/(numAvail-numParam-1);
