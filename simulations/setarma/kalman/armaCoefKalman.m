function [arParam,maParam,wnVariance,residuals,aic,aicc,errFlag] = ...
    armaCoefKalman(traj,arParam0,maParam0)
%ARMACOEFKALMAN fits an ARMA(p,q) model to a time series which could have missing data points using Kalman prediction and filtering.
%
%SYNOPSIS [arParam,maParam,wnVariance,residuals,aic,aicc,errFlag] = ...
%    armaCoefKalman(traj,arParam0,maParam0)
%
%INPUT  traj        : Trajectory to be modeled (with measurement uncertainties).
%                     Missing points should be indicated with NaN.
%       arParam0    : Initial guess of autoregressive coefficients (row vector).
%       maParam0    : Initial guess of moving average coefficients (row vector).
%
%OUTPUT arParam    : Estimated AR parameters.
%       maParam    : Estimated MA parameters.
%       wnVariance : Estimated variance of white noise.
%       residuals  : Estimated variance of observational error.
%       aic        : Akaike's Information Criterion.
%       aicc       : Akaike's Information Criterion with bias Correction.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper. The main
%        difference is that I do not estimate the observational error
%        variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent).
%
%Khuloud Jaqaman, July 2004

%initialize output
arParam = [];
maParam = [];
wnVariance = [];
residuals = [];
aic = [];
aicc = [];
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
    if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
        traj = [traj zeros(trajLength,1)];
    else
        disp('--armaCoefKalman: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if ~isempty(arParam0)
    [nRow,arOrder] = size(arParam0);
    if nRow ~= 1
        disp('--armaCoefKalman: "arParam0" should be a row vector!');
        errFlag = 1;
    end
    r = abs(roots([-arParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Initial AR part should be causal (i.e. all roots of AR polynomial should be greater than 0)!');
        errFlag = 1;
    end
end
if ~isempty(maParam0)
    [nRow,maOrder] = size(maParam0);
    if nRow ~= 1
        disp('--armaCoefKalman: "maParam0" should be a row vector!');
        errFlag = 1;
    end
    r = abs(roots([maParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Initial MA part should be invertible (i.e. all roots of MA polynomial should be greater than 0)!');
        errFlag = 1;
    end
end
if errFlag
    disp('--armaCoefKalman: Please fix input data!');
    return
end

%initial values of parameters to be estimated
param0 = [arParam0 maParam0];

%obtain number of available observations and their indices
available = find(~isnan(traj(:,1)));
numAvail = length(available);

%define optimization options.
options = optimset('Display','off','maxIter',1e6,'maxFunEvals',1e6');

%minimize -2ln(likelihood) to get ARMA coefficients and variance of observational error
[params,minFunc,exitFlag,output] = fmincon(@armaCoefKalmanObj,param0,[],[],[],[],...
    [-50*ones(1,arOrder-1) -0.99 -50*ones(1,maOrder-1) -0.99],...
    [50*ones(1,arOrder-1) 0.99 50*ones(1,maOrder-1) 0.99],...
    @armaCoefKalmanConst,options,arOrder,traj,available);

%assign parameters obtained through minimization
arParam = params(1:arOrder);
maParam = params(arOrder+1:end);

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
[innovation,innovationVar,errFlag] = armaKalmanInnov(traj,arParam,maParam);

%calculate white noise variance
wnVariance = mean(innovation(available).^2./innovationVar(available));

%calculate -2ln(likelihood) 
neg2LnLikelihood = sum(log(innovationVar(available))) ...
    + numAvail*log(sum(innovation(available).^2./innovationVar(available)));

%get number of parameters estimated: arOrder AR coefficients, maOrder MA
%coefficients and white noise variance
numParam = arOrder + maOrder + 1;

%evaluate Akaike's Information Criterion
aic = neg2LnLikelihood + 2*numParam;

%evaluate the bias-corrected Akaike's Information Criterion
aicc = neg2LnLikelihood + 2*numParam*numAvail/(numAvail-numParam-1);
