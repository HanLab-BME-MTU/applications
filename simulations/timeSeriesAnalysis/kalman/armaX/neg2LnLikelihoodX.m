function neg2LnLikelihoodV = neg2LnLikelihoodX(param,prob)
%NEG2LNLIKELIHOODX calculates -2ln(likelihood) of the fit of an ARMAX model to time series.
%
%SYNOPSIS neg2LnLikelihoodV = neg2LnLikelihoodX(param,prob)
%
%INPUT  param       : Set of parameters related to ARMAX coefficients.
%       prob        : structure in Tomlab format containing variables
%                     needed for calculations.
%          .user.arOrder   : Order of AR part of process.
%          .user.maOrder   : Order of MA part of process.
%          .user.trajOut   : Observations of time series to be fitted. 
%                            Either an array of structures 
%                            trajOut(1:nTraj).observations, where 
%                            "observations" is a 2D array of measurements
%                            and uncertainties, or a 2D array representing
%                            one single trajectory. 
%                            Missing points should be indicated with NaN.
%          .user.trajIn    : Observations of input time series. 
%                            Either an array of structures 
%                            trajIn(1:nTraj).observations, where 
%                            "observations" is a 2D array of measurements
%                            and uncertainties, or a 2D array representing
%                            one single trajectory. 
%                            Missing points should be indicated with NaN.
%          .user.numAvail  : Total number of available observations.
%          .user.constParam: Parameters to be constrained. Structure with
%                            2 fields:
%                         .ar: AR parameters.
%                         .ma: MA parameters.
%
%OUTPUT neg2LnLikelihoodV: Value of -2ln(likelihood).
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper.
%
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neg2LnLikelihoodV = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of input arguments are used
if nargin ~= nargin('neg2LnLikelihoodX')
    disp('--neg2LnLikelihoodX: Incorrect number of input arguments!');
    return
end

%get variables from structure "prob"
arOrder  = prob.user.arOrder;
maOrder  = prob.user.maOrder;
trajOut  = prob.user.trajOut;
trajIn   = prob.user.trajIn;
numAvail = prob.user.numAvail;

%check "trajOut" and turn it into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    clear tmp
elseif ~isfield(trajOut,'observations')
    disp('--neg2LnLikelihoodX: Please input trajOut in fields ''observations''!')
    errFlag = 1;
end

%check "trajIn" and turn into struct if necessary
if ~isstruct(trajIn)
    tmp = trajIn;
    clear trajIn
    trajIn.observations = tmp;
    clear tmp
elseif ~isfield(trajIn,'observations')
    disp('--neg2LnLikelihoodX: Please input trajIn in fields ''observations''!')
    errFlag = 1;
end

%assign parameters
arParamP = param(1:arOrder);
maParamP = param(arOrder+1:arOrder+maOrder);
xParam   = param(arOrder+maOrder+1:end);

%get AR and MA coefficients from the partial AR and MA coefficients, respectively
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Likelihood calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum1 = 0; %1st sum in Eq. 3.15
sum2 = 0; %2nd sum in Eq. 3.15

%go over all trajectories to get innovations and their variances
for i = 1:length(trajOut)

    %get the innovations, their variances and process white noise
    %using Kalman prediction and filtering
    [innovation,innovationVar,wnVector,errFlag] = ...
        armaXKalmanInnov(trajOut(i).observations,trajIn(i).observations,...
        arParam,maParam,xParam);

    %1st sum in Eq. 3.15
    sum1 = sum1 + nansum(log(innovationVar));
    %2nd sum in Eq. 3.15
    sum2 = sum2 + nansum(innovation.^2./innovationVar);

end

%construct -2ln(likelihood)
neg2LnLikelihoodV = (sum1 + numAvail*log(sum2))/1000;


%%%%% ~~ the end ~~ %%%%%
