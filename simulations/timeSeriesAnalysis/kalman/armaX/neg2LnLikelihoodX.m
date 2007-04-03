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
% % % % %          .user.xLag      : Lag in dependence on exogenous variable. (NOT USED CURRENTLY)
%          .user.trajOut   : Observations of time series to be fitted. 
%                            An array of structures:
%                    .observations: 2D array of measurements and
%                                   uncertainties. Missing points should be 
%                                   indicated with NaN.
%                    .weight      : The weight with which each movie
%                                   belongs to the group.
%          .user.trajIn    : Observations of input time series. 
%                            An array of structures 
%                            trajIn(1:nTraj).observations, where 
%                            "observations" is a 2D array of measurements
%                            and uncertainties. Must not have any missing
%                            points. Enter [] in "observations" field if 
%                            there is no input series.
%          .user.numAvail  : Total number of available observations.
%          .user.constParam: Parameters to be constrained. Structure with
%                            2 fields:
%                         .ar: AR parameters.
%                         .ma: MA parameters.
%                         .x : X parameters.
%
%OUTPUT neg2LnLikelihoodV: Value of -2ln(likelihood).
%
%REMARKS The algorithm implemented here is an ARMAX generalized version of
%        the algorithm presented in R. H. Jones,"Maximum Likelihood Fitting 
%        of ARMA Models to Time Series with Missing Observations", 
%        Technometrics 22: 389-395 (1980). All equation numbers used are 
%        those in the paper.
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
% xLag     = prob.user.xLag;
trajOut  = prob.user.trajOut;
trajIn   = prob.user.trajIn;
numAvail = prob.user.numAvail;

%assign parameters
arParamP = param(1:arOrder)';
maParamP = param(arOrder+1:arOrder+maOrder)';
xParam   = param(arOrder+maOrder+1:end)';
% if ~isempty(xLag) && xLag ~= 0
%     xParam = [zeros(xLag,1); xParam];
% end

%get AR and MA coefficients from the partial AR and MA coefficients, respectively
if ~isempty(arParamP)
    [arParam,errFlag] = levinsonDurbinExpoAR(arParamP);
else
    arParam = [];
end
if ~isempty(maParamP)
    [maParam,errFlag] = levinsonDurbinExpoMA(maParamP);
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
    [innovation,innovationVar,wnVector,dummy1,dummy2,errFlag] = ...
        armaxKalmanInnov(trajOut(i).observations,trajIn(i).observations,...
        arParam,maParam,xParam);

    %1st sum in Eq. 3.15
    sum1 = sum1 + trajOut(i).weight*nansum(log(innovationVar));
    %2nd sum in Eq. 3.15
    sum2 = sum2 + trajOut(i).weight*nansum(innovation.^2./innovationVar);

end

%construct -2ln(likelihood) (Eq. 3.15)
neg2LnLikelihoodV = (sum1 + numAvail*log(sum2));


%%%%% ~~ the end ~~ %%%%%
