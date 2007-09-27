function constParamV = minKalmanConstraint(param,prob)
%MINKALMANCONSTRAINT determines the values of the parameters to be constrained
%
%SYNOPSIS constParamV = minKalmanConstraint(param,prob)
%
%INPUT  param       : Set of ARMAX coefficients (ARMA coefficients as partial).
%       prob        : structure in Tomlab format containing variables
%                     needed for calculations.
%          .user.arOrder   : Order of AR part of process.
%          .user.maOrder   : Order of MA part of process.
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
%OUTPUT constParamV : Vector of Values of constrained parameters.
%
%Khuloud Jaqaman, January 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constParamV = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of input arguments is used
if nargin ~= nargin('minKalmanConstraint')
    disp('--minKalmanConstraint: Incorrect number of input arguments!');
    return
end

%get variables from structure "prob"
arOrder = prob.user.arOrder;
maOrder = prob.user.maOrder;
constArIndx = prob.user.constParam.ar;
constMaIndx = prob.user.constParam.ma;
constXIndx  = prob.user.constParam.x;

%assign parameters
arParamP = param(1:arOrder);
maParamP = param(arOrder+1:arOrder+maOrder);
xParam   = param(arOrder+maOrder+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%get current values of constrained parameters
constParamV = [arParam(constArIndx) maParam(constMaIndx) xParam(constXIndx+1)]';


%%%%% ~~ the end ~~ %%%%%
