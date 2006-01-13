function constParamV = minKalmanConstraint(param,prob)
%MINKALMANCONSTRAINT determines the values of the parameters to be constrained
%
%SYNOPSIS constParamV = minKalmanConstraint(param,prob)
%
%INPUT  param       : Set of partial ARMA coefficients.
%       prob        : structure in Tomlab format containing variables
%                     needed for calculations.
%          .user.arOrder     : Order of autoregressive part of process.
%          .user.trajectories: Observations of time series to be fitted. 
%                              Either an array of structures 
%                              traj(1:nTraj).observations, where 
%                              "observations" is a 2D array of measurements
%                              and uncertainties, or a 2D array representing
%                              one single trajectory. 
%                              Missing points should be indicated with NaN.
%          .user.numAvail    : Total number of available observations.
%          .user.constParam  : Parameters to be constrained. Structure with
%                              2 fields:
%                           .ar : AR parameters.
%                           .ma : MA parameters.
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

%check if correct number of input arguments are used
if nargin ~= nargin('minKalmanConstraint')
    disp('--minKalmanConstraint: Incorrect number of input arguments!');
    return
end

%get variables from structure "prob"
arOrder = prob.user.arOrder;
constArIndx = prob.user.constParam.ar;
constMaIndx = prob.user.constParam.ma;

%assign parameters
arParamP = param(1:arOrder);
maParamP = param(arOrder+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%get current values of constrained parameters
constParamV = [arParam(constArIndx) maParam(constMaIndx)]';


%%%%% ~~ the end ~~ %%%%%
