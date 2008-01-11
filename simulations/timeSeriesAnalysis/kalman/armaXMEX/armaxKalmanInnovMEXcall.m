function [innovations,innovationVariances,whiteNoise,dummy1,dummy2,errFlag] = ...
    armaxKalmanInnovMEXcall(trajOut,trajIn,arParam,maParam,xParam,wnVariance)

%ARMAXKALMANINNOV does forward Kalman prediction and filtering of a time series using an ARMAX model
%
%This is simply a port to C which does the same thing faster.
%
%SYNOPSIS [innovation,innovationVar,wnVector,stateVec,stateCov,errFlag] = ...
%    armaxKalmanInnov(trajOut,trajIn,arParam,maParam,xParam,wnVariance)
%
%INPUT  trajOut   : Trajectory to be modeled (with measurement uncertainties).
%                   Missing points should be indicated with NaN.
%       trajIn    : Input series. Must not have any missing points.
%                   Enter [] if there is no input series.
%       arParam   : Autoregressive coefficients (row vector).
%       maParam   : Moving average coefficients (row vector).
%       xParam    : Coefficients indicating dependence on input (row vector).
%       wnVariance: White noise Variance. Optional. Default: 1.
%
%OUTPUT innovation   : Vector of differences between predicted and observed data, or innovations.
%       innovationVar: Vector of innovation variances.
%       wnVector     : Estimated white noise in the process.
%       stateVec     : Forward-predicted state vector at missing time points.
%       stateCov     : Forward-predicted covariance matrix at missing time points.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is an ARMAX generalized version of
%        the algorithm presented in R. H. Jones,"Maximum Likelihood Fitting 
%        of ARMA Models to Time Series with Missing Observations", 
%        Technometrics 22: 389-395 (1980). All equation numbers used are 
%        those in the paper. However, I do not estimate the observational 
%        error variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent).
%
%Khuloud Jaqaman, January 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

innovations = [];
innovationVariances = [];
whiteNoise = [];
dummy1 = [];
dummy2 = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < 5
    disp('--armaxKalmanInnov: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%find trajectory length and number of missing observations
trajLength = size(trajOut,1);

%make sure that input and output series have the same length
if ~isempty(trajIn)
    if size(trajIn,1)~=trajLength
        disp('--armaxKalmanInnov: Input and output series must have the same length!');
        errFlag = 1;
        return
    end
end

%assign 1 to WN variance if not supplied
if nargin < 6 || isempty(wnVariance)
    prob.wnVariance = 1;
else
    prob.wnVariance = wnVariance;
end

%find arOrder, maOrder and xOrder
prob.arOrder = length(arParam);
prob.maOrder = length(maParam);
prob.xOrder  = length(xParam) - 1;

paramV = cat(2,arParam,maParam,xParam)';

if (prob.xOrder < 1) || isempty(trajIn);
   trajIn = zeros(trajLength,2);   
end


prob.TRAJ = cat(3, trajIn,trajOut);
prob.TRAJ = permute( prob.TRAJ, [1 3 2]);

%Call the mex function which actually calculates the innovations

[innovations,innovationVariances,whiteNoise] = armaxKalmanInnovMEX(paramV,prob);
