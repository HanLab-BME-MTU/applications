function [stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,maContribution)
%COVKALMANINIT calculates the covariance matrix of the initial state used in Kalman recursions.
%
%SYNOPSIS [stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,maContribution)
%
%INPUT arParam       : Autoregressive coefficients (row vector).
%      maParam       : Moving average coefficients (row vector).
%      maContribution: Vector G (Eq. 2.15, 2.16) which is equivalent to the
%                      covariance of the process with the white noise
%                      (column vector).
%
%OUTPUT stateCovMat00: Covariance matrix of initial state.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2004

%initialize output
stateCovMat00 = [];
errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < nargin('armaCovMat00')
    disp('--armaCovMat00: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

