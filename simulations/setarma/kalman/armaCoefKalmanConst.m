function [c,ceq] = armaCoefKalmanConst(param,arOrder,traj,available)
%ARMACOEFKALMANCONST imposes causality and invertibility constraints on ARMA coefficients
%
%SYNOPSIS [c,ceq] = armaCoefKalmanConst(param,arOrder,traj,available)
%
%INPUT  param    : Set of parameters in ARMA model (concat. of arParam and maParam)
%       arOrder  : Order of autoregressive part of process.
%       traj     : Observed trajectory.
%       available: Indices of available observations.
%
%OUTPUT c   : Constraint that roots should be larger than one.
%       ceq : Empty, but must be included due to matlab requirements.
%
%Khuloud Jaqaman, July 2004

%initialize output
c = [];
ceq = [];
errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('armaCoefKalmanConst')
    disp('--armaCoefKalmanConst: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end
    
%distribute parameters
arParam = param(1:arOrder);
maParam = param(arOrder+1:end);

%roots of AR polynomial
rAr = abs(roots([-arParam(end:-1:1) 1]))';

%roots of MA polynomial
rMa = abs(roots([maParam(end:-1:1) 1]))';

%nonlinear constraint (inequality)
c = 1.00001 - [rAr rMa];
