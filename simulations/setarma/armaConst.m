function [c,ceq] = armaConst(param,arOrder,maOrder,traj)
%ARMACONST imposes causality and invertibility constraints
%
%SYNOPSIS [c,ceq] = armaConst(param,traj,arOrder,maOrder)
%
%INPUT  param   : Set of parameters in ARMA model (concat. of arParam and maParam)
%       arOrder : Order of autoregressive part of process.
%       maOrder : Order of moving average part of process.
%       traj    : Observed trajectory.innovPredict
%
%OUTPUT c   : Constraint that roots should be larger than one.
%       ceq : Empty, but must be included due to matlab requirements.
%
%Khuloud Jaqaman, February 2004

%check if correct number of arguments were used when function was called
if nargin ~= nargin('armaConst')
    disp('--armaConst: Incorrect number of input arguments!');
    errFlag  = 1;
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

%distribute parameters
arParam = param(1:arOrder);
maParam = param(arOrder+1:end);

%roots of AR polynomial
if arOrder ~= 0
    rAr = abs(roots([-arParam(end:-1:1) 1]));
else
    rAr = [];
end

%roots of MA polynomial
if maOrder ~= 0
    rMa = abs(roots([maParam(end:-1:1) 1]));
else
    rMa = [];
end

%nonlinear constraint (inequality)
c = 1.001 - [rAr rMa];

%nonlinear constrain equality
ceq = [];
