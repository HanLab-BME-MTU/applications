function [c,ceq] = mlestimConst(param,arOrder,maOrder,traj)
%MLESTIMCONST imposes causality and invertibility constraints
%
%SYNOPSIS [c,ceq] = mlestimConst(param,traj,arOrder,maOrder)
%
%INPUT  param   : Set of parameters in ARMA model (concat. of arParam and maParam)
%       arOrder : Order of autoregressive part of process.
%       maOrder : Order of moving average part of process.
%       traj    : Observed trajectory.
%
%OUTPUT c   : Constraint that roots should be larger than one.
%       ceq : Empty, but must be included due to matlab requirements.
%
%Khuloud Jaqaman, February 2004

%check if correct number of arguments were used when function was called
if nargin ~= nargin('mlestimConst')
    disp('--mlestimConst: Incorrect number of input arguments!');
    errFlag  = 1;
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

%check input data
[nRow,nCol] = size(param);
if nRow ~= 1
    disp('--mlestimConst: "param" should be a row vector!');
    errFlag = 1;
end
if nCol ~= arOrder+maOrder
    disp('--mlestimConst: Wrong length of vector "param"!');
    errFlag = 1;
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

%nonlinear constraint equality
ceq = [];
