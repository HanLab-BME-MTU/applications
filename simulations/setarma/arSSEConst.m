function [c,ceq] = arSSEConst(unknown,arOrder,traj)
%arSSEConst imposes causality constraint on AR coefficients
%
%SYNOPSIS [c,ceq] = arSSEConst(unknown,arOrder,traj)
%
%INPUT  unknown : Set of AR coefficients and measurement error-free values. 
%       arOrder : Order of autoregressive model.
%       traj    : Observed trajectory.
%
%OUTPUT c   : Constraint that roots should be larger than one.
%       ceq : Empty, but must be included due to matlab requirements.
%
%Khuloud Jaqaman, February 2004

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arSSEConst')
    disp('--arSSEConst: Incorrect number of input arguments!');
    errFlag  = 1;
    c = [];
    ceq = [];
    return
end

%check input data
[nRow,nCol] = size(unknown);
if nRow ~= arOrder+length(traj(:,1))
    disp('--arSSEConst: Wrong length of vector "unknown"!');
    errFlag = 1;
end
if nCol ~= 1
    disp('--arSSEConst: "param" should be a column vector!');
    errFlag = 1;
end

%distribute parameters
arParam = unknown(1:arOrder);

%roots of AR polynomial
rAr = abs(roots([-arParam(end:-1:1) 1]))';

%nonlinear constraint (inequality)
c = 1.001 - rAr;

%nonlinear constraint equality
ceq = [];