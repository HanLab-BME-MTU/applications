function [c,ceq] = armaCoefKalmanConst(param,arOrder,maOrder,traj,available)
%ARMACOEFKALMANCONST imposes causality and invertibility constraints on ARMA coefficients
%
%SYNOPSIS [c,ceq] = armaCoefKalmanConst(param,arOrder,maOrder,traj,available)
%
%INPUT  param    : Set of parameters in ARMA model (concat. of arParam and maParam)
%       arOrder  : Order of autoregressive part of process.
%       maOrder  : Order of moving average part of process.
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

%check input data
[nRow,nCol] = size(param);
if nRow ~= 1
    disp('--armaCoefKalmanConst: "param" should be a row vector!');
    errFlag = 1;
end
if nCol ~= arOrder+maOrder+1
    disp('--armaCoefKalmanConst: Wrong length of vector "param"!');
    errFlag = 1;
end
if errFlag
    disp('--armaCoefKalmanConst: Please fix input data!');
    return
end
    
%distribute parameters
arParam = param(1:arOrder);
maParam = param(arOrder+1:end-1);
obsVariance = param(end);

%roots of AR polynomial
rAr = abs(roots([-arParam(end:-1:1) 1]))';

%roots of MA polynomial
rMa = abs(roots([maParam(end:-1:1) 1]))';

%nonlinear constraint (inequality)
c = 1.00001 - [rAr rMa];
