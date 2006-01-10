function [c,ceq] = arlsestim1Const(unknown,arOrder,traj,noiseSigma,multi)
%ARLSESTIM1CONST imposes causality constraint on AR coefficients
%
%SYNOPSIS [c,ceq] = arlsestim1Const(unknown,arOrder,traj,noiseSigma,multi)
%
%INPUT  unknown    : Set of AR coefficients and measurement error-free values. 
%       arOrder    : Order of autoregressive model.
%       traj       : Observed trajectory.
%       noiseSigma : Standard deviation of white noise. 
%       multi      : 0/1 if fmincon/fminimax are to be used, respectively.
%       Note that "traj", "noiseSigma" and "multi" are not needed in this
%       function, but must be passed because of fmincon/fminimax requirements.
%
%OUTPUT c   : Constraint that roots should be larger than one.
%       ceq : Empty, but must be included due to matlab requirements.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestim1Const')
    disp('--arlsestim1Const: Incorrect number of input arguments!');
    errFlag  = 1;
    c = [];
    ceq = [];
    return
end

%check input data
[nRow,nCol] = size(unknown);
if nRow ~= arOrder+length(traj(:,1))
    disp('--arlsestim1Const: Wrong length of vector "unknown"!');
    errFlag = 1;
end
if nCol ~= 1
    disp('--arlsestim1Const: "param" should be a column vector!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim1Const: please fix input data!');
    return
end

%distribute parameters
arParam = unknown(1:arOrder)';

%roots of AR polynomial
rAr = abs(roots([-arParam(end:-1:1) 1]))';

%nonlinear constraint (inequality)
c = 1.00001 - rAr;

%nonlinear constraint equality
ceq = [];
