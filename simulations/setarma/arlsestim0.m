function [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,arOrder,method,tol)
%ARLSESTIM0 estimates parameters of an AR model using least square fitting
%
%SYNOPSIS [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,arOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainty).
%       arOrder      : Order of proposed AR model.
%       method (opt) : Solution method: 'dir' for direct least square minimization 
%                      using the matlab "\", 'iter' for iterative refinement using 
%                      the function "lsIterRefn". The default is 'dir'.
%       tol (opt)    : Tolerance at which calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT arParam     : Estimated parameters in model.
%       noiseSigma  : Estimated standard deviation of white noise.
%       varCovMat   : Variance-covariance matrix of estimated parameters
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2004

errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < 2
    disp('--arlsestim0: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    noiseSigma = [];
    return
end

%check input data
[trajLength,nCol] = size(traj);
if trajLength < arOrder
    disp('--arlsestim0: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--arlsestim0: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if arOrder < 1
    disp('--arlsestim0: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim0: Please fix input data!');
    arParam = [];
    noiseSigma = [];
    varCovMat = [];
    return
end

%check optional parameters
if nargin >= 3
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--arlsestim0: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--arlsestim0: Warning: "tol" should be positive! A value of 0.001 assigned!');
                tol = 0.001;
            end
        else
            tol = 0.001;
        end
    end
    
else
    method = 'dir';
end

%weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
%[(trajLength-arOrder) by arOrder matrix]
lhsMat = zeros(trajLength-arOrder,arOrder);
for i=1:arOrder
    lhsMat(:,i) = traj(arOrder+1-i:end-i,1)./traj(arOrder+1:end,2);
end

%weighted vector of points used to determine parameters (on RHS of equation)
%[(trajLength-arOrder) by 1 vector]
rhsVec = traj(arOrder+1:end,1)./traj(arOrder+1:end,2);

%estimated values of parameters
switch method
    case 'dir' %use MATLAB "\"
        arParam = (lhsMat\rhsVec)';
    case 'iter' %use iterative refinement method
        [arParam,errFlag] = lsIterRefn(lhsMat,rhsVec,tol);
        arParam = arParam';
end
if errFlag
    disp('--arlsestim0: "lsIterRefn" did not function normally!');
    arParam = [];
    noiseSigma = [];
    varCovMat = [];
    return
end

%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestim0: Warning: Predicted model not causal!');
end

%vector of weighted residuals
epsilon = lhsMat*arParam' - rhsVec;

%variance-covariance matrix
varCovMat = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(trajLength-2*arOrder));

%standard deviation of white noise (note that nonweighted residuals should
%be used in this calculation)
noiseSigma = std(epsilon.*traj(arOrder+1:end,2));
