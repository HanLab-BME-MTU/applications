function [arParam,varCovMat,residuals,noiseSigma,errFlag] = arlsestim0(traj,arOrder,method,tol)
%ARLSESTIM0 fits an AR model to a trajectory (which could have missing data points) using least squares.
%
%SYNOPSIS [arParam,varCovMat,residuals,noiseSigma,errFlag] = arlsestim0(traj,arOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainty).
%                      Missing points should be indicated with NaN.
%       arOrder      : Order of proposed AR model.
%       method (opt) : Solution method: 'dir' (default) for direct least square 
%                      minimization using the matlab "\", 'iter' for iterative 
%                      refinement using the function "lsIterRefn".
%       tol (opt)    : Tolerance at which calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT arParam      : Estimated parameters in model.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < 2
    disp('--arlsestim0: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    varCovMat = [];
    residuals = [];
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
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
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

%get indices of available points
indxP = find(~isnan(traj(:,1))); 

%find data points to be used in fitting
fitSet = indxP(find((indxP(arOrder+1:end)-indxP(1:end-arOrder))==arOrder)+arOrder);
fitLength = length(fitSet);

%construct weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
%[size: fitLength by arOrder]
lhsMat = zeros(fitLength,arOrder);
for i=1:arOrder
    lhsMat(:,i) = traj(fitSet-i,1)./traj(fitSet,2);
end

%construct weighted vector of points used to determine parameters (on RHS of equation)
%[size: fitLength by 1]
rhsVec = traj(fitSet,1)./traj(fitSet,2);

%estimate AR coefficients
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
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    return
end


%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestim0: Warning: Predicted model (arParam0) not causal!');
end

%get vector of weighted residuals
residuals = NaN*ones(trajLength,1);
residuals(fitSet) = lhsMat*arParam' - rhsVec;

%calculate variance-covariance matrix
varCovMat = inv(lhsMat'*lhsMat)*(residuals(fitSet)'*residuals(fitSet)/(fitLength-arOrder));

%get nonweighted residuals
residuals(fitSet) = residuals(fitSet).*traj(fitSet,2);

%get standard deviation of white noise
noiseSigma = std(residuals(fitSet));
