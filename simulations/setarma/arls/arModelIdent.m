function [arParam,varCovMat,residuals,noiseSigma,fitSet,aicV,errFlag] =...
    arModelIdent(traj,arOrderTest,method,tol)
%ARMODELIDENT fits an AR model, i.e. determines its order and coefficients, to a time series which could have missing data points.
%
%SYNOPSIS [arParam,varCovMat,residuals,noiseSigma,fitSet,aicV,errFlag] =...
%    arModelIdent(traj,arOrderTest,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       arOrderTest  : Row vector of possible orders of proposed AR model 
%       method (opt) : Solution method: 'dir' (default) for direct least square 
%                      minimization using the matlab "\", 'iter' for iterative 
%                      refinement using the function "lsIterRefn".
%       tol (opt)    : Tolerance at which calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT arParam      : Estimated parameters in each regime.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise in each regime.
%       fitSet       : Set of points used for data fitting. Each column in 
%                      matrix corresponds to a certain regime. 
%       aicV         : Value of Akaike's Information Criterion for the best model.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < 2
    disp('--arModelIdent: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    fitSet = [];
    aicV = []
    return
end

%check input data
dummy = size(arOrderTest,1);
if dummy ~= 1
    disp('--arModelIdent: arOrderTest should be a row vector!');
    errFlag;
end

%exit if there are problems with input data
if errFlag
    disp('--arModelIdent: Please fix input data!');
    arParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    fitSet = [];
    aicV = []
    return
end

%check optional parameters
if nargin >= 3
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--arModelIdent: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--arModelIdent: Warning: "tol" should be positive! A value of 0.001 assigned!');
                tol = 0.001;
            end
        else
            tol = 0.001;
        end
    end
    
else
    method = 'dir';
    tol = [];
end

%initial value of Akaikes Information Criterion (AIC)
aicV = 1e20; %ridiculously large number

for arOrder1 = arOrderTest %go over all AR orders
    
    %estimate coeffients and residuals
    [arParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,errFlag] =...
        arlsestim0(traj,arOrder1,method,tol);
    if errFlag
        disp('--arModelIdent: arlsestim0 did not function properly!');
        arParam = [];
        varCovMat = [];
        residuals = [];
        noiseSigma = [];
        fitSet = [];
        aicV = []
        return
    end
    
    %calculate AIC
    aicV1 = fitSet1(1,:)*log(noiseSigma1.^2) + 2*arOrder1 + 2;
    
    %compare current AIC to AIC in previous orders trial
    %if it is smaller, then update results
    if aicV1 < aicV
        arParam = arParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        arOrder = arOrder1;
        aicV = aicV1;
    end
    
end %(for i = 1:numComb)

%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arModelIdent: Warning: Predicted model not causal!');
end
