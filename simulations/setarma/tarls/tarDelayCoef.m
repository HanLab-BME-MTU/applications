function [delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] = tarDelayCoef(...
    traj,vThresholds,delayTest,tarOrder,method,tol)
%TARDELAYCOEF fits a TAR model of specified segmentation, AR orders and thresholds (i.e. determines its delay parameter and coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] = tarDelayCoef(...
%    traj,vThresholds,delayTest,tarOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       vThresholds  : Column vector of thresholds, sorted in increasing order.
%       delayTest    : Row vector of values of delay parameter.
%       tarOrder     : Order of proposed TAR model in each regime.
%       method (opt) : Solution method: 'dir' (default) for direct least square 
%                      minimization using the matlab "\", 'iter' for iterative 
%                      refinement using the function "lsIterRefn".
%       tol (opt)    : Tolerance at which calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT delay        : Time lag (delay parameter) of value compared to vThresholds.
%       tarParam     : Estimated parameters in each regime.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise in each regime.
%       fitSet       : Set of points used for data fitting. Each column in 
%                      matrix corresponds to a certain regime. 
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

%initialize output
errFlag = 0;
delay = [];
tarParam = [];
varCovMat = [];
residuals = [];
noiseSigma = [];
fitSet = [];

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--tarDelayCoef: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(delayTest)
    dummy = size(delayTest,1);
    if dummy ~= 1
        disp('--tarDelayCoef: "delayTest" must be a row vector!');
        errFlag = 1;
    end
else
    delayTest = 1;
end
if errFlag
    disp('--tarDelayCoef: Please fix input data!');
    return
end

%check optional parameters
if nargin >= 5
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--tarDelayCoef: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--tarDelayCoef: Warning: "tol" should be positive! A value of 0.001 assigned!');
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

%initial sum of squares of residuals
sumSqResid = 1e20; %ridiculously large number

for delay1 = delayTest %go over all suggested delay parameters
    
    %estimate coeffients and residuals
    [tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,errFlag] = tarCoefEstim(...
        traj,vThresholds,delay1,tarOrder,method,tol);
    if errFlag
        disp('--tarDelayCoef: tarCoefEstim did not function properly!');
        return
    end
    
    %get sum over squares of all residuals
    sumSqResid1 = fitSet1(1,:)*noiseSigma1.^2;
    
    %compare current sum over squared residuals to sum in previous delay parameter trial
    %if it is smaller, then update results
    if sumSqResid1 < sumSqResid
        delay = delay1;
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        sumSqResid = sumSqResid1;
    end
    
end %(for delay = delayTest)

% %check for causality of estimated model
% for level = 1:length(vThresholds)+1
%     r = abs(roots([-tarParam(level,tarOrder(level):-1:1) 1]));
%     if ~isempty(find(r<=1.00001))
%         disp('--tarDelayCoef: Warning: Predicted model not causal!');
%     end
% end
