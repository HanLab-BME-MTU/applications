function [tarParam,varCovMat,residuals,noiseSigma,fitSet,delay,vThresholds,...
        aicV,errFlag] = tarModelIdent(traj,modelParam,method,tol)
%TARMODELIDENT fits a TAR model, i.e. determines its segmentation, thresholds, delay parameter, AR orders and coefficients, to a time series which could have missing data points.
%
%SYNOPSIS [tarParam,varCovMat,residuals,noiseSigma,fitSet,delay,vThresholds,...
%        aicV,errFlag] = tarModelIdent(traj,vThreshTest,delayTest,tarOrderTest,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       modelParam   : structure array of models to  be tested. Each entry consists of 
%                      3 elements: vThreshTest, delayTest, tarOrderTest.
%                      See "tarOrderThreshDelayCoef" for descriptions of these variables.
%       method (opt) : Solution method: 'dir' (default) for direct least square 
%                      minimization using the matlab "\", 'iter' for iterative 
%                      refinement using the function "lsIterRefn".
%       tol (opt)    : Tolerance at which calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT tarParam     : Estimated parameters in each regime.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise in each regime.
%       fitSet       : Set of points used for data fitting. Each column in 
%                      matrix corresponds to a certain regime. 
%       delay        : Time lag (delay parameter) of value compared to vThresholds.
%       vThresholds  : Column vector of estimated thresholds, sorted in increasing order.
%       aicV         : Value of Akaike's Information Criterion for the best model.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

tic

errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < 2
    disp('--tarOrderThreshDelayCoef: Incorrect number of input arguments!');
    errFlag  = 1;
    tarParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    fitSet = [];
    delay = [];
    vThresholds = [];
    aicV = []
    return
end

%check optional parameters
if nargin >= 3
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--tarOrderThreshDelayCoef: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--tarOrderThreshDelayCoef: Warning: "tol" should be positive! A value of 0.001 assigned!');
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

for i = 1:length(modelParam) %go over all suggested models
    
    %get model characteristics
    vThreshTest1 = modelParam(i).vThreshTest;
    delayTest1 = modelParam(i).delayTest;
    tarOrderTest1 = modelParam(i).tarOrderTest;
    
    %estimate coeffients, residuals, delay parameter and thresholds
    [tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,delay1,vThresholds1,...
            aicV1,errFlag] = tarOrderThreshDelayCoef(traj,vThreshTest1,delayTest1,...
        tarOrderTest1,method,tol);
    if errFlag
        disp('--tarModelIdent: tarOrderThreshDelayCoef did not function properly!');
        tarParam = [];
        varCovMat = [];
        residuals = [];
        noiseSigma = [];
        fitSet = [];
        delay = [];
        vThresholds = [];
        aicV = [];
        return
    end
    
    %compare current AIC to AIC in previous trials
    %if it is smaller, then update results
    if aicV1 < aicV
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        delay = delay1;
        vThresholds = vThresholds1;
        aicV = aicV1;
    end
    
end %(for i = 1:numComb)

%check for causality of estimated model
for level = 1:length(vThresholds)+1
    tarOrder = length(find(~isnan(tarParam(level,:))));
    r = abs(roots([-tarParam(level,tarOrder:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--tarOrderThreshDelayCoef: Warning: Predicted model not causal!');
    end
end

toc
