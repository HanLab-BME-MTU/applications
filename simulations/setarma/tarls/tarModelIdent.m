function [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
        aicV,errFlag] = tarModelIdent(traj,modelParam,method,tol)
%TARMODELIDENT fits a TAR model, i.e. determines its segmentation, thresholds, delay parameter, AR orders and coefficients, to a time series which could have missing data points.
%
%SYNOPSIS [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
%        aicV,errFlag] = tarModelIdent(traj,modelParam,method,tol)
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
%OUTPUT vThresholds  : Column vector of estimated thresholds, sorted in increasing order.
%       delay        : Time lag (delay parameter) of value compared to vThresholds.
%       tarParam     : Estimated parameters in each regime.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise in each regime.
%       fitSet       : Set of points used for data fitting. Each column in 
%                      matrix corresponds to a certain regime. 
%       aicV         : Value of Akaike's Information Criterion for the best model.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

tic

%initialize output
errFlag = 0;
vThresholds = [];
delay = [];
tarParam = [];
varCovMat = [];
residuals = [];
noiseSigma = [];
fitSet = [];
aicV = [];

%check if correct number of arguments was used when function was called
if nargin < 2
    disp('--tarModelIdent: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check optional parameters
if nargin >= 3
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--tarModelIdent: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--tarModelIdent: Warning: "tol" should be positive! A value of 0.001 assigned!');
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

%assign initial value of Akaikes Information Criterion (AIC)
aicV = 1e10; %(ridiculously large number)

for i = 1:length(modelParam) %go over all suggested models
    
    %get model characteristics
    vThreshTest1 = modelParam(i).vThreshTest;
    delayTest1 = modelParam(i).delayTest;
    tarOrderTest1 = modelParam(i).tarOrderTest;
    
    %estimate thresholds, delay parameter, coeffients and residuals
    [vThresholds1,delay1,tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,...
            aicV1,errFlag] = tarOrderThreshDelayCoef(traj,vThreshTest1,delayTest1,...
        tarOrderTest1,method,tol);
    if errFlag
        disp('--tarModelIdent: tarOrderThreshDelayCoef did not function properly!');
        return
    end
    
    %compare current AIC to AIC in previous trials
    %if it is smaller, then update results
    if aicV1 < aicV
        vThresholds = vThresholds1;
        delay = delay1;
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        aicV = aicV1;
    end
    
end %(for i = 1:length(modelParam))

%check whether regimes found are truly distinct. If any two regimes are not
%significantly different from each other, merge them into one

%get number of thresholds
nThresholds = length(vThresholds); 
%get order in each regime
for level=1:nThresholds+1
    tarOrder(level) = length(find(~isnan(tarParam(level,:))));
end

%compare each regime to the one after it
newThresh = [];
for level = 1:nThresholds
    
    %find larger order in the two regimes
    maxOrder = max(tarOrder(level),tarOrder(level+1));
    
    %get estimated parameters
    param1M = tarParam(level,1:maxOrder)';
    param1M(find(isnan(param1M))) = 0;
    param2M = tarParam(level+1,1:maxOrder)';
    param2M(find(isnan(param2M))) = 0;
    
    %get corresponding variance-covariance matrix
    param1V = varCovMat(1:maxOrder,1:maxOrder,level);
    param1V(find(isnan(param1V))) = 0;
    param2V = varCovMat(1:maxOrder,1:maxOrder,level+1);
    param2V(find(isnan(param2V))) = 0;
    
    %calculate vector of differences in coefficients
    diffM = param1M - param2M;
    
    %calculate variance-covariance matrix of difference vector
    diffV = [eye(maxOrder) -eye(maxOrder)]*[param1V zeros(maxOrder); ...
            zeros(maxOrder) param2V]*[eye(maxOrder) -eye(maxOrder)]';
    
    %compute testStatistic whose cumulative distribution function is to be
    %compared with (1-significance level) (notice that testStatistic has a
    %chi2 distribution)
    testStatistic = diffM'*(diffV\diffM)/maxOrder;
    
    %if the cumulative density of this difference is greater than 99%, start a new level
    if chi2cdf(testStatistic,maxOrder) > 0.99
        newThresh(end+1) = level; %i.e. new threshold separates between "level" and "level+1"
    end
    
end %(for level = 1:nThresholds)

%find new orders in new regimes. When levels are combined, give the new
%level an order equal to the smallest order of the combined levels.
if isempty(newThresh) %if all levels are merged
    newOrder(1) = min(tarOrder);
else %if there are at least 2 levels
    newOrder(1) = min(tarOrder(1:newThresh(1)));
    for i=2:length(newThresh)
        newOrder(i) = min(tarOrder(newThresh(i-1)+1:newThresh(i)));
    end
    i = length(newThresh); %if there is only 1 threshold (in which case loop is skipped)
    newOrder(end+1) = min(tarOrder(newThresh(i)+1:end));
end

%assign new thresholds
vThresholds = vThresholds(newThresh);
nThresholdsN = length(vThresholds);

if nThresholdsN < nThresholds %if some levels were merged
    
    %assign new orders
    nThresholds = nThresholdsN;
    tarOrder = newOrder';
    
    %estimate the coefficients one last time using the new thresholds and
    %orders, plus the previously obtained delay parameter
    [tarParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] = tarCoefEstim(...
        traj,vThresholds,delay,tarOrder,method,tol);
    
    %get the new model's AIC value
    [aicV,errFlag] = tarAic(traj(:,1),vThresholds,delay,tarParam);
    
end

%check for causality of estimated model
for level = 1:nThresholds+1
    r = abs(roots([-tarParam(level,tarOrder(level):-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--tarModelIdent: Warning: Predicted model not causal!');
    end
end

toc
