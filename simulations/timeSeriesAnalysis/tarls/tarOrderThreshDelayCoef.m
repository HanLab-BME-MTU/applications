function [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
        aicV,errFlag] = tarOrderThreshDelayCoef(traj,vThreshTest,delayTest,...
    tarOrderTest,method,tol)
%TARORDERTHRESHDELAYCOEF fits a TAR model of specified segmentation (i.e. determines its AR orders, thresholds, delay parameter and coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
%        aicV,errFlag] = tarOrderThreshDelayCoef(traj,vThreshTest,delayTest,...
%    tarOrderTest,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       vThreshTest  : Matrix containing possible values of thresholds, which are 
%                      sorted in increasing order in a column for each threshold.
%                      Note that the min of a threshold should be larger than the max 
%                      of the previous threshold. Extra entries at the end of a 
%                      column should be indicated with NaN.
%       delayTest    : Row vector of possible values of delay parameter.
%       tarOrderTest : Matrix of possible orders of proposed TAR model in each regime,
%                      where each column corresponds to a different regime. Unnecessary 
%                      entries at the end of a column should be indicated with NaN. 
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
if nargin < 4
    disp('--tarOrderThreshDelayCoef: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check optional parameters
if nargin >= 5
    
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

%get all possible combinations of AR orders
nRegimes = size(tarOrderTest,2); %number of thresholds
for i = 1:nRegimes %number of possible orders in each regime
    numValues(i) = length(find(~isnan(tarOrderTest(:,i))));
end
numComb = prod(numValues); %number of combinations
orderComb = zeros(numComb,nRegimes); %matrix of all combinations
orderComb(:,1) = repeatEntries(tarOrderTest(1:numValues(1),1),numComb/prod(numValues(1))); %entries for 1st regime
for i = 2:nRegimes %entries for rest of regimes
    orderComb(:,i) = repmat(repeatEntries(tarOrderTest(1:numValues(i),i),...
        numComb/prod(numValues(1:i))),prod(numValues(1:i-1)),1);
end

%initial value of Akaikes Information Criterion (AIC)
aicV = 1e20; %ridiculously large number

for i = 1:numComb %go over all order combinations
    
    %get orders for current run
    tarOrder1 = orderComb(i,:)';
    
    %estimate coeffients, residuals, delay parameter and thresholds
    [vThresholds1,delay1,tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,errFlag] =...
        tarThreshDelayCoef(traj,vThreshTest,delayTest,tarOrder1,method,tol);
    if errFlag
        disp('--tarOrderThreshDelayCoef: tarThreshDelayCoef did not function properly!');
        return
    end
    
    %calculate AIC
    aicV1 = fitSet1(1,:)*log(noiseSigma1.^2) + sum(2*tarOrder1+2);
    
    %compare current AIC to AIC in previous orders trial
    %if it is smaller, then update results
    if aicV1 < aicV
        vThresholds = vThresholds1;
        delay = delay1;
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        tarOrder = tarOrder1;
        aicV = aicV1;
    end
    
end %(for i = 1:numComb)

%check for causality of estimated model
for level = 1:nRegimes
    r = abs(roots([-tarParam(level,tarOrder(level):-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--tarOrderThreshDelayCoef: Warning: Predicted model not causal!');
    end
end
