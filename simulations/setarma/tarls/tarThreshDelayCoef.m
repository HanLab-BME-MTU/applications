function [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
        errFlag] = tarThreshDelayCoef(traj,vThreshTest,delayTest,tarOrder,method,tol)
%TARTHRESHDELAYCOEF fits a TAR model of specified segmentation and AR orders (i.e. determines its thresholds, delay parameter and coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [vThresholds,delay,tarParam,varCovMat,residuals,noiseSigma,fitSet,...
%        errFlag] = tarThreshDelayCoef(traj,vThreshTest,delayTest,tarOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       vThreshTest  : Matrix containing possible values of thresholds, which are 
%                      sorted in increasing order in a column for each threshold.
%                      Note that the min of a threshold should be larger than the max 
%                      of the previous threshold. Extra entries at the end of a 
%                      column should be indicated with NaN.
%       delayTest    : Row vector of values of delay parameter.
%       tarOrder     : Order of proposed TAR model in each regime.
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

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--tarThreshDelayCoef: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(vThreshTest)
    if min(vThreshTest(2:end)-vThreshTest(1:end-1)) <= 0
        disp('--tarThreshDelayCoef: Entries in "vThreshTest" should be sorted in increasing order, column by column!');
        errFlag = 1;
    end
end
if errFlag
    disp('--tarThreshDelayCoef: Please fix input data!');
    return
end

%check optional parameters
if nargin >= 5
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--tarThreshDelayCoef: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--tarThreshDelayCoef: Warning: "tol" should be positive! A value of 0.001 assigned!');
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

%get all possible combinations for threshold values
nThresholds = size(vThreshTest,2); %number of thresholds
if nThresholds ~= 0
    for i = 1:nThresholds %number of possible values of each threshold
        numValues(i) = length(find(~isnan(vThreshTest(:,i))));
    end
    numComb = prod(numValues); %number of combinations
    threshComb = zeros(numComb,nThresholds); %matrix of all combinations
    threshComb(:,1) = repeatEntries(vThreshTest(1:numValues(1),1),numComb/prod(numValues(1))); %entries for 1st threshold
    for i = 2:nThresholds %entries for rest of thresholds
        threshComb(:,i) = repmat(repeatEntries(vThreshTest(1:numValues(i),i),...
            numComb/prod(numValues(1:i))),prod(numValues(1:i-1)),1);
    end
else
   numComb = 1;
end

%initial sum of squares of residuals
sumSqResid = 1e20; %ridiculously large number

for i = 1:numComb %go over all threshold combinations
    
    %get thresholds for current run
    if nThresholds ~= 0
        vThresholds1 = threshComb(i,:)';
    else
        vThresholds1 = [];
    end
    
    %estimate coeffients, residuals and delay parameter
    [delay1,tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,errFlag] = ...
        tarDelayCoef(traj,vThresholds1,delayTest,tarOrder,method,tol);
    if errFlag
        disp('--tarThreshDelayCoef: tarDelayCoef did not function properly!');
        return
    end
    
    %get sum over squares of all residuals
    sumSqResid1 = fitSet1(1,:)*noiseSigma1.^2;
    
    %compare current sum over squared residuals to sum in previous thresholds trial
    %if it is smaller, then update results
    if sumSqResid1 < sumSqResid
        vThresholds = vThresholds1;
        delay = delay1;
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        sumSqResid = sumSqResid1;
    end
    
end %(for i = 1:numComb)

% %check for causality of estimated model
% for level = 1:nThresholds+1
%     r = abs(roots([-tarParam(level,tarOrder(level):-1:1) 1]));
%     if ~isempty(find(r<=1.00001))
%         disp('--tarThreshDelayCoef: Warning: Predicted model not causal!');
%     end
% end
