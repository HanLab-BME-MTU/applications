function [tarParam,varCovMat,residuals,noiseSigma,fitSet,delay,vThresholds,...
        errFlag] = tarThreshDelayCoef(traj,vThreshTest,delayTest,tarOrder,method,tol)
%TARTHRESHDELAYCOEF fits a TAR model of specified segmentation, AR orders and thresholds (i.e. determines its delay parameter and coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [tarParam,varCovMat,residuals,noiseSigma,fitSet,delay,errFlag] = tarThreshDelayCoef(...
%    traj,vThresholds,delayTest,tarOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       vThreshTest  : 2-column matrix containing min and max values of thresholds
%                      Thresholds are sorted in increasing order in a row.
%                      Note that the max of a threshold should be larger
%                      than its min and smaller than the min of the next threshold.
%       delayTest    : Values of delay parameter to be tested.
%       tarOrder     : Order of proposed TAR model in each regime.
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
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--tarThreshDelayCoef: Incorrect number of input arguments!');
    errFlag  = 1;
    tarParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    fitSet = [];
    delay = [];
    vThresholds = [];
    return
end

%check input data
[nThresholds,dummy] = size(vThreshTest,2);
if dummy ~= 2
    disp('--tarThreshDelayCoef: "vThreshTest" must have two columns!');
    errFlag = 1;
else
    if min(vThreshTest(1:end,1)-vThreshTest(1:end,1)) <= 0
        disp('--tarThreshDelayCoef: Entries in 2nd column in "vThreshTest" should be larger than corresponding entries in 1st column!');
        errFlag = 1;
    end
    if min(vThreshTest(2:end,1)-vThreshTest(1:end-1,2)) < 0
        disp('--tarThreshDelayCoef: Min possible value of a certain threshold should not be smaller than the max possible value of a previous threshold!');
        errFlag = 1;
    end
end
if errFlag
    disp('--tarThreshDelayCoef: Please fix input data!');
    tarParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    fitSet = [];
    delay = [];
    vThresholds = [];
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

numPermute = nThresholds^rangeDiv
permutation = zeros(numPermute,nThresholds);
permutation(:,1) = repeatEntries([1:rangeDiv]',nThresholds);
permutation(:,2) = 


%initial sum of squares of residuals
sumSqResid = 1e20; %ridiculously large number


for delay1 = delayTest %go over all suggested delay parameters
    
    %estimate coeffients and residuals
    [tarParam1,varCovMat1,residuals1,noiseSigma1,fitSet1,errFlag] = tarlsestim0(...
        traj,vThresholds,delay1,tarOrder,method,tol);
    if errFlag
        disp('--tarThreshDelayCoef: tarlsestim0 did not function properly!');
        tarParam = [];
        varCovMat = [];
        residuals = [];
        noiseSigma = [];
        fitSet = [];
        delay = [];
        vThresholds = [];
        return
    end
    
    %get sum over squares of all residuals
    sumSqResid1 = fitSet1(1,:)*noiseSigma1';
    
    %compare current sum over squared residuals to sum in previous delay parameter trial
    %if it is smaller, then update results
    if sumSqResid1 < sumSqResid
        tarParam = tarParam1;
        varCovMat = varCovMat1;
        residuals = residuals1;
        noiseSigma = noiseSigma1;
        fitSet = fitSet1;
        delay = delay1;
        sumSqResid = sumSqResid1;
    end
    
end %(for delay = delayTest)