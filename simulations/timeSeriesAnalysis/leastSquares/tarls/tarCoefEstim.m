function [tarParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] = tarCoefEstim(...
    traj,vThresholds,delay,tarOrder,method,tol)
%TARCOEFESTIM uses least squares to fit a TAR model of specified segmentation, AR orders, thresholds and delay parameter (i.e. determines its coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [tarParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] = tarCoefEstim(...
%    traj,vThresholds,delay,tarOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainties).
%                      Missing points should be indicated with NaN.
%       vThresholds  : Column vector of thresholds, sorted in increasing order.
%       delay        : Time lag of value compared to vThresholds.
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
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

%initialize output
errFlag = 0;
tarParam = [];
varCovMat = [];
residuals = [];
noiseSigma = [];
fitSet = [];

%check if correct number of arguments was used when function was called
if nargin < 4
    disp('--tarCoefEstim: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
[trajLength,nCol] = size(traj);
if trajLength < max(tarOrder)
    disp('--tarCoefEstim: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--tarCoefEstim: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if ~isempty(vThresholds)
    [nThresholds,nCol] = size(vThresholds);
    if nCol ~= 1
        disp('tarCoefEstim: "vThresholds" should be a column vector!');
    else
        if min(vThresholds(2:end)-vThresholds(1:end-1)) <= 0
            disp('--tarCoefEstim: Entries in "vThresholds" should be sorted in increasing order, with no two elements alike!');
            errFlag = 1;
        end
    end
    if delay <= 0
        disp('--tarCoefEstim: "delay" should be a positive integer!');
        errFlag = 1;
    end
else
    nThresholds = 0;
    delay = 1;
end
[dummy,nCol] = size(tarOrder);
if nCol ~= 1
    disp('--tarCoefEstim: tarOrder should be a column vector!');
    errFlag = 1;
end
if dummy ~= nThresholds + 1
    disp('--tarCoefEstim: Wrong number of entries in "tarOrder"!');
    errFlag = 1;
else
    if min(tarOrder) < 1
        disp('--tarCoefEstim: All entries in "tarOrder" should be >= 1!');
        errFlag = 1;
    end
end
if errFlag
    disp('--tarCoefEstim: Please fix input data!');
    return
end

%check optional parameters
if nargin >= 5
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--tarCoefEstim: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 4
            if tol <= 0
                disp('--tarCoefEstim: Warning: "tol" should be positive! A value of 0.001 assigned!');
                tol = 0.001;
            end
        else
            tol = 0.001;
        end
    end
    
else
    method = 'dir';
end

%put -/+ infinity at ends of thresholds vector
vThresholds = [-Inf; vThresholds; Inf];

%find data points to be used in fitting and classify them in the different
%regimes (steps 1-6)

%1. get indices of available points
indx = find(~isnan(traj(:,1)));

%2. discard first max(tarOrder) trajectory points
indx2 = find(indx>max(tarOrder));
indx2 = indx2(find(indx2>max(tarOrder)));

%3. find all points whose "delay"-time-steps predecessors are available
indx2 = indx2(find(indx(indx2)>delay));
indx2 = indx2(find(~isnan(traj(indx(indx2)-delay,1))));

%4. classify points into the different regimes
indxClass = zeros(trajLength,nThresholds+1); %initialize vector
for level = 1:nThresholds+1 %go over all levels (regimes)
    temp = indx2(find((traj(indx(indx2)-delay,1)>vThresholds(level)) + ...
        (traj(indx(indx2)-delay,1)<=vThresholds(level+1)) == 2)); %find all points in this level
    fitLength(level) = length(temp);
    indxClass(1:fitLength(level),level) = temp;
end
%remove empty entries in indxClass
indxClass = indxClass(1:max(fitLength),:);
%delete variable indx2
clear indx2;

%5. find all points whose tarOrder(level) previous points are available
fitSet = zeros(size(indxClass));
for level = 1:nThresholds+1
    temp = indx(indxClass(1:fitLength(level),level))-indx(indxClass(1:fitLength(level),level)-tarOrder(level));
    temp = indx(indxClass(find(temp==tarOrder(level)),level));
    fitLength(level) = length(temp);
    fitSet(1:fitLength(level),level) = temp;
end
%remove empty entries in fitSet
fitSet = fitSet(1:max(fitLength),:);
%delete variable temp
clear temp indxClass;

%exit if any of the levels does not have a sufficient number of point
if ~isempty(find(fitLength <= tarOrder'))
    disp('--tarCoefEstim: Bad segmentation of data - cannot determine coefficients!');
    errFlag = 1;
    return
end

%initialize residuals vector (needed later)
residuals = NaN*ones(trajLength,1);

%reserve memory for tarParam and varCovMat
tarParam = NaN*ones(nThresholds+1,max(tarOrder));
varCovMat = NaN*ones(max(tarOrder),max(tarOrder),nThresholds+1);

%estimate parameters in each level
for level = 1:nThresholds+1
    
    %construct weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
    %[size: fitLength by tarOrder(level)]
    lhsMat = zeros(fitLength(level),tarOrder(level));
    for i=1:tarOrder(level)
        lhsMat(:,i) = traj(fitSet(1:fitLength(level),level)-i,1)...
            ./traj(fitSet(1:fitLength(level),level),2);
    end
    
    %construct weighted vector of points used to determine parameters (on RHS of equation)
    %[size: fitLength by 1]
    rhsVec = traj(fitSet(1:fitLength(level),level),1)...
        ./traj(fitSet(1:fitLength(level),level),2);
    
    %estimate AR coefficients
    switch method
        case 'dir' %use MATLAB "\"
            tarParam(level,1:tarOrder(level)) = (lhsMat\rhsVec)';
        case 'iter' %use iterative refinement method
            [tarParam1,errFlag] = lsIterRefn(lhsMat,rhsVec,tol);
            tarParam(level,1:tarOrder(level)) = tarParam1';
    end
    if errFlag
        disp('--tarCoefEstim: "lsIterRefn" did not function normally!');
        return
    end
    
%     %check for causality of estimated model
%     r = abs(roots([-tarParam(level,tarOrder(level):-1:1) 1]));
%     if ~isempty(find(r<=1.00001))
%         disp('--tarCoefEstim: Warning: Predicted model not causal!');
%     end
    
    %get vector of weighted residuals
    residuals(fitSet(1:fitLength(level),level)) = lhsMat*tarParam(level,1:tarOrder(level))' - rhsVec;
    
    %calculate variance-covariance matrix
    varCovMat(1:tarOrder(level),1:tarOrder(level),level) = inv(lhsMat'*lhsMat)*...
        (residuals(fitSet(1:fitLength(level),level))'...
        *residuals(fitSet(1:fitLength(level),level))/(fitLength(level)-tarOrder(level)));
    
    %get nonweighted residuals
    residuals(fitSet(1:fitLength(level),level)) = residuals(fitSet(1:fitLength(level),level))...
        .*traj(fitSet(1:fitLength(level),level),2);
    
    %get standard deviation of white noise
    noiseSigma(level) = std(residuals(fitSet(1:fitLength(level),level)));
    
end %(for levels=1:nThresholds+1)

%add to the beginning of each column the number of data points in that regime
fitSet = [fitLength; fitSet];

%transform noiseSigma into a column vector
noiseSigma = noiseSigma';
