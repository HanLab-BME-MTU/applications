function [arParam,maParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] =...
    armaCoefEstim(traj,arOrder,maOrder,method,tol)
%ARMACOEFESTIM uses least squares to fit an ARMA model of specified AR and MA orders (i.e. determines its coefficients) to a time series which could have missing data points.
%
%SYNOPSIS [arParam,varCovMat,residuals,noiseSigma,fitSet,errFlag] =...
%    arCoefEstim(traj,arOrder,maOrder,method,tol)
%
%INPUT  traj         : Trajectory to be modeled (with measurement uncertainty).
%                      Missing points should be indicated with NaN.
%       arOrder      : Order of AR part of proposed model.
%       maOrder      : Order of MA part of proposed model.
%       method (opt) : Solution method: 'dir' (default) for direct least square 
%                      minimization using the matlab "\", 'iter' for iterative 
%                      refinement using the function "lsIterRefn".
%       tol (opt)    : Tolerance at which iterative refinement calculation is stopped.
%                      Needed only when method 'iter' is used. If not
%                      supplied, 0.001 is assumed. 
%
%OUTPUT arParam      : Estimated AR parameters in model.
%       maParam      : Estimated MA parameters in model.
%       varCovMat    : Variance-covariance matrix of estimated parameters.
%       residuals    : Difference between measurements and model predictions.
%       noiseSigma   : Estimated standard deviation of white noise.
%       fitSet       : Set of points used for data fitting.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%COMMENTS Algorithm is based on Hannan-Rissanen algorithm
%
%Khuloud Jaqaman, May 2004

%initialize output
errFlag = 0;
arParam = [];
maParam = [];
varCovMat = [];
residuals = [];
noiseSigma = [];

%check if correct number of arguments was used when function was called
if nargin < 3
    disp('--armaCoefEstim: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
[trajLength,nCol] = size(traj);
if trajLength < arOrder+maOrder+1
    disp('--armaCoefEstim: Length of trajectory should be larger than sum of model orders + 1!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--armaCoefEstim: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if arOrder < 0
    disp('--armaCoefEstim: Variable "arOrder" should be >= 0!');
    errFlag = 1;
end
if maOrder < 1
    disp('--armaCoefEstim: Variable "maOrder" should be >= 1!');
    errFlag = 1;
end  
if errFlag
    disp('--armaCoefEstim: Please fix input data!');
    return
end

%check optional parameters
if nargin >= 4
    
    if ~strncmp(method,'dir',3) && ~strncmp(method,'iter',4) 
        disp('--armaCoefEstim: Warning: Wrong input for "method". "dir" assumed!');
        method = 'dir';
    end
    
    if strncmp(method,'iter',4)
        if nargin == 5
            if tol <= 0
                disp('--armaCoefEstim: Warning: "tol" should be positive! A value of 0.001 assigned!');
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
indx = find(~isnan(traj(:,1))); 

%STEP 1: Estimate noise by fitting to a high-order AR model

%determine order of intermediate AR model
fraction = 1-length(indx)/trajLength; %fraction of missing data points
arOrderI = round(20*exp(-4*fraction)); %formula based on some experimentation ...

%find data points to be used in fitting
fitSet = indx(find((indx(arOrderI+1:end)-indx(1:end-arOrderI))==arOrderI)+arOrderI);
fitLength = length(fitSet);

%construct weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
%[size: fitLength by arOrderI]
lhsMat = zeros(fitLength,arOrderI);
for i=1:arOrderI
    lhsMat(:,i) = traj(fitSet-i,1)./traj(fitSet,2);
end

%construct weighted vector of points used to determine parameters (on RHS of equation)
%[size: fitLength by 1]
rhsVec = traj(fitSet,1)./traj(fitSet,2);

%estimate AR coefficients of intermediate model
switch method
    case 'dir' %use MATLAB "\"
        arParamI = (lhsMat\rhsVec)';
    case 'iter' %use iterative refinement method
        [arParamI,errFlag] = lsIterRefn(lhsMat,rhsVec,tol);
        arParamI = arParamI';
end
if errFlag
    disp('--armaCoefEstim: "lsIterRefn" did not function normally!');
    return
end

%compute the noise at each available point
noise = NaN*ones(trajLength,1);
noise(fitSet) = 


%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--armaCoefEstim: Warning: Predicted model (arParam0) not causal!');
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

%add to the beginning of each column the number of data points in that regime
fitSet = [fitLength; fitSet];
