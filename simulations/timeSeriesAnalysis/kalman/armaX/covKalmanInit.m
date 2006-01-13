function [stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
    arOrder,maOrder,maxOrder)
%COVKALMANINIT calculates the covariance matrix of the initial state used in Kalman recursions.
%
%SYNOPSIS [stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
%    arOrder,maOrder,maxOrder)
%
%INPUT  arParam       : Autoregressive coefficients (row vector).
%       maParam       : Moving average coefficients (row vector).
%       procErrCov: Vector G (Eq. 2.15, 2.16) which is equivalent to the
%                       covariance of the process with the white noise
%                       (column vector).
%       arOrder       : Order of AR part of proposed ARMA model.
%       maOrder       : Order of MA part of proposed ARMA model.
%       maxOrder      : max(arOrder,maOrder+1).
%
%OUTPUT stateCovMat00: Covariance matrix of initial state.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980), Section 4.
%        All equation numbers used are those in Jones' paper.
%
%Khuloud Jaqaman, July 2004

%initialize output
stateCovMat00 = [];
errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < nargin('covKalmanInit')
    disp('--covKalmanInit: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%solve the set of linear equations (4.11) to obtain the trajectory autocovariances

%first, calculate vector on right hand side of Eq. 4.11
maParamMod = [1 maParam];
rhsVec = zeros(arOrder+1,1);
for i = 1:min(maOrder,arOrder)+1
    rhsVec(i) = maParamMod(i:maOrder+1)*procErrCov(1:maOrder+2-i);
end

%second, construct left hand side matrix, without taking into account that
%C(-i) = C(i) where C = correlation
lhsMat = zeros(arOrder+1,2*arOrder+1);
for i = 1:arOrder+1
    lhsMat(i,i:arOrder+i) = [-arParam(end:-1:1) 1];
end

%then, combine matrix elements using the fact that C(-i) = C(i), to get
%final form of matrix on left hand side
for i = 1:arOrder
    lhsMat(:,end-i+1) = lhsMat(:,end-i+1) + lhsMat(:,i);
end
lhsMat = lhsMat(:,arOrder+1:end);

%now solve the system of linear equations for the correlations
correlations = lhsMat\rhsVec;

%if maOrder > arOrder, calculate the rest of the correlations (Eq. 4.9)
if maOrder > arOrder
    correlations = [correlations; zeros(maOrder-arOrder,1)];
    for i = arOrder+2:maOrder+1
        if ~isempty(arParam)
            correlations(i) = arParam*correlations(i-1:-1:i-arOrder);
        end
        if ~isempty(maParam)
            correlations(i) = correlations(i) + maParam(i-1:maOrder)...
                *procErrCov(1:maOrder+2-i);
        end
    end
end

%Evaluate covariance matrix
stateCovMat00 = zeros(maxOrder);
stateCovMat00(1,:) = correlations(1:maxOrder)';
for i = 2:maxOrder
    for j = i:maxOrder
        stateCovMat00(i,j) = correlations(j-i+1) - procErrCov(1:i-1)'...
            *procErrCov(1+j-i:j-1);
    end
end
stateCovMat00 = stateCovMat00 + stateCovMat00' - diag(diag(stateCovMat00));


%%%%% ~~ the end ~~ %%%%%
