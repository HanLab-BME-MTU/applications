function [paramM,errFlag] = extractSignificantParam(param,varCovMat,...
    sigmaPost,numDegFree,alpha)
%EXTRACTSIGNIFICANTPARAM retains the significant part of parameters via a posteriori parameter vector orthogonalization
%
%SYNOPSIS [paramM,errFlag] = extractSignificantParam(param,varCovMat,...
%    sigmaPost,numDegFree,alpha)
%
%INPUT  param     : Parameters to be tested.
%       varCovMat : Variance-covariance matrix of parameters.
%       sigmaPost : A posteriori error standard deviation.
%       numDegFree: Number of degrees of freedom in parameter estimation.
%       alpha     : significance level for hypothesis testing.
%
%OUTPUT paramM    : Parameters with only significant parts.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARK This algorithm is taken from Gaudenz's thesis, p. 133.
%
%Khuloud Jaqaman, April 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramM = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin ~= nargin('extractSignificantParam')
    disp('--extractSignificantParam: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%check size of varCovMat
[numParam,dummy] = size(varCovMat);
if dummy ~= numParam
    disp('--extractSignificantParam: varCovMat should be a square matrix!');
    errFlag = 1;
end

%exit is there are problems with input data
if errFlag
    disp('--extractSignificantParam: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determinability calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the cofactor matrix
cofactMat = varCovMat/sigmaPost^2;

%get eigenvalues and eigenvectors
[eigenVec,eigenVal] = eig(cofactMat);

%transform the parameters
paramTrans = inv(sqrt(eigenVal))*eigenVec'*param;

%get the vector of test statistics
testStatistic = abs(paramTrans/sigmaPost);

%find the p-values of obtaining a value more extreme than the observed
%test-statistics
pValues = 1-tcdf(testStatistic,numDegFree);

%find parameters that are not significantly different from zero
indx = find(pValues > alpha/2);

%Place zero instead of the values of insignificant parameters
paramTrans(indx) = 0;

%get the modified parameter set
paramM = eigenVec*sqrt(eigenVal)*paramTrans;

%%%%% ~~ the end ~~ %%%%%

