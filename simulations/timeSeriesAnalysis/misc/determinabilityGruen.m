function [rho,errFlag] = determinabilityGruen(varCovMat)
%RELDETERMINABILITY calculates the determinability of parameters in a model
%
%SYNOPSIS [rho,errFlag] = determinabilityGruen(varCovMat)
%
%INPUT  varCovMat: Variance-covariance matrix of estimated parameters.
%
%OUTPUT rho      : Vector of relative parameter determinability.
%       errFlag  : 0 if function executes normally, 1 otherwise.
%
%REMARK Formula used here is taken from Gaudenz's dissertation, p. 132,
%       after Eq. A.17.
%
%Khuloud Jaqaman, September 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin ~= 1
    disp('--determinabilityGruen: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%check size of varCovMat
[numParam,dummy] = size(varCovMat);
if dummy ~= numParam
    disp('--determinabilityGruen: varCovMat should be a square matrix!');
    errFlag = 1;
end

%exit is there are problems with input data
if errFlag
    disp('--determinabilityGruen: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determinability calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize relative determinability measure
rho = NaN*ones(1,numParam);

%calculate trace of matrix for denominator
denom = trace(varCovMat);

for i=1:numParam
    rho(i) = sum(varCovMat(i,:).^2)/varCovMat(i,i)/denom;
end


%%%%% ~~ the end ~~ %%%%%

