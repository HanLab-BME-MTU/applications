function [H,errFlag] = armaCoefComp(armaCoef1,varCovMat1,armaCoef2,...
    varCovMat2,significance)
%ARMACOEFCOMP tests whether the the ARMA coefficients of 2 models are identical

%SYNOPSIS [H,errFlag] = armaCoefComp(armaCoef1,varCovMat1,armaCoef2,...
%    varCovMat2,significance)
%
%INPUT  armaCoef1   : Structure containing ARMA coefficients of 1st model:
%           .arParam: AR coefficients (row vector).
%           .maParam: MA coefficients (row vector).
%       varCovMat1  : Variance-covariance matrix of ARMA coeffients of 1st
%                     model.
%       armaCoef2   : Same as armaCoef1, but for 2nd model.
%       varCovMat2  : Same as varCovMat1, but for 2nd model.
%       significance: Significance level of hypothesis test. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARK 
%
%Khuloud Jaqaman, September 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin < 4
    disp('--armaCoefComp: You must input ARMA coefficients and the variance-covariance matrices!');
    errFlag = 1;
    return
end

[nRow,arOrder1] = size(armaCoef1.arParam);
if nRow ~= 1
    disp('--armaCoefComp: armaCoef1.arParam should be a row vector!');
    errFlag = 1;
end
[nRow,maOrder1] = size(armaCoef1.maParam);
if nRow ~= 1
    disp('--armaCoefComp: armaCoef1.maParam should be a row vector!');
    errFlag = 1;
end
[nRow,nCol] = size(varCovMat1);
if nRow ~= arOrder1 + maOrder1  || nCol ~= arOrder1 + maOrder1
    disp('--armaCoefComp: varCovMat1 should be a square matrix of side length equal to the sum of the AR and MA orders of model 1!');
    errFlag = 1;
end

[nRow,arOrder2] = size(armaCoef2.arParam);
if nRow ~= 1
    disp('--armaCoefComp: armaCoef2.arParam should be a row vector!');
    errFlag = 1;
end
[nRow,maOrder2] = size(armaCoef2.maParam);
if nRow ~= 1
    disp('--armaCoefComp: armaCoef2.maParam should be a row vector!');
    errFlag = 1;
end
[nRow,nCol] = size(varCovMat2);
if nRow ~= arOrder2 + maOrder2  || nCol ~= arOrder2 + maOrder2
    disp('--armaCoefComp: varCovMat2 should be a square matrix of side length equal to the sum of the AR and MA orders of model 2!');
    errFlag = 1;
end

%check significance and assign default value, if needed
if nargin < 5
    significance = 0.05;
else
    if significance < 0 || significance > 1
        disp('--armaCoefComp: "significance should be between 0 and 1!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaCoefComp: Please fix input data!');
    return
end

%make sure that variance-covariance matrices are not identically zero
if isempty(nonzeros(varCovMat1))
    varCovMat1 = 1e-10*eye(arOrder1+maOrder1);
end
if isempty(nonzeros(varCovMat2))
    varCovMat2 = 1e-10*eye(arOrder2+maOrder2);
end

%compare arOrder1 to arOrder2
if arOrder1 < arOrder2 %add zeroes to armaCoef1.arParam and to varCovMat1
    armaCoef1.arParam = [armaCoef1.arParam zeros(1,arOrder2-arOrder1)];
    varCovMat1 = blkdiag(varCovMat1(1:arOrder1,1:arOrder1),...
        min(abs(nonzeros(varCovMat1)))*eye(arOrder2-arOrder1),...
        varCovMat1(arOrder1+1:end,arOrder1+1:end));
elseif arOrder1 > arOrder2 %add zeroes to armaCoef2.arParam and to varCovMat2
    armaCoef2.arParam = [armaCoef2.arParam zeros(1,arOrder1-arOrder2)];
    varCovMat2 = blkdiag(varCovMat2(1:arOrder2,1:arOrder2),...
        min(abs(nonzeros(varCovMat2)))*eye(arOrder1-arOrder2),...
        varCovMat2(arOrder2+1:end,arOrder2+1:end));
end

%compare maOrder1 to maOrder2
if maOrder1 < maOrder2 %add zeroes to armaCoef1.maParam and to varCovMat1
    armaCoef1.maParam = [armaCoef1.maParam zeros(1,maOrder2-maOrder1)];
    varCovMat1 = blkdiag(varCovMat1,min(abs(nonzeros(varCovMat1)))*...
        eye(maOrder2-maOrder1));
elseif maOrder1 > maOrder2
    armaCoef2.maParam = [armaCoef2.maParam zeros(1,maOrder1-maOrder2)];
    varCovMat2 = blkdiag(varCovMat2,min(abs(nonzeros(varCovMat2)))*...
        eye(maOrder1-maOrder2));
end

%combine AR and MA coefficients for each model
armaParam1 = [armaCoef1.arParam armaCoef1.maParam];
armaParam2 = [armaCoef2.arParam armaCoef2.maParam];

%get "combined" ARMA order
combOrder = length(armaParam1); %=length(armaParam2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hypothesis testing

%Null hypothesis: The two models are the same, i.e. their coefficients are identical
%Alternative hypothesis: Model coefficients are not identical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate vector of differences in coefficients
diffM = armaParam1 - armaParam2;

%calculate variance-covariance matrix of difference vector
diffV = [eye(combOrder) -eye(combOrder)]*[varCovMat1 zeros(combOrder); ...
    zeros(combOrder) varCovMat2]*[eye(combOrder) -eye(combOrder)]';

%compute testStatistic
testStatistic = diffM*(diffV\diffM')/combOrder;

%get the p-value of the test statistic assuming a chi2 distribution
pValue = 1 - chi2cdf(testStatistic,combOrder);

if pValue < significance %if p-value is smaller than probability of type I error
    H = 1; %reject null hypothesis that series is IID
else %if p-value is larger than probability of type I error
    H = 0; %cannot reject null hypothesis
end


%%%%% ~~ the end ~~ %%%%%

