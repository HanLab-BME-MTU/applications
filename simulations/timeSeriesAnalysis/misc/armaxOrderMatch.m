function [paramVec1,paramVec2,varCovMat1,varCovMat2,errFlag] = ...
    armaxOrderMatch(fitResults1,fitResults2)
%ARMAXORDERMATCH compensates for order mismatch between 2 ARMAX models.
%
%SYNOPSIS [paramVec1,paramVec2,varCovMat1,varCovMat2,errFlag] = ...
%    armaxOrderMatch(fitResults1,fitResults2)
% 
%INPUT  fitResults1: Structure output from armaxFitKalman. Must contain at
%                    least the following fields:
%           .arParamK     : AR coefficients (row vector).
%           .maParamK     : MA coefficients (row vector).
%           .xParamK      : X coefficients (row vector).
%           .varCovMatF   : Variance-covariance matrix of coefficients.
%
%       fitResults2: same as fitResults1.
%
%OUTPUT paramVec1  : Vector of parameters in model 1, with 0s in place of
%                    parameters that exist in model 2 only.
%       paramVec2  : Vector of parameters in model 2, with 0s in place of
%                    parameters that exist in model 1 only.
%       varCovMat1 : Variance-covariance matrix of paramVec1.
%       varCovMat2 : Variance-covariance matrix of paramVec2.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramVec1 = [];
paramVec2 = [];
varCovMat1 = [];
varCovMat2 = [];
errFlag    =  0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin ~= 2 
    disp('--armaxOrderMatch: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%check 1st model

%check its ARMAX coefficients
if isfield(fitResults1,'arParamK')
    [nRow,arOrder1] = size(fitResults1.arParamK);
else
    disp('--armaxOrderMatch: fitResults1 must have the field arParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'maParamK')
    [nRow,maOrder1] = size(fitResults1.maParamK);
else
    disp('--armaxOrderMatch: fitResults1 must have the field maParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'xParamK')
    [nRow,xOrder1] = size(fitResults1.xParamK);
    xOrder1 = xOrder1 - 1;
else
    disp('--armaxOrderMatch: fitResults1 must have the field xParamK!');
    errFlag = 1;
end

%get the length of the coefficient vector
numParam1 = arOrder1 + maOrder1 + xOrder1 + 1;

%check its variance-covariance matrix
if isfield(fitResults1,'varCovMatF')
    [nRow,nCol] = size(fitResults1.varCovMatF);
    if nRow ~= nCol || nRow ~= numParam1
        disp('--armaxOrderMatch: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
        errFlag = 1;
    end
else
    disp('--armaxOrderMatch: fitResults1 must have the field varCovMatF!');
    errFlag = 1;
end

%check 2nd model

%check its ARMAX coefficients
if isfield(fitResults2,'arParamK')
    [nRow,arOrder2] = size(fitResults2.arParamK);
else
    disp('--armaxOrderMatch: fitResults2 must have the field arParamK!');
    errFlag = 1;
end
if isfield(fitResults2,'maParamK')
    [nRow,maOrder2] = size(fitResults2.maParamK);
else
    disp('--armaxOrderMatch: fitResults2 must have the field maParamK!');
    errFlag = 1;
end
if isfield(fitResults2,'xParamK')
    [nRow,xOrder2] = size(fitResults2.xParamK);
    xOrder2 = xOrder2 - 1;
else
    disp('--armaxOrderMatch: fitResults2 must have the field xParamK!');
    errFlag = 1;
end

%get the length of the coefficient vector
numParam2 = arOrder2 + maOrder2 + xOrder2 + 1;

%check its variance-covariance matrix
if isfield(fitResults2,'varCovMatF')
    [nRow,nCol] = size(fitResults2.varCovMatF);
    if nRow ~= nCol || nRow ~= numParam2
        disp('--armaxOrderMatch: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
        errFlag = 1;
    end
else
    disp('--armaxOrderMatch: fitResults2 must have the field varCovMatF!');
    errFlag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model order matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get largest orders
arOrderMax = max(arOrder1,arOrder2);
maOrderMax = max(maOrder1,maOrder2);
xOrderMax = max(xOrder1,xOrder2);
numParamMax = arOrderMax + maOrderMax + xOrderMax + 1;

%find order difference
arOrderDiff = arOrder2 - arOrder1;
maOrderDiff = maOrder2 - maOrder1;
xOrderDiff  = xOrder2  - xOrder1;

%fill 1st parameter vector
paramVec1 = [fitResults1.arParamK(1,:) zeros(1,max(0,arOrderDiff)) ...
    fitResults1.maParamK(1,:) zeros(1,max(0,maOrderDiff)) ...
    fitResults1.xParamK(1,:) zeros(1,max(0,xOrderDiff))];

%fill 2nd parameter vector
paramVec2 = [fitResults2.arParamK(1,:) zeros(1,max(0,-arOrderDiff)) ...
    fitResults2.maParamK(1,:) zeros(1,max(0,-maOrderDiff)) ...
    fitResults2.xParamK(1,:) zeros(1,max(0,-xOrderDiff))];

%fill 1st variance-covariance matrix
tmp = [fitResults1.varCovMatF(:,1:arOrder1) zeros(numParam1,max(0,arOrderDiff)) ...
    fitResults1.varCovMatF(:,arOrder1+1:arOrder1+maOrder1) zeros(numParam1,max(0,maOrderDiff)) ...
    fitResults1.varCovMatF(:,arOrder1+maOrder1+1:end) zeros(numParam1,max(0,xOrderDiff))];
varCovMat1 = [tmp(1:arOrder1,:); zeros(max(0,arOrderDiff),numParamMax); ...
    tmp(arOrder1+1:arOrder1+maOrder1,:); zeros(max(0,maOrderDiff),numParamMax); ...
    tmp(arOrder1+maOrder1+1:end,:); zeros(max(0,xOrderDiff),numParamMax)];

%fill 2nd variance-covariance matrix
tmp = [fitResults2.varCovMatF(:,1:arOrder2) zeros(numParam2,max(0,-arOrderDiff)) ...
    fitResults2.varCovMatF(:,arOrder2+1:arOrder2+maOrder2) zeros(numParam2,max(0,-maOrderDiff)) ...
    fitResults2.varCovMatF(:,arOrder2+maOrder2+1:end) zeros(numParam2,max(0,-xOrderDiff))];
varCovMat2 = [tmp(1:arOrder2,:); zeros(max(0,-arOrderDiff),numParamMax); ...
    tmp(arOrder2+1:arOrder2+maOrder2,:); zeros(max(0,-maOrderDiff),numParamMax); ...
    tmp(arOrder2+maOrder2+1:end,:); zeros(max(0,-xOrderDiff),numParamMax)];


%%%%% ~~ the end ~~ %%%%%
