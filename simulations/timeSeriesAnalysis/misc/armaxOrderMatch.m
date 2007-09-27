function [paramVec1,paramVec2,varCovMat1,varCovMat2,errFlag] = ...
    armaxOrderMatch(fitResults1,fitResults2,lengthSeries1,lengthSeries2)
%ARMAXORDERMATCH compensates for order mismatch between 2 ARMAX models.
%
%SYNOPSIS [paramVec1,paramVec2,varCovMat1,varCovMat2,errFlag] = ...
%    armaxOrderMatch(fitResults1,fitResults2,lengthSeries1,lengthSeries2)
%
%INPUT  fitResults1: Structure output from armaxFitKalman. Must contain at
%                    least the following fields:
%           .arParamK     : AR coefficients (row vector).
%           .maParamK     : MA coefficients (row vector).
%           .xParamK      : X coefficients (row vector).
%           .varCovMatF   : Variance-covariance matrix of coefficients.
%
%       fitResults2: same as fitResults1, or a vector with
%               [arOrder, maOrder, xOrder]. In the case of vector input,
%               there will be empty paramVec2 and varCovMat2
%   OPTIONAL
%       lengthSeries1: Length series corresponding to fitResults1.
%       lengthSeries2: Length series correponsding to fitResults2.
%                      Both variables look like the structure input to armaxFitKalman.
%                      If lengthSeries1 and lengthSeries2 are supplied, the
%                      variance-covariance matrices after zero-padding are
%                      calculated from the data. If not, then they will be
%                      simply padded with zeros.
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
if nargin < 2
    disp('--armaxOrderMatch: Wrong number of input arguments!');
    errFlag = 1;
    return
end

if isempty(fitResults1) || isempty(fitResults2) ||...
        length(fitResults1)>1 || (length(fitResults2)>1 && isstruct(fitResults2))
    disp('--armaxOrderMatch: empty or too long input arguments!')
    errFlag = 1;
    return
end

%check 1st model

%check its ARMAX coefficients
if isfield(fitResults1,'arParamK')
    arOrder1 = size(fitResults1.arParamK,2);
else
    disp('--armaxOrderMatch: fitResults1 must have the field arParamK!');
    errFlag = 1;
    return
end
if isfield(fitResults1,'maParamK')
    maOrder1 = size(fitResults1.maParamK,2);
else
    disp('--armaxOrderMatch: fitResults1 must have the field maParamK!');
    errFlag = 1;
    return
end
if isfield(fitResults1,'xParamK')
    xOrder1 = size(fitResults1.xParamK,2);
    xOrder1 = xOrder1 - 1;
else
    disp('--armaxOrderMatch: fitResults1 must have the field xParamK!');
    errFlag = 1;
    return
end

%get the length of the coefficient vector
numParam1 = arOrder1 + maOrder1 + xOrder1 + 1;

%check its variance-covariance matrix
if isfield(fitResults1,'varCovMatF')
    [nRow,nCol] = size(fitResults1.varCovMatF);
    if nRow ~= nCol || nRow ~= numParam1
        disp('--armaxOrderMatch: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
        errFlag = 1;
        return
    end
else
    disp('--armaxOrderMatch: fitResults1 must have the field varCovMatF!');
    errFlag = 1;
    return
end

%check 2nd model

if ~isstruct(fitResults2)
    % check whether we are just augmenting against a maximum order
    if length(fitResults2) == 3
        arOrder2 = fitResults2(1);
        maOrder2 = fitResults2(2);
        xOrder2 = fitResults2(3);

        twoResults = false;
    else
        disp('--armaxOrderMatch: fitResult2 is either a structure or a vector of length 3!')
        errFlag = 1;
    end

else

    % there are two fitResults to augment
    twoResults = true;

    %check its ARMAX coefficients
    if isfield(fitResults2,'arParamK')
        arOrder2 = size(fitResults2.arParamK,2);
    else
        disp('--armaxOrderMatch: fitResults2 must have the field arParamK!');
        errFlag = 1;
        return
    end
    if isfield(fitResults2,'maParamK')
        maOrder2 = size(fitResults2.maParamK,2);
    else
        disp('--armaxOrderMatch: fitResults2 must have the field maParamK!');
        errFlag = 1;
        return
    end
    if isfield(fitResults2,'xParamK')
        xOrder2 = size(fitResults2.xParamK,2);
        xOrder2 = xOrder2 - 1;
    else
        disp('--armaxOrderMatch: fitResults2 must have the field xParamK!');
        errFlag = 1;
        return
    end
end % if isstruct(fitResults2)

%get the length of the coefficient vector
numParam2 = arOrder2 + maOrder2 + xOrder2 + 1;

if twoResults
    
    %check its variance-covariance matrix
    if isfield(fitResults2,'varCovMatF')
        [nRow,nCol] = size(fitResults2.varCovMatF);
        if nRow ~= nCol || nRow ~= numParam2
            disp('--armaxOrderMatch: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
            errFlag = 1;
            return
        end
    else
        disp('--armaxOrderMatch: fitResults2 must have the field varCovMatF!');
        errFlag = 1;
        return
    end

end % if twoResults

%check the length series
if nargin < 3
    lengthSeries1 = [];
end
if nargin < 4
    lengthSeries2 = [];
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

if twoResults
    %fill 2nd parameter vector
    paramVec2 = [fitResults2.arParamK(1,:) zeros(1,max(0,-arOrderDiff)) ...
        fitResults2.maParamK(1,:) zeros(1,max(0,-maOrderDiff)) ...
        fitResults2.xParamK(1,:) zeros(1,max(0,-xOrderDiff))];
end

%determine whether to pad with zeros or whether to use the Fisher
%infotmation matrix
padWithZeros = isempty(lengthSeries1) || isempty(lengthSeries2);

%modify the variance-covariance matrices
if ~padWithZeros

    %calculate the Fisher information matrix of the first padded model
    [fishInfoMat,errFlag] = armaxFisherInfoMatrix(lengthSeries1,[],...
        paramVec1(1:arOrderMax),paramVec1(arOrderMax+1:arOrderMax+maOrderMax),...
        paramVec1(arOrderMax+maOrderMax+1:end),fitResults1.wnVariance);

    %get the variance-covariance matrix of the first padded model
    varCovMat1 = inv(fishInfoMat/fitResults1.numObserve)/fitResults1.numObserve;

    %resort to padding with zeroes of this method fails
    padWithZeros = padWithZeros + (norm(varCovMat1,'fro') > 2);

    if twoResults

        %calculate the Fisher information matrix of the second padded model
        [fishInfoMat,errFlag] = armaxFisherInfoMatrix(lengthSeries2,[],...
            paramVec2(1:arOrderMax),paramVec2(arOrderMax+1:arOrderMax+maOrderMax),...
            paramVec2(arOrderMax+maOrderMax+1:end),fitResults2.wnVariance);

        %get the variance-covariance matrix of the second padded model
        varCovMat2 = inv(fishInfoMat/fitResults2.numObserve)/fitResults2.numObserve;

        %resort to padding with zeroes of this method fails
        padWithZeros = padWithZeros + (norm(varCovMat2,'fro') > 2);

    end

end

if padWithZeros

    %fill 1st variance-covariance matrix
    tmp = [fitResults1.varCovMatF(:,1:arOrder1) zeros(numParam1,max(0,arOrderDiff)) ...
        fitResults1.varCovMatF(:,arOrder1+1:arOrder1+maOrder1) zeros(numParam1,max(0,maOrderDiff)) ...
        fitResults1.varCovMatF(:,arOrder1+maOrder1+1:end) zeros(numParam1,max(0,xOrderDiff))];
    varCovMat1 = [tmp(1:arOrder1,:); zeros(max(0,arOrderDiff),numParamMax); ...
        tmp(arOrder1+1:arOrder1+maOrder1,:); zeros(max(0,maOrderDiff),numParamMax); ...
        tmp(arOrder1+maOrder1+1:end,:); zeros(max(0,xOrderDiff),numParamMax)];

    if twoResults
        %fill 2nd variance-covariance matrix
        tmp = [fitResults2.varCovMatF(:,1:arOrder2) zeros(numParam2,max(0,-arOrderDiff)) ...
            fitResults2.varCovMatF(:,arOrder2+1:arOrder2+maOrder2) zeros(numParam2,max(0,-maOrderDiff)) ...
            fitResults2.varCovMatF(:,arOrder2+maOrder2+1:end) zeros(numParam2,max(0,-xOrderDiff))];
        varCovMat2 = [tmp(1:arOrder2,:); zeros(max(0,-arOrderDiff),numParamMax); ...
            tmp(arOrder2+1:arOrder2+maOrder2,:); zeros(max(0,-maOrderDiff),numParamMax); ...
            tmp(arOrder2+maOrder2+1:end,:); zeros(max(0,-xOrderDiff),numParamMax)];
    end

end

%%%%% ~~ the end ~~ %%%%%


