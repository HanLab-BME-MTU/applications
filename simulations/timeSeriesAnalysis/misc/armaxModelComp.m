function [pValueCoef,pValueVar,errFlag] = armaxModelComp(fitResults1,...
    fitResults2,compOpt)
%ARMAXCOEFCOMP tests whether two ARMA models are different
%
%SYNOPSIS [pValueCoef,pValueVar,errFlag] = armaxModelComp(fitResults1,...
%    fitResults2,compOpt)
%
%INPUT  
%   Mandatory
%       fitResults1 : Structure output from armaxFitKalman. Must contain at
%                     least the following fields:
%           .arParamK     : AR coefficients (row vector).
%           .maParamK     : MA coefficients (row vector).
%           .xParamK      : X coefficients (row vector).
%           .varCovMatF   : Variance-covariance matrix of coefficients.
%           .wnVariance   : Estimated variance of white noise in process.
%           .numObserve   : Number of observations used in estimation.
%   Optional
%       fitResults2 : Same as fitResults1, but for 2nd model. Default: All
%                     coefficients are zero and white noise variance is
%                     one.
%       compOpt     : Type of coefficient comparison to be performed:
%                     -'global': Ensemble comparison of all coefficients, 
%                      taking into account coefficient covariances.
%                     -'element': Separate comparison for each coefficient,
%                      ignoring covariances.
%                     -'AR/MA/X': Separate comparison for AR coefficients,
%                       MA cofficients and X coefficients, taking into 
%                       account covariances within each set of coefficients
%                       but ignoring covariances between sets.
%                     Default: 'global'.
%
%OUTPUT pValueCoef  : Probability that difference between coefficients is 
%                     >= difference observed assuming that the null 
%                     hypothesis is true.
%       pValueVar   : Probability that difference between white noise
%                     variances >= difference observed assuming that the
%                     null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pValueCoef = [];
pValueVar  = [];
errFlag    =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin < 1 
    disp('--armaxModelComp: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%check 1st model

%check its ARMAX coefficients
if isfield(fitResults1,'arParamK')
    [nRow,arOrder1] = size(fitResults1.arParamK);
else
    disp('--armaxModelComp: fitResults1 must have the field arParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'maParamK')
    [nRow,maOrder1] = size(fitResults1.maParamK);
else
    disp('--armaxModelComp: fitResults1 must have the field maParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'xParamK')
    [nRow,xOrder1] = size(fitResults1.xParamK);
    xOrder1 = xOrder1 - 1;
else
    disp('--armaxModelComp: fitResults1 must have the field xParamK!');
    errFlag = 1;
end

%get the length of the coefficient vector
numParam1 = arOrder1 + maOrder1 + xOrder1 + 1;

%check its variance-covariance matrix
if isfield(fitResults1,'varCovMatF')
    [nRow,nCol] = size(fitResults1.varCovMatF);
    if nRow ~= nCol || nRow ~= numParam1
        disp('--armaxModelComp: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
        errFlag = 1;
    end
else
    disp('--armaxModelComp: fitResults1 must have the field varCovMatF!');
    errFlag = 1;
end

%check its white noise variance
if ~isfield(fitResults1,'wnVariance')
    disp('--armaxModelComp: fitResults1 must have the field wnVariance!');
    errFlag = 1;
end

%check the number of observations used in its estimation
if isfield(fitResults1,'numObserve')
    numDegFree1 = fitResults1.numObserve - numParam1 - 1;
else
    disp('--armaxModelComp: fitResults1 must have the field numObserve!');
    errFlag = 1;
end
    
%exit if there are problems in mandatory input
if errFlag
    disp('--armaxModelComp: Please fix input data!');
    return
end

%assign default values of optional arguments
arOrder2_def = arOrder1;
maOrder2_def = maOrder1;
xOrder2_def = xOrder1;
numParam2_def = numParam1;
numDegFree2_def = numDegFree1;
fitResults2_def = struct('arParamK',zeros(1,arOrder2_def),...
    'maParamK',zeros(1,maOrder2_def),...
    'xOrderK',zeros(1,xOrder2_def+1),...
    'varCovMatF',1e-10*eye(numParam2_def),'wnVariance',1);
compOpt_def = 'global';

if nargin < 2 || isempty(fitResults2) %if 2nd model was not input

    %assign 2nd model to default
    fitResults2 = fitResults2_def;
    arOrder2 = arOrder2_def;
    maOrder2 = maOrder2_def;
    xOrder2 = xOrder2_def;
    numParam2 = numParam2_def;
    numDegFree2 = numDegFree2_def;
    
else %if user specified a 2nd model
    
    %check 2nd model

    %check its ARMAX coefficients
    if isfield(fitResults2,'arParamK')
        [nRow,arOrder2] = size(fitResults2.arParamK);
    else
        disp('--armaxModelComp: fitResults2 must have the field arParamK!');
        errFlag = 1;
    end
    if isfield(fitResults2,'maParamK')
        [nRow,maOrder2] = size(fitResults2.maParamK);
    else
        disp('--armaxModelComp: fitResults2 must have the field maParamK!');
        errFlag = 1;
    end
    if isfield(fitResults2,'xParamK')
        [nRow,xOrder2] = size(fitResults2.xParamK);
        xOrder2 = xOrder2 - 1;
    else
        disp('--armaxModelComp: fitResults2 must have the field xParamK!');
        errFlag = 1;
    end

    %get the length of the coefficient vector
    numParam2 = arOrder2 + maOrder2 + xOrder2 + 1;

    %check its variance-covariance matrix
    if isfield(fitResults2,'varCovMatF')
        [nRow,nCol] = size(fitResults2.varCovMatF);
        if nRow ~= nCol || nRow ~= numParam2
            disp('--armaxModelComp: fitResults1.varCovMatF should be a square matrix of side length equal to ARorder+MAorder+Xorder+1!');
            errFlag = 1;
        end
    else
        disp('--armaxModelComp: fitResults2 must have the field varCovMatF!');
        errFlag = 1;
    end

    %check its white noise variance
    if ~isfield(fitResults2,'wnVariance')
        disp('--armaxModelComp: fitResults2 must have the field wnVariance!');
        errFlag = 1;
    end

    %check the number of observations used in its estimation
    if isfield(fitResults2,'numObserve')
        numDegFree2 = fitResults2.numObserve - numParam2 - 1;
    else
        disp('--armaxModelComp: fitResults2 must have the field numObserve!');
        errFlag = 1;
    end

end %(if nargin < 2 || isempty(fitResults2) ... else ...)

if nargin  < 3 || isempty(compOpt) %if compOpt was not input

    compOpt = compOpt_def; %assign default

else %if user specified compOpt

    %check that it's a valid option
    if ~strcmp(compOpt,'global') && ~strcmp(compOpt,'element') && ...
            ~strcmp(compOpt,'AR/MA/X')
        disp('--armaxModelComp: "compOpt" should be either "global", "element" or "AR/MA/X"!');
        errFlag = 1;
    end

end %(if nargin < 3 || isempty(compOpt) ... else ...)

%exit if there are problems in input data
if errFlag
    disp('--armaxModelComp: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coefficient comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get coefficient vectors and their variance-covariance matrices after
%matching model orders
[armaxParam1,armaxParam2,varCovMatT1,varCovMatT2,errFlag] = ...
    armaxOrderMatch(fitResults1,fitResults2);
numParamMax = length(armaxParam1);

%get the smaller number of degrees of freedom
numDegFreeMin = min(numDegFree1,numDegFree2);

switch compOpt

    case 'global' %compare all coefficients at the same time

        %compute testStatistic
        [testStatistic,errFlag] = globalCoefTestStat(armaxParam1,armaxParam2,...
            varCovMatT1,varCovMatT2);
        
        %get the p-value of the test statistic assuming an F-distribution
        pValueCoef = 1 - fcdf(testStatistic,numParamMax,numDegFreeMin);

    case 'element' %compare one coefficient at a time

        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %calculate variance-covariance matrix of difference vector
        diffV = diag(varCovMatT1 + varCovMatT2);

        %calculate the test statistic
        testStatistic = diffM./sqrt(diffV);

        %get the p-value assuming that each elemenet in testStatistic
        %follows a student t-distribution
        pValueCoef = 1 - tcdf(abs(testStatistic),numDegFreeMin);

    case 'AR/MA/X' %compare AR, MA and X coefficients alone
     
        H = NaN*ones(1,3);
        pValueCoef = H;
        
        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %AR test
        if arOrderMax ~= 0

            %calculate variance-covariance matrix of AR difference vector
            diffV = varCovMatT1(1:arOrderMax,1:arOrderMax) + ...
                varCovMatT2(1:arOrderMax,1:arOrderMax);

            %compute testStatistic
            testStatistic = diffM(1:arOrderMax)*(diffV\diffM(1:arOrderMax)')/arOrderMax;

            %get the p-value of the test statistic assuming an F-distribution
            pValueCoef(1) = 1 - fcdf(testStatistic,arOrderMax,numDegFreeMin);

        end

        %MA test
        if maOrderMax ~= 0

            %calculate variance-covariance matrix of MA difference vector
            diffV = varCovMatT1(arOrderMax+1:arOrderMax+maOrderMax,...
                arOrderMax+1:arOrderMax+maOrderMax) + ...
                varCovMatT2(arOrderMax+1:arOrderMax+maOrderMax,...
                arOrderMax+1:arOrderMax+maOrderMax);

            %compute testStatistic
            testStatistic = diffM(arOrderMax+1:arOrderMax+maOrderMax)*...
                (diffV\diffM(arOrderMax+1:arOrderMax+maOrderMax)')/maOrderMax;

            %get the p-value of the test statistic assuming an F-distribution
            pValueCoef(2) = 1 - fcdf(testStatistic,maOrderMax,numDegFreeMin);

        end

        %X test
        if xOrderMax ~= -1

            %calculate variance-covariance matrix of X difference vector
            diffV = varCovMatT1(arOrderMax+maOrderMax+1:end,...
                arOrderMax+maOrderMax+1:end) + ...
                varCovMatT2(arOrderMax+maOrderMax+1:end,...
                arOrderMax+maOrderMax+1:end);

            %compute testStatistic
            testStatistic = diffM(arOrderMax+maOrderMax+1:end)*...
                (diffV\diffM(arOrderMax+maOrderMax+1:end)')/(xOrderMax+1);

            %get the p-value of the test statistic assuming a Fisher distribution
            pValueCoef(3) = 1 - fcdf(testStatistic,xOrderMax+1,numDegFreeMin);

        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%White noise variance comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate test-statistic
testStatistic = fitResults1.wnVariance/fitResults2.wnVariance;

if testStatistic > 1
    pValueVar = 1 - fcdf(testStatistic,numDegFree1,numDegFree2);
else
    pValueVar = fcdf(testStatistic,numDegFree1,numDegFree2);
end


%%%%% ~~ the end ~~ %%%%%

