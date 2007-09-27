function [mLogPValueCoef,mLogPValueVar,errFlag] = armaxModelComp(fitResults1,...
    fitResults2,compOpt,lengthSeries1,lengthSeries2)
%ARMAXCOEFCOMP tests whether two ARMA models are different
%
%SYNOPSIS [mLogPValueCoef,mLogPValueVar,errFlag] = armaxModelComp(fitResults1,...
%    fitResults2,compOpt,lengthSeries1,lengthSeries2)
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
%       lengthSeries1: Length series described by fitResults1.
%       lengthSeries2: Length series described by fitResults2.
%                      Both variables look like the structure input to armaxFitKalman.
%                      If lengthSeries1 and lengthSeries2 are supplied, the
%                      variance-covariance matrices after zero-padding are
%                      calculated from the data. If not, then they will be
%                      simply padded with zeros.
%
%OUTPUT mLogPValueCoef: Probability that difference between coefficients is 
%                       >= difference observed assuming that the null 
%                       hypothesis is true.
%       mLogPValueVar : Probability that difference between white noise
%                       variances >= difference observed assuming that the
%                       null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mLogPValueCoef = [];
mLogPValueVar  = [];
errFlag = 0;

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
    arOrder1 = size(fitResults1.arParamK,2);
else
    disp('--armaxModelComp: fitResults1 must have the field arParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'maParamK')
    maOrder1 = size(fitResults1.maParamK,2);
else
    disp('--armaxModelComp: fitResults1 must have the field maParamK!');
    errFlag = 1;
end
if isfield(fitResults1,'xParamK')
    xOrder1 = size(fitResults1.xParamK,2);
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
lengthSeries1_def = [];
lengthSeries2_def = [];

if nargin < 2 || isempty(fitResults2) %if 2nd model was not input

    %assign 2nd model to default
    fitResults2 = fitResults2_def;
    %     arOrder2 = arOrder2_def;
    %     maOrder2 = maOrder2_def;
    %     xOrder2 = xOrder2_def;
    %     numParam2 = numParam2_def;
    numDegFree2 = numDegFree2_def;
    
else %if user specified a 2nd model
    
    %check 2nd model

    %check its ARMAX coefficients
    if isfield(fitResults2,'arParamK')
        arOrder2 = size(fitResults2.arParamK,2);
    else
        disp('--armaxModelComp: fitResults2 must have the field arParamK!');
        errFlag = 1;
    end
    if isfield(fitResults2,'maParamK')
        maOrder2 = size(fitResults2.maParamK,2);
    else
        disp('--armaxModelComp: fitResults2 must have the field maParamK!');
        errFlag = 1;
    end
    if isfield(fitResults2,'xParamK')
        xOrder2 = size(fitResults2.xParamK,2);
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

%check whether the length series were input
if nargin < 4 || isempty(lengthSeries1)
    lengthSeries1 = lengthSeries1_def;
end
if nargin < 5 || isempty(lengthSeries2)
    lengthSeries2 = lengthSeries2_def;
end
  
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
    armaxOrderMatch(fitResults1,fitResults2,lengthSeries1,lengthSeries2);
numParamMax = length(armaxParam1);

%get the smaller number of degrees of freedom
numDegFreeMin = min(numDegFree1,numDegFree2);

switch compOpt

    case 'global' %compare all coefficients at the same time

        %compute testStatistic
        [testStatistic,errFlag] = globalCoefTestStat(armaxParam1,armaxParam2,...
            varCovMatT1,varCovMatT2);
        
        %get -log10(p-value) of the test statistic assuming a chi-square distribution
%         mLogPValueCoef = -log10(fcdf(1/testStatistic,numDegFreeMin,numParamMax));
        mLogPValueCoef = chi2cdfExtrapolateLog(testStatistic,numParamMax);

    case 'element' %compare one coefficient at a time

        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %calculate variance-covariance matrix of difference vector
        diffV = diag(varCovMatT1 + varCovMatT2);

        %calculate the test statistic
        testStatistic = diffM./sqrt(diffV);

        %get -log(p-value) assuming that each elemenet in testStatistic
        %follows a student t-distribution
        mLogPValueCoef = -log10(1-tcdf(abs(testStatistic),numDegFreeMin));

    case 'AR/MA/X' %compare AR, MA and X coefficients alone
     
        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %AR test
        if arOrderMax ~= 0

            %calculate variance-covariance matrix of AR difference vector
            diffV = varCovMatT1(1:arOrderMax,1:arOrderMax) + ...
                varCovMatT2(1:arOrderMax,1:arOrderMax);

            %compute testStatistic
            testStatistic = diffM(1:arOrderMax)*(diffV\diffM(1:arOrderMax)');
%             /arOrderMax;
            
            %get -log(p-value) of the test statistic assuming a chi-square distribution
%             mLogPValueCoef(1) = -log10(fcdf(1/testStatistic,numDegFreeMin,arOrderMax));
            mLogPValueCoef(1) = chi2cdfExtrapolateLog(testStatistic,arOrderMax);

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
                (diffV\diffM(arOrderMax+1:arOrderMax+maOrderMax)');
%             /maOrderMax;

            %get -log(p-value) of the test statistic assuming a chi-square distribution
%             mLogPValueCoef(2) = -log10(fcdf(1/testStatistic,numDegFreeMin,maOrderMax));
            mLogPValueCoef(2) = chi2cdfExtrapolateLog(testStatistic,maOrderMax);

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
                (diffV\diffM(arOrderMax+maOrderMax+1:end)');
%             /(xOrderMax+1);

            %get -log(p-value) of the test statistic assuming a chi-square distribution
%             mLogPValueCoef(3) = -log10(fcdf(1/testStatistic,numDegFreeMin,xOrderMax+1));
            mLogPValueCoef(3) = chi2cdfExtrapolateLog(testStatistic,xOrderMax+1);

        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%White noise variance comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate test-statistic
testStatistic = fitResults1.wnVariance/fitResults2.wnVariance;

%get -log(p-value) of the test statistic assuming an F-distribution
mLogPValueVar = fcdfExtrapolateLog(testStatistic,numDegFree1,numDegFree2);

%shift by -log2 in order to assign a distance of 0 to variances that are exactly the same
%do not allow negative values, which happen due to numerical approximations
%when variances are very, very close to each other.
mLogPValueVar = max(0,mLogPValueVar-log10(2));


%%%%% ~~ the end ~~ %%%%%

