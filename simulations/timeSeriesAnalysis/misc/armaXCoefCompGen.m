function [H,pValue,errFlag] = armaXCoefCompGen(fitResults1,fitResults2,...
    compOpt,significance)
%ARMAXCOEFCOMP tests whether the ARMA coefficients of 2 models are different
%
%SYNOPSIS [H,pValue,errFlag] = armaXCoefCompGen(fitResults1,fitResults2,...
%    compOpt,significance)
%
%INPUT  
%   Mandatory
%       fitResults1 : Structure output from armaXFitKalman. Must contain at
%                     least the following fields:
%           .arParamK     : AR coefficients (row vector).
%           .maParamK     : MA coefficients (row vector).
%           .xParamK      : X coefficients (row vector).
%           .varCovMatF   : Variance-covariance matrix of coefficients.
%   Optional
%       armaXCoef2  : Same as armaXCoef1, but for 2nd model. If it's a vector
%                     of zeros (i.e. if comparing 1st model to zero), it
%                     can be entered as []. Default: vector of zeros
%       varCovMat2  : Same as varCovMat1, but for 2nd model. 
%                     Default: Diagonal Matrix with 10^-10 on the diagonals.
%       numDegFree2 : Number of degrees of freedom in estimating armaXCoef2.
%                     Defaults: numDegFree1.
%       compOpt     : Type of comparison to be performed:
%                     -'global': Ensemble comparison of all coefficients, 
%                      taking into account coefficient covariances.
%                     -'element': Separate comparison for each coefficient,
%                      ignoring covariances.
%                     -'AR/MA/X': Separate comparison for AR coefficients,
%                       MA cofficients and X coefficients, taking into 
%                       account covariances within each set of coefficients
%                       but ignoring covariances between sets.
%                     Default: 'global'.
%       significance: Significance level of hypothesis test. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       pValue  : probability that difference between the two sets of
%                 coefficients is >= difference observed assuming that the
%                 null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = [];
pValue = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin < 3 
    disp('--armaXCoefCompGen: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%check 1st model
if ~isempty(armaXCoef1.arParam)
    [nRow,arOrder1] = size(armaXCoef1.arParam);
    if nRow ~= 1
        disp('--armaXCoefCompGen: armaXCoef1.arParam should be a row vector!');
        errFlag = 1;
    end
else
    arOrder1 = 0;
end
if ~isempty(armaXCoef1.maParam)
    [nRow,maOrder1] = size(armaXCoef1.maParam);
    if nRow ~= 1
        disp('--armaXCoefCompGen: armaXCoef1.maParam should be a row vector!');
        errFlag = 1;
    end
else
    maOrder1 = 0;
end
if ~isempty(armaXCoef1.xParam)
    [nRow,xOrder1] = size(armaXCoef1.xParam);
    if nRow ~= 1
        disp('--armaXCoefCompGen: armaXCoef1.xParam should be a row vector!');
        errFlag = 1;
    end
    xOrder1 = xOrder1 - 1;
else
    xOrder1 = -1;
end
numParam1 = arOrder1 + maOrder1 + xOrder1 + 1;
if numParam1 == 0 %exit if there are no parameters to test
    disp('--armaXCoefCompGen: Input for armaXCoef1 not valid!');
    errFlag = 1;
    return
end

%check its variance-covariance matrix
[nRow,nCol] = size(varCovMat1);
if nRow ~= nCol || nRow ~= numParam1
    disp('--armaXCoefCompGen: varCovMat1 should be a square matrix of side length equal to AR order + MA order + X order +1!');
    errFlag = 1;
end

%make sure that the number of degrees of freedom in estimation is positive
if numDegFree1 <= 0
    disp('--armaXCoefCompGen: numDegFree1 should be a positive integer!');
    errFlag = 1;
end
    
%exit if there are problems in mandatory input
if errFlag
    disp('--armaXCoefCompGen: Please fix input data!');
    return
end

%assign default values of optional arguments
armaXCoef2_def = struct('arParam',zeros(1,arOrder1),'maParam',zeros(1,maOrder1),'xOrder',zeros(1,xOrder1+1));
varCovMat2_def = 1e-10*eye(numParam1);;
numDegFree2_def = numDegFree1;
compOpt_def = 'global';
significance_def = 0.05;

if nargin < 4 || isempty(armaXCoef2) %if 2nd model was not input

    %assign 2nd model to zero (same order as 1st model)
    arOrder2 = arOrder1;
    maOrder2 = maOrder1;
    xOrder2 = xOrder1;
    numParam2 = numParam1;
    armaXCoef2 = armaXCoef2_def;

else %if user specified a 2nd model
    
    %check 2nd model
    if ~isempty(armaXCoef2.arParam)
        [nRow,arOrder2] = size(armaXCoef2.arParam);
        if nRow ~= 1
            disp('--armaXCoefCompGen: armaXCoef2.arParam should be a row vector!');
            errFlag = 1;
        end
    else
        arOrder2 = 0;
    end
    if ~isempty(armaXCoef2.maParam)
        [nRow,maOrder2] = size(armaXCoef2.maParam);
        if nRow ~= 1
            disp('--armaXCoefCompGen: armaXCoef2.maParam should be a row vector!');
            errFlag = 1;
        end
    else
        maOrder2 = 0;
    end
    if ~isempty(armaXCoef2.xParam)
        [nRow,xOrder2] = size(armaXCoef2.xParam);
        if nRow ~= 1
            disp('--armaXCoefCompGen: armaXCoef2.xParam should be a row vector!');
            errFlag = 1;
        end
        xOrder2 = xOrder2 - 1;
    else
        xOrder2 = -1;
    end
    numParam2 = arOrder2 + maOrder2 + xOrder2 + 1;

end %(if nargin < 4 || isempty(armaXCoef2) ... else ...)

if nargin < 5 || isempty(varCovMat2) %if var-cov matrix of 2nd model was not input

    varCovMat2 = varCovMat2_def;

else %if user specified a var-cov matrix for 2nd model

    [nRow,nCol] = size(varCovMat2);
    if nRow ~= nCol || nRow ~= numParam2
        disp('--armaXCoefCompGen: varCovMat2 should be a square matrix of side length equal to AR order + MA order + X order +1!');
        errFlag = 1;
    end
    
end %(if nargin < 5 || isempty(varCovMat2) ... else ...)

%check number of degrees of freedom used in estimating the second model
if nargin < 6 || isempty(numDegFree2)
    numDegFree2 = numDegFree2_def;
elseif numDegFree2 <= 0
    disp('--armaXCoefCompGen: numDegFree2 should be a positive integer!');
    errFlag = 1;
end %(if nargin < 6 || isempty(numDegFree2) ... elseif ...)

if nargin  < 7 || isempty(compOpt) %if compOpt was not input

    compOpt = compOpt_def; %assign default

else %if user specified compOpt

    %check that it's a valid option
    if ~strcmp(compOpt,'global') && ~strcmp(compOpt,'element') && ...
            ~strcmp(compOpt,'AR/MA/X')
        disp('--armaXCoefCompGen: "compOpt" should be either "global", "element" or "AR/MA/X"!');
        errFlag = 1;
    end

end %(if nargin < 7 || isempty(compOpt) ... else ...)

if nargin < 8 || isempty(significance) %if user did not specify significance level of test
    
    significance = significance_def; %assign default

else %if user specified significance level
    
    %check that it's in the correct range
    if significance < 0 || significance > 1
        disp('--armaXCoefCompGen: "significance level should be between 0 and 1!');
        errFlag = 1;
    end

end %(if nargin < 8 || isempty(significance) ... else ...)

%exit if there are problems in input data
if errFlag
    disp('--armaXCoefCompGen: Please fix input data!');
    return
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
armaXParam1 = [armaXCoef1.arParam zeros(1,max(0,arOrderDiff)) ...
    armaXCoef1.maParam zeros(1,max(0,maOrderDiff)) ...
    armaXCoef1.xParam zeros(1,max(0,xOrderDiff))];

%fill 2nd parameter vector
armaXParam2 = [armaXCoef2.arParam zeros(1,max(0,-arOrderDiff)) ...
    armaXCoef2.maParam zeros(1,max(0,-maOrderDiff)) ...
    armaXCoef2.xParam zeros(1,max(0,-xOrderDiff))];

%fill 1st variance-covariance matrix
tmp = [varCovMat1(:,1:arOrder1) zeros(numParam1,max(0,arOrderDiff)) ...
    varCovMat1(:,arOrder1+1:arOrder1+maOrder1) zeros(numParam1,max(0,maOrderDiff)) ...
    varCovMat1(:,arOrder1+maOrder1+1:end) zeros(numParam1,max(0,xOrderDiff))];
varCovMatT1 = [tmp(1:arOrder1,:); zeros(max(0,arOrderDiff),numParamMax); ...
    tmp(arOrder1+1:arOrder1+maOrder1,:); zeros(max(0,maOrderDiff),numParamMax); ...
    tmp(arOrder1+maOrder1+1:end,:); zeros(max(0,xOrderDiff),numParamMax)];

%fill 2nd variance-covariance matrix
tmp = [varCovMat2(:,1:arOrder2) zeros(numParam2,max(0,-arOrderDiff)) ...
    varCovMat2(:,arOrder2+1:arOrder2+maOrder2) zeros(numParam2,max(0,-maOrderDiff)) ...
    varCovMat2(:,arOrder2+maOrder2+1:end) zeros(numParam2,max(0,-xOrderDiff))];
varCovMatT2 = [tmp(1:arOrder2,:); zeros(max(0,-arOrderDiff),numParamMax); ...
    tmp(arOrder2+1:arOrder2+maOrder2,:); zeros(max(0,-maOrderDiff),numParamMax); ...
    tmp(arOrder2+maOrder2+1:end,:); zeros(max(0,-xOrderDiff),numParamMax)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hypothesis testing
%Null hypothesis: The two models are the same, i.e. their coefficients are identical
%Alternative hypothesis: Model coefficients are not identical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the number of points in the smaller of the two data sets
numDegFreeMin = min(numDegFree1,numDegFree2);

switch compOpt

    case 'global' %compare all coefficients at the same time

        %calculate vector of differences in coefficients
        diffM = armaXParam1 - armaXParam2;

        %calculate variance-covariance matrix of difference vector
        diffV = varCovMatT1 + varCovMatT2;

        %compute testStatistic
        testStatistic = diffM*(diffV\diffM')/numParamMax;

        %get the p-value of the test statistic assuming an F-distribution
        pValue = 1 - fcdf(testStatistic,numParamMax,numDegFreeMin);

        %compare p-value to significance level
        if pValue < significance %if p-value is smaller than probability of type I error
            H = 1; %reject null hypothesis that the two models are identical
        else %if p-value is larger than probability of type I error
            H = 0; %cannot reject null hypothesis
        end

    case 'element' %compare one coefficient at a time

        %calculate vector of differences in coefficients
        diffM = armaXParam1 - armaXParam2;

        %calculate variance-covariance matrix of difference vector
        diffV = diag(varCovMatT1 + varCovMatT2);

        %calculate the test statistic
        testStatistic = diffM./sqrt(diffV);

        %get the p-value assuming that each elemenet in testStatistic
        %follows a student t-distribution
        pValue = 1 - tcdf(abs(testStatistic),numDegFreeMin);

        %compare p-value to significance level
        for i=1:length(pValue)
            if pValue(i) < significance/2 %if p-value is smaller than probability of type I error
                H(i) = 1; %reject hypothesis that elements are identical
            else %if p-value is larger than probability of type I error
                H(i) = 0; %cannot reject hypothesis
            end
        end

    case 'AR/MA/X' %compare AR, MA and X coefficients alone
     
        H = NaN*ones(1,3);
        pValue = H;
        
        %calculate vector of differences in coefficients
        diffM = armaXParam1 - armaXParam2;

        %AR test
        if arOrderMax ~= 0

            %calculate variance-covariance matrix of AR difference vector
            diffV = varCovMatT1(1:arOrderMax,1:arOrderMax) + ...
                varCovMatT2(1:arOrderMax,1:arOrderMax);

            %compute testStatistic
            testStatistic = diffM(1:arOrderMax)*(diffV\diffM(1:arOrderMax)')/arOrderMax;

            %get the p-value of the test statistic assuming an F-distribution
            pValue(1) = 1 - fcdf(testStatistic,arOrderMax,numDegFreeMin);

            %compare p-value to significance level
            if pValue(1) < significance %if p-value is smaller than probability of type I error
                H(1) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(1) = 0; %cannot reject null hypothesis
            end

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
            pValue(2) = 1 - fcdf(testStatistic,maOrderMax,numDegFreeMin);

            %compare p-value to significance level
            if pValue(2) < significance %if p-value is smaller than probability of type I error
                H(2) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(2) = 0; %cannot reject null hypothesis
            end

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
            pValue(3) = 1 - fcdf(testStatistic,xOrderMax+1,numDegFreeMin);

            %compare p-value to significance level
            if pValue(3) < significance %if p-value is smaller than probability of type I error
                H(3) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(3) = 0; %cannot reject null hypothesis
            end

        end

end


%%%%% ~~ the end ~~ %%%%%

