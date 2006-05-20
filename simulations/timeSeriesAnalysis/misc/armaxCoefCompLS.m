function [H,pValue,errFlag] = armaxCoefCompLS(armaxCoef1,armaxCoef2,...
    varCovMat1,varCovMat2,compOpt,significance)
%ARMAXCOEFCOMP tests whether the ARMA coefficients of 2 models are different
%
%SYNOPSIS [H,pValue,errFlag] = armaxCoefCompLS(armaxCoef1,armaxCoef2,...
%    varCovMat1,varCovMat2,compOpt,significance)
%
%INPUT  
%   Mandatory
%       armaxCoef1  : Structure containing ARMAX coefficients of 1st model:
%           .arParam      : AR coefficients (row vector).
%           .maParam      : MA coefficients (row vector).
%           .xParam       : X coefficients (row vector).
%   Optional
%       armaxCoef2  : Same as armaxCoef1, but for 2nd model. If it's a vector
%                     of zeros (i.e. if comparing 1st model to zero), it
%                     can be entered as []. Default: vector of zeros
%       varCovMat1  : Variance-covariance matrix corresponding to 1st model:
%           .cofactorMat  : Cofactor matrix for largest of AR, MA and X 
%                           orders among the 2 models.
%           .posterioriVar: A posteriori estimate of residuals' variance 
%                           for fitting data with 1st model.
%                     Default: Diagonal Matrix with 10^-10 on the diagonals.
%       varCovMat2  : Same as varCovMat1, but for 2nd model. 
%                     Default: Diagonal Matrix with 10^-10 on the diagonals.
%       compOpt     : Type of comparison to be performed:
%                     -'global': Simultaneous comparison of all coefficients, 
%                      taking into account coefficients' covariances.
%                     -'element': Separate comparison for each coefficient,
%                      ignoring their covariances.
%                     -'AR/MA/X': Separate comparison for AR coefficients,
%                       MA cofficients and X coefficients, taking into 
%                       account covariances within each set of coefficients
%                       but ignoring covariances between the sets.
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
if nargin < 1 || isempty(armaxCoef1) %if 1st model was not input
    disp('--armaxCoefCompLS: You must input at least 1 set of ARMA coefficients!');
    errFlag = 1;
    return
end

%assign default values of some optional arguments
compOpt_def = 'global';
significance_def = 0.05;

%check 1st model
if ~isempty(armaxCoef1.arParam)
    [nRow,arOrder1] = size(armaxCoef1.arParam);
    if nRow ~= 1
        disp('--armaxCoefCompLS: armaxCoef1.arParam should be a row vector!');
        errFlag = 1;
    end
else
    arOrder1 = 0;
end
if ~isempty(armaxCoef1.maParam)
    [nRow,maOrder1] = size(armaxCoef1.maParam);
    if nRow ~= 1
        disp('--armaxCoefCompLS: armaxCoef1.maParam should be a row vector!');
        errFlag = 1;
    end
else
    maOrder1 = 0;
end
if ~isempty(armaxCoef1.xParam)
    [nRow,xOrder1] = size(armaxCoef1.xParam);
    if nRow ~= 1
        disp('--armaxCoefCompLS: armaxCoef1.xParam should be a row vector!');
        errFlag = 1;
    end
    xOrder1 = xOrder1 - 1;
else
    xOrder1 = -1;
end
if arOrder1 == 0 && maOrder1 == 0 && xOrder1 == 0 %exit if no model is of order 0
    disp('--armaxCoefCompLS: Input for armaxCoef1 not valid!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(armaxCoef2) %if 2nd model was not input

    %assign 2nd model to zero (same order as 1st model)
    arOrder2 = arOrder1;
    maOrder2 = maOrder1;
    xOrder2 = xOrder1;
    armaxCoef2.arParam = zeros(1,arOrder1);
    armaxCoef2.maParam = zeros(1,maOrder1);
    armaxCoef2.xParam = zeros(1,xOrder1+1);

else %if user specified a 2nd model
    
    %check 2nd model
    if ~isempty(armaxCoef2.arParam)
        [nRow,arOrder2] = size(armaxCoef2.arParam);
        if nRow ~= 1
            disp('--armaxCoefCompLS: armaxCoef2.arParam should be a row vector!');
            errFlag = 1;
        end
    else
        arOrder2 = 0;
    end
    if ~isempty(armaxCoef2.maParam)
        [nRow,maOrder2] = size(armaxCoef2.maParam);
        if nRow ~= 1
            disp('--armaxCoefCompLS: armaxCoef2.maParam should be a row vector!');
            errFlag = 1;
        end
    else
        maOrder2 = 0;
    end
    if ~isempty(armaxCoef2.xParam)
        [nRow,xOrder2] = size(armaxCoef2.xParam);
        if nRow ~= 1
            disp('--armaxCoefCompLS: armaxCoef2.xParam should be a row vector!');
            errFlag = 1;
        end
        xOrder2 = xOrder2 - 1;
    else
        xOrder2 = -1;
    end

end %(if nargin < 2 || isempty(armaxCoef2) ... else ...)

%get largest orders
arOrderL = max(arOrder1,arOrder2);
maOrderL = max(maOrder1,maOrder2);
xOrderL = max(xOrder1,xOrder2);
combOrder = arOrderL + maOrderL + xOrderL + 1;

if nargin < 3 || isempty(varCovMat1) %if var-cov matrix of 1st model was not input

    varCovMat1.cofactorMat = 1e-10*eye(combOrder);
    varCovMat1.posterioriVar = 1;

else %if user specified a var-cov matrix for 1st model

    [nRow,nCol] = size(varCovMat1.cofactorMat);
    if nRow ~= nCol || nRow ~= combOrder
        disp('--armaxCoefCompLS: varCovMat1.cofactorMat should be a square matrix of side length equal to largest AR order + largest MA order!');
        errFlag = 1;
    end
    if varCovMat1.posterioriVar <= 0
        disp('--armaxCoefCompLS: varCovMat1.posterioriVar should be positive!');
        errFlag = 1;
    end

end %(if nargin < 3 || isempty(varCovMat1) ... else ...)

if nargin < 4 || isempty(varCovMat2) %if var-cov matrix of 2nd model was not input

    varCovMat2.cofactorMat = 1e-10*eye(combOrder);
    varCovMat2.posterioriVar = 1;

else %if user specified a var-cov matrix for 2nd model

    [nRow,nCol] = size(varCovMat2.cofactorMat);
    if nRow ~= nCol || nRow ~= combOrder
        disp('--armaxCoefCompLS: varCovMat2.cofactorMat should be a square matrix of side length equal to largest AR order + largest MA order!');
        errFlag = 1;
    end
    if varCovMat2.posterioriVar <= 0
        disp('--armaxCoefCompLS: varCovMat2.posterioriVar should be positive!');
        errFlag = 1;
    end

end %(if nargin < 4 || isempty(varCovMat2) ... else ...)

if nargin  < 5 || isempty(compOpt) %if compOpt was not input

    compOpt = compOpt_def; %assign default

else %if user specified compOpt

    %check that it's a valid option
    if ~strcmp(compOpt,'global') && ~strcmp(compOpt,'element') && ...
            ~strcmp(compOpt,'AR/MA/X')
        disp('--armaxCoefCompLS: "compOpt" should be either "global", "element" or "AR/MA/X"!');
        errFlag = 1;
    end

end %(if nargin  < 5 || isempty(compOpt) ... else ...)

if nargin < 6 || isempty(significance) %if user did not specify significance level of test
    
    significance = significance_def; %assign default

else %if user specified significance level
    
    %check that it's in the correct range
    if significance < 0 || significance > 1
        disp('--armaxCoefCompLS: "significance level should be between 0 and 1!');
        errFlag = 1;
    end

end %(if nargin < 6 || isempty(significance) ... else ...)

%exit if there are problems in input data
if errFlag
    disp('--armaxCoefCompLS: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model order matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize constraint matrices
constMat1 = [];
constMat2 = [];

%compare AR orders, add zeros to AR coefficients and devise constraint 
%matrices to update variance-covariance matrix
arOrderDiff = arOrder2 - arOrder1;
if arOrderDiff < 0
    armaxCoef2.arParam = [armaxCoef2.arParam zeros(1,-arOrderDiff)];
    constMat2 = [zeros(-arOrderDiff,arOrder2) eye(-arOrderDiff) ...
        zeros(-arOrderDiff,maOrderL+xOrderL+1)];
elseif arOrderDiff > 0
    armaxCoef1.arParam = [armaxCoef1.arParam zeros(1,arOrderDiff)];
    constMat1 = [zeros(arOrderDiff,arOrder1) eye(arOrderDiff) ...
        zeros(arOrderDiff,maOrderL+xOrderL+1)];
end

%compare MA orders, add zeros to MA coefficients and devise constraint 
%matrices to update variance-covariance matrix
maOrderDiff = maOrder2 - maOrder1;
if maOrderDiff < 0
    armaxCoef2.maParam = [armaxCoef2.maParam zeros(1,-maOrderDiff)];
    constMat2 = [constMat2; zeros(-maOrderDiff,arOrderL) ...
        zeros(-maOrderDiff,maOrder2) eye(-maOrderDiff) ...
        zeros(-maOrderDiff,xOrderL)];
elseif maOrderDiff > 0
    armaxCoef1.maParam = [armaxCoef1.maParam zeros(1,maOrderDiff)];
    constMat1 = [constMat1; zeros(maOrderDiff,arOrderL) ...
        zeros(maOrderDiff,maOrder1) eye(maOrderDiff) ...
        zeros(maOrderDiff,xOrderL)];
end

%compare X orders, add zeros to X coefficients and devise constraint
%matrices to update variance-covariance matrix
xOrderDiff = xOrder2 - xOrder1;
if xOrderDiff < 0
    armaxCoef2.xParam = [armaxCoef2.xParam zeros(1,-xOrderDiff)];
    constMat2 = [constMat2; zeros(-xOrderDiff,arOrderL+maOrderL) ...
        zeros(-xOrderDiff,xOrder2+1) eye(-xOrderDiff)];
elseif xOrderDiff > 0
    armaxCoef1.xParam = [armaxCoef1.xParam zeros(1,xOrderDiff)];
    constMat1 = [constMat1; zeros(xOrderDiff,arOrderL+maOrderL) ...
        zeros(xOrderDiff,xOrder1+1) eye(xOrderDiff)];
end

%combine AR, MA and X coefficients for each model
armaxParam1 = [armaxCoef1.arParam armaxCoef1.maParam armaxCoef1.xParam];
armaxParam2 = [armaxCoef2.arParam armaxCoef2.maParam armaxCoef2.xParam];

%get new variance-covariance matrices, if needed
if ~isempty(constMat1)
    cofactorMat = varCovMat1.cofactorMat;
    varCovMatT1 = (cofactorMat - cofactorMat*constMat1'*(inv(constMat1*...
        cofactorMat*constMat1'))*constMat1*cofactorMat)...
        *varCovMat1.posterioriVar;
else
    varCovMatT1 = varCovMat1.cofactorMat*varCovMat1.posterioriVar;
end
if ~isempty(constMat2)
    cofactorMat = varCovMat2.cofactorMat;
    varCovMatT2 = (cofactorMat - cofactorMat*constMat2'*(inv(constMat2*...
        cofactorMat*constMat2'))*constMat2*cofactorMat)...
        *varCovMat2.posterioriVar;
else
    varCovMatT2 = varCovMat2.cofactorMat*varCovMat2.posterioriVar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hypothesis testing
%Null hypothesis: The two models are the same, i.e. their coefficients are identical
%Alternative hypothesis: Model coefficients are not identical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch compOpt

    case 'global' %compare all coefficients at the same time

        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %calculate variance-covariance matrix of difference vector
        diffV = varCovMatT1 + varCovMatT2;
        %         diffV = [eye(combOrder) -eye(combOrder)]*[varCovMatT1 zeros(combOrder); ...
        %             zeros(combOrder) varCovMatT2]*[eye(combOrder) -eye(combOrder)]';
        %the above is the analytical simplification of this formula!!!

        %compute testStatistic
        testStatistic = diffM*(diffV\diffM')/combOrder;

        %get the p-value of the test statistic assuming a Fisher distribution
        pValue = 1 - fcdf(testStatistic,combOrder,2000);

        %compare p-value to significance
        if pValue < significance %if p-value is smaller than probability of type I error
            H = 1; %reject null hypothesis that the two models are identical
        else %if p-value is larger than probability of type I error
            H = 0; %cannot reject null hypothesis
        end

    case 'element' %compare one coefficient at a time

        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %calculate variance-covariance matrix of difference vector
        varVec1 = diag(varCovMatT1)';
        varVec2 = diag(varCovMatT2)';
        diffV = varVec1 + varVec2;

        %calculate the test statistic
        testStatistic = diffM./sqrt(diffV);

        %get the p-value assuming that each elemenet in testStatistic
        %is normally distributed with mean zero and variance 1
        pValue = 1 - normcdf(abs(testStatistic),0,1);

        %compare p-value to significance
        for i=1:length(pValue)
            if pValue(i) < significance/2 %if p-value is smaller than probability of type I error
                %             if pValue(i) < significance %if p-value is smaller than probability of type I error
                H(i) = 1; %reject hypothesis that elements are identical
            else %if p-value is larger than probability of type I error
                H(i) = 0; %cannot reject hypothesis
            end
        end

    case 'AR/MA/X' %compare AR, MA and X coefficients alone
     
        H = NaN*ones(1,3);
        pValue = H;
        
        %calculate vector of differences in coefficients
        diffM = armaxParam1 - armaxParam2;

        %AR test
        if arOrderL ~= 0

            %calculate variance-covariance matrix of AR difference vector
            diffV = varCovMatT1(1:arOrderL,1:arOrderL) + ...
                varCovMatT2(1:arOrderL,1:arOrderL);

            %compute testStatistic
            testStatistic = diffM(1:arOrderL)*(diffV\diffM(1:arOrderL)')/arOrderL;

            %get the p-value of the test statistic assuming a Fisher distribution
            pValue(1) = 1 - fcdf(testStatistic,arOrderL,2000);

            %compare p-value to significance
            if pValue(1) < significance %if p-value is smaller than probability of type I error
                H(1) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(1) = 0; %cannot reject null hypothesis
            end

        end

        %MA test
        if maOrderL ~= 0

            %calculate variance-covariance matrix of MA difference vector
            diffV = varCovMatT1(arOrderL+1:arOrderL+maOrderL,...
                arOrderL+1:arOrderL+maOrderL) + ...
                varCovMatT2(arOrderL+1:arOrderL+maOrderL,...
                arOrderL+1:arOrderL+maOrderL);

            %compute testStatistic
            testStatistic = diffM(arOrderL+1:arOrderL+maOrderL)*...
                (diffV\diffM(arOrderL+1:arOrderL+maOrderL)')/maOrderL;

            %get the p-value of the test statistic assuming a Fisher distribution
            pValue(2) = 1 - fcdf(testStatistic,maOrderL,2000);

            %compare p-value to significance
            if pValue(2) < significance %if p-value is smaller than probability of type I error
                H(2) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(2) = 0; %cannot reject null hypothesis
            end

        end

        %X test
        if xOrderL ~= -1

            %calculate variance-covariance matrix of X difference vector
            diffV = varCovMatT1(arOrderL+maOrderL+1:end,...
                arOrderL+maOrderL+1:end) + ...
                varCovMatT2(arOrderL+maOrderL+1:end,...
                arOrderL+maOrderL+1:end);

            %compute testStatistic
            testStatistic = diffM(arOrderL+maOrderL+1:end)*...
                (diffV\diffM(arOrderL+maOrderL+1:end)')/(xOrderL+1);

            %get the p-value of the test statistic assuming a Fisher distribution
            pValue(3) = 1 - fcdf(testStatistic,xOrderL+1,2000);

            %compare p-value to significance
            if pValue(3) < significance %if p-value is smaller than probability of type I error
                H(3) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(3) = 0; %cannot reject null hypothesis
            end

        end

end


%%%%% ~~ the end ~~ %%%%%

