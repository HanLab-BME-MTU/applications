function [H,pValue,errFlag] = armaCoefComp(armaCoef1,armaCoef2,...
    varCovMat1,varCovMat2,compOpt,significance)
%ARMACOEFCOMP tests whether the ARMA coefficients of 2 models are different
%
%SYNOPSIS [H,pValue,errFlag] = armaCoefComp(armaCoef1,armaCoef2,...
%    varCovMat1,varCovMat2,compOpt,significance)
%
%INPUT  armaCoef1   : Structure containing ARMA coefficients of 1st model:
%           .arParam      : AR coefficients (row vector).
%           .maParam      : MA coefficients (row vector).
%       armaCoef2   : Same as armaCoef1, but for 2nd model. Optional. If 
%                     comparing 1st model to zero, can be skipped or entered 
%                     as [] if other input parameters follow. Default:
%                     Vector of zeroes.
%       varCovMat1  : Variance-covariance matrix corresponding to 1st model:
%           .cofactorMat  : Cofactor matrix for largest of AR and MA orders 
%                           among the 2 models.
%           .posterioriVar: A posteriori estimate of residuals' variance 
%                           for fitting data with 1st model.
%                     It will be converted to the appropriate variance-
%                     covariance matrix given the orders of the 2 models of 
%                     interest, using least squares fitting with constraints.
%                     If it is identically zero, [] can be entered.
%       varCovMat2  : Same as varCovMat1, but for 2nd model. Optional. If 
%                     comparing 1st model to zero, can be skipped. Default:
%                     Diagonal Matrix with 10^-10 on the diagonals.
%       compOpt     : Type of comparison to be performed:
%                     -'global': Simultaneous comparison of all coefficients, 
%                      taking into account coefficients' covariances.
%                     -'element': Separate comparison for each coefficient,
%                      ignoring their covariances.
%                     -'AR/MA': Separate comparison for AR coefficients and
%                       MA cofficients, taking into account covariances
%                       between AR coefficients and between MA coefficients, 
%                       but ignoring covariances between then two sets.
%                      Optional. Default: 'global'.
%       significance: Significance level of hypothesis test. Optional. Default: 0.05.
%
%OUTPUT H       : 1 if hypothesis can be rejected, 0 otherwise.
%       pValue  : probability that difference between the two sets of
%                 coefficients is >= difference observed assuming that the
%                 null hypothesis is true.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, September 2004

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
if nargin < 1 || isempty(armaCoef1) %if 1st model was not input

    disp('--armaCoefComp: You must input at least 1 set of ARMA coefficients!');
    errFlag = 1;
    return

else %if 1st model was input
    
    %check 1st model
    if ~isempty(armaCoef1.arParam)
        [nRow,arOrder1] = size(armaCoef1.arParam);
        if nRow ~= 1
            disp('--armaCoefComp: armaCoef1.arParam should be a row vector!');
            errFlag = 1;
        end
    else
        arOrder1 = 0;
    end
    if ~isempty(armaCoef1.maParam)
        [nRow,maOrder1] = size(armaCoef1.maParam);
        if nRow ~= 1
            disp('--armaCoefComp: armaCoef1.maParam should be a row vector!');
            errFlag = 1;
        end
    else
        maOrder1 = 0;
    end
    if arOrder1 == 0 && maOrder1 == 0 %exit if no model is of order 0
        disp('--armaCoefComp: Input for armaCoef1 not valid!');
        errFlag = 1;
        return
    end

end %(if nargin < 1 || isempty(armaCoef1))

if nargin < 2 || isempty(armaCoef2) %if 2nd model was not input

    %assign 2nd model to zero (same order as 1st model)
    arOrder2 = arOrder1;
    maOrder2 = maOrder1;
    armaCoef2.arParam = zeros(1,arOrder1);
    armaCoef2.maParam = zeros(1,maOrder1);

else %if user specified a 2nd model
    
    %check 2nd model
    if ~isempty(armaCoef2.arParam)
        [nRow,arOrder2] = size(armaCoef2.arParam);
        if nRow ~= 1
            disp('--armaCoefComp: armaCoef2.arParam should be a row vector!');
            errFlag = 1;
        end
    else
        arOrder2 = 0;
    end
    if ~isempty(armaCoef2.maParam)
        [nRow,maOrder2] = size(armaCoef2.maParam);
        if nRow ~= 1
            disp('--armaCoefComp: armaCoef2.maParam should be a row vector!');
            errFlag = 1;
        end
    else
        maOrder2 = 0;
    end

end %(if nargin < 2 || isempty(armaCoef2))

%get largest orders
arOrderL = max(arOrder1,arOrder2);
maOrderL = max(maOrder1,maOrder2);
combOrder = arOrderL + maOrderL;

if nargin < 3 || isempty(varCovMat1) %if var-cov matrix of 1st model was not input

    varCovMat1.cofactorMat = 1e-10*eye(combOrder);
    varCovMat1.posterioriVar = 1;

else %if user specified a var-cov matrix for 1st model

    [nRow,nCol] = size(varCovMat1.cofactorMat);
    if nRow ~= nCol || nRow ~= combOrder
        disp('--armaCoefComp: varCovMat1.cofactorMat should be a square matrix of side length equal to largest AR order + largest MA order!');
        errFlag = 1;
    end
    if varCovMat1.posterioriVar <= 0
        disp('--armaCoefComp: varCovMat1.posterioriVar should be positive!');
        errFlag = 1;
    end

end %(if nargin < 3 || isempty(varCovMat1))

if nargin < 4 || isempty(varCovMat2) %if var-cov matrix of 2nd model was not input

    varCovMat2.cofactorMat = 1e-10*eye(combOrder);
    varCovMat2.posterioriVar = 1;

else %if user specified a var-cov matrix for 2nd model

    [nRow,nCol] = size(varCovMat2.cofactorMat);
    if nRow ~= nCol || nRow ~= combOrder
        disp('--armaCoefComp: varCovMat2.cofactorMat should be a square matrix of side length equal to largest AR order + largest MA order!');
        errFlag = 1;
    end
    if varCovMat2.posterioriVar <= 0
        disp('--armaCoefComp: varCovMat2.posterioriVar should be positive!');
        errFlag = 1;
    end

end %(if nargin < 4 || isempty(varCovMat2))

if nargin  < 5 || isempty(compOpt) %if compOpt was not input

    compOpt = 'global'; %assign default

else %if user specified compOpt

    %check that it's a valid option
    if ~strcmp(compOpt,'global') && ~strcmp(compOpt,'element') && ...
            ~strcmp(compOpt,'AR/MA')
        disp('--armaCoefComp: "compOpt" should be either "global", "element" or "AR/MA"!');
        errFlag = 1;
    end

end %(if nargin  < 5 || isempty(compOpt))

if nargin < 6 || isempty(significance) %if user did not specify significance level of test
    
    significance = 0.05; %assign default

else %if user specified significance level
    
    %check that it's in the correct range
    if significance < 0 || significance > 1
        disp('--armaCoefComp: "significance level should be between 0 and 1!');
        errFlag = 1;
    end

end %(if nargin < 6 || isempty(significance))

%exit if there are problems in input data
if errFlag
    disp('--armaCoefComp: Please fix input data!');
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
    armaCoef2.arParam = [armaCoef2.arParam zeros(1,-arOrderDiff)];
    constMat2 = [zeros(-arOrderDiff,arOrder2) eye(-arOrderDiff) ...
        zeros(-arOrderDiff,maOrderL)];
elseif arOrderDiff > 0
    armaCoef1.arParam = [armaCoef1.arParam zeros(1,arOrderDiff)];
    constMat1 = [zeros(arOrderDiff,arOrder1) eye(arOrderDiff) ...
        zeros(arOrderDiff,maOrderL)];
end

%compare MA orders, add zeros to MA coefficients and devise constraint 
%matrices to update variance-covariance matrix
maOrderDiff = maOrder2 - maOrder1;
if maOrderDiff < 0
    armaCoef2.maParam = [armaCoef2.maParam zeros(1,-maOrderDiff)];
    constMat2 = [constMat2; zeros(-maOrderDiff,arOrderL) ...
        zeros(-maOrderDiff,maOrder2) eye(-maOrderDiff)];
elseif maOrderDiff > 0
    armaCoef1.maParam = [armaCoef1.maParam zeros(1,maOrderDiff)];
    constMat1 = [constMat1; zeros(maOrderDiff,arOrderL) ...
        zeros(maOrderDiff,maOrder1) eye(maOrderDiff)];
end

%combine AR and MA coefficients for each model
armaParam1 = [armaCoef1.arParam armaCoef1.maParam];
armaParam2 = [armaCoef2.arParam armaCoef2.maParam];

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
        diffM = armaParam1 - armaParam2;

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
        diffM = armaParam1 - armaParam2;

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

    case 'AR/MA' %compare AR coefficients alone and MA coefficients alone

        H = NaN*ones(1,2);
        
        %calculate vector of differences in coefficients
        diffM = armaParam1 - armaParam2;

        %AR test
        if arOrderL ~= 0

            %calculate variance-covariance matrix of AR difference vector
            diffV = varCovMatT1(1:arOrderL,1:arOrderL) + ...
                varCovMatT2(1:arOrderL,1:arOrderL);

            %compute testStatistic
            testStatistic = diffM(1:arOrderL)*(diffV\diffM(1:arOrderL)')/arOrderL;

            %get the p-value of the test statistic assuming a Fisher distribution
            pValue = 1 - fcdf(testStatistic,arOrderL,2000);

            %compare p-value to significance
            if pValue(1) < significance %if p-value is smaller than probability of type I error
                H(1) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(1) = 0; %cannot reject null hypothesis
            end

        else
            pValue(1) = NaN;
        end

        %MA test
        if maOrderL ~= 0

            %calculate variance-covariance matrix of MA difference vector
            diffV = varCovMatT1(1+arOrderL:end,1+arOrderL:end) + ...
                varCovMatT2(1+arOrderL:end,1+arOrderL:end);

            %compute testStatistic
            testStatistic = diffM(arOrderL+1:end)*(diffV\diffM(...
                arOrderL+1:end)')/maOrderL;

            %get the p-value of the test statistic assuming a Fisher distribution
            pValue(2) = 1 - fcdf(testStatistic,maOrderL,2000);

            %compare p-value to significance
            if pValue(2) < significance %if p-value is smaller than probability of type I error
                H(2) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(2) = 0; %cannot reject null hypothesis
            end

        else
            pValue(2) = NaN;
        end

end


%%%%% ~~ the end ~~ %%%%%

