function [H,errFlag] = armaCoefComp(armaCoef1,varCovMat1,armaCoef2,...
    varCovMat2,compOpt,significance)
%ARMACOEFCOMP tests whether the the ARMA coefficients of 2 models are identical
%
%SYNOPSIS [H,errFlag] = armaCoefComp(armaCoef1,varCovMat1,armaCoef2,...
%    varCovMat2,compOpt,significance)
%
%INPUT  armaCoef1   : Structure containing ARMA coefficients of 1st model:
%           .arParam: AR coefficients (row vector).
%           .maParam: MA coefficients (row vector).
%       varCovMat1  : Variance-covariance matrix of ARMA coeffients of 1st
%                     model. If it is identically zero, [] can be entered.
%       armaCoef2   : Same as armaCoef1, but for 2nd model. Optional. If 
%                     comparing 1st model to zero, can be skipped. Default:
%                     Vector of zeroes.
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
%       errFlag : 0 if function executes normally, 1 otherwise.
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
if nargin < 2
    disp('--armaCoefComp: You must input at least 1 set of ARMA coefficients and their variance-covariance matrix!');
    errFlag = 1;
    return
end

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
if arOrder1 == 0 && maOrder1 == 0 %exit if no model was entered
    disp('--armaCoefComp: Input for armaCoef1 not valid!');
    errFlag = 1;
    return
end
if isempty(varCovMat1) %assign default variance-covariance matrix if empty
    varCovMat1 = 1e-10*eye(arOrder1+maOrder1);
else %check matrix dimensions
    [nRow,nCol] = size(varCovMat1);
    if nRow ~= nCol || nRow ~= arOrder1 + maOrder1
        disp('--armaCoefComp: varCovMat1 should be a square matrix of side length equal to AR order + MA order!');
        errFlag = 1;
    end
end

if nargin < 4 || isempty(armaCoef2) %if 2nd model was not input

    %assign 2nd model to zero
    arOrder2 = arOrder1;
    maOrder2 = maOrder1;
    armaCoef2.arParam = zeros(1,arOrder1);
    armaCoef2.maParam = zeros(1,maOrder1);
    varCovMat2 = 1e-10*eye(arOrder1+maOrder1);

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
    
    if arOrder2 == 0 && maOrder2 == 0 %if no model was input, assume default
        arOrder2 = arOrder1;
        maOrder2 = maOrder1;
        armaCoef2.arParam = zeros(1,arOrder1);
        armaCoef2.maParam = zeros(1,maOrder1);
        varCovMat2 = 1e-10*eye(arOrder1+maOrder1);
    else

        if isempty(varCovMat2) %assign default variance-covariance matrix if empty
            varCovMat2 = 1e-10*eye(arOrder2+maOrder2);
        else %check matrix dimensions
            [nRow,nCol] = size(varCovMat2);
            if nRow ~= nCol || nRow ~= arOrder2 + maOrder2
                disp('--armaCoefComp: varCovMat2 should be a square matrix of side length equal to AR order + MA order!');
                errFlag = 1;
            end
        end
    
    end

end %(if nargin < 4)

if nargin  < 5 %if compOpt was not input

    compOpt = 'global'; %assign default

else %if user specified compOpt

    %check that it's a valid option
    if ~strcmp(compOpt,'global') && ~strcmp(compOpt,'element') && ...
            ~strcmp(compOpt,'AR/MA')
        disp('--armaCoefComp: "compOpt" should be either "global", "element" or "AR/MA"!');
        errFlag = 1;
    end

end %(if nargin  < 5)

if nargin < 6 %if user did not specify test significance
    
    significance = 0.05; %assign default

else %if user specified significance
    
    %check that it's in the correct range
    if significance < 0 || significance > 1
        disp('--armaCoefComp: "significance should be between 0 and 1!');
        errFlag = 1;
    end

end %(if nargin < 6)

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
arOrderL = arOrder1; %get largest of the two AR orders
if arOrder1 < arOrder2 %add zeroes to armaCoef1.arParam and to varCovMat1
    armaCoef1.arParam = [armaCoef1.arParam zeros(1,arOrder2-arOrder1)];
    varCovMat1 = blkdiag(varCovMat1(1:arOrder1,1:arOrder1),...
        min(abs(nonzeros(varCovMat1)))*eye(arOrder2-arOrder1),...
        varCovMat1(arOrder1+1:end,arOrder1+1:end));
    arOrderL = arOrder2;
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

%get "combined" ARMA order and largest of the two MA orders
combOrder = length(armaParam1); %=length(armaParam2)
maOrderL = combOrder - arOrderL;

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
        diffV = [eye(combOrder) -eye(combOrder)]*[varCovMat1 zeros(combOrder); ...
            zeros(combOrder) varCovMat2]*[eye(combOrder) -eye(combOrder)]';

        %compute testStatistic
        testStatistic = diffM*(diffV\diffM')/combOrder;

        %get the p-value of the test statistic assuming a chi2 distribution
        pValue = 1 - chi2cdf(testStatistic,combOrder);

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
        varVec1 = diag(varCovMat1)';
        varVec2 = diag(varCovMat2)';
        diffV = varVec1 + varVec2;

        %calculate the test statistic
        testStatistic = diffM./sqrt(diffV);
        %         testStatistic = diffM.^2./diffV;

        %get the p-value assuming that each elemenet in testStatistic
        %is normally distributed with mean zero and variance 1
        pValue = 1 - normcdf(abs(testStatistic),0,1);
        %         pValue = 1 - chi2cdf(testStatistic,1);

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
            diffV = [eye(arOrderL) -eye(arOrderL)]*[varCovMat1(1:arOrderL,...
                1:arOrderL) zeros(arOrderL); zeros(arOrderL) varCovMat2(...
                1:arOrderL,1:arOrderL)]*[eye(arOrderL) -eye(arOrderL)]';

            %compute testStatistic
            testStatistic = diffM(1:arOrderL)*(diffV\diffM(1:arOrderL)')/arOrderL;

            %get the p-value of the test statistic assuming a chi2 distribution
            pValue = 1 - chi2cdf(testStatistic,arOrderL);

            %compare p-value to significance
            if pValue < significance %if p-value is smaller than probability of type I error
                H(1) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(1) = 0; %cannot reject null hypothesis
            end
            
        end

        %MA test
        if maOrderL ~= 0

            %calculate variance-covariance matrix of MA difference vector
            diffV = [eye(maOrderL) -eye(maOrderL)]*[varCovMat1(1+arOrderL:end,...
                1+arOrderL:end) zeros(maOrderL); zeros(maOrderL) varCovMat2(...
                1+arOrderL:end,1+arOrderL:end)]*[eye(maOrderL) -eye(maOrderL)]';

            %compute testStatistic
            testStatistic = diffM(arOrderL+1:end)*(diffV\diffM(...
                arOrderL+1:end)')/maOrderL;

            %get the p-value of the test statistic assuming a chi2 distribution
            pValue = 1 - chi2cdf(testStatistic,maOrderL);

            %compare p-value to significance
            if pValue < significance %if p-value is smaller than probability of type I error
                H(2) = 1; %reject null hypothesis that the two models are identical
            else %if p-value is larger than probability of type I error
                H(2) = 0; %cannot reject null hypothesis
            end

        end

end


%%%%% ~~ the end ~~ %%%%%

