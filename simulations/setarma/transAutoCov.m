function [kappa,errFlag] = transAutoCov(arOrder,maOrder,arParam,maParam,...
    variance,maxLagI,maxLagJ)
%TRANSAUTOCOV calculates the autocovariance function for the innovations algorithm (see Remark)
%
%SYNOPSIS [kappa,errFlag] = transAutoCov(arOrder,maOrder,arParam,maParam,...
%    variance,maxLagI,maxLagJ)
%
%INPUT  arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       variance    : Variance of white noise in model. 
%       maxLag      : Maximum lag at which autocovariance function is
%                     calculated.
%
%OUTPUT kappa       : Autocovariance matrix for lags 0 to maxLag.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARK This function computed the autocovariance function of the transformed process 
%defined in Eq. 3.3.1 of "Introduction to Time Series and Forecasting" by
%Brockwell and Davis. This transformation simplifies the use of the
%innovations algorithm for 1-step prediction of an ARMA(p,q) process.

%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('transAutoCov')
    disp('--transAutoCov: Incorrect number of input arguments!');
    errFlag  = 1;
    kappa = [];
    return
end

%check input data
if arOrder < 0
    disp('--transAutoCov: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--transAutoCov: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--transAutoCov: Please fix input data!');
    kappa = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam);
    if nRow ~= 1
        disp('--transAutoCov: "arParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= arOrder
            disp('--transAutoCov: Wrong length of "arParam"!');
            errFlag = 1;
        end
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--transAutoCov: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam);
    if nRow ~= 1
        disp('--transAutoCov: "maParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= maOrder
            disp('--transAutoCov: Wrong length of "maParam"!');
            errFlag = 1;
        end
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--transAutoCov: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if variance < 0
    disp('--transAutoCov: White noise variance should be nonnegative!');
    errFlag;
end
if maxLagI < 0
    disp('--transAutoCov: "maxLagI" should be a nonnegative integer!');
    errFlag = 1;
end
if maxLagJ < 0
    disp('--transAutoCov: "maxLagJ" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--transAutoCov: Please fix input data!');
    kappa = [];
    return
end


maxLag = max(maxLagI,maxLagJ);
[gamma,gammaV,errFlag] = autoCov(arOrder,maOrder,arParam,maParam,variance,maxLag);
gammaV = [gammaV(end:-1:2); gammaV(1); gammaV(2:end)];

maxOrder = max(arOrder,maOrder);

kappa = zeros(maxLagI+1,maxLagJ+1);

for lagI = 0:maxLagI
    for lagJ = 0:maxLagJ
        
        minLocal = min(lagI,lagJ);
        maxLocal = max(lagI,lagJ);
        lagDiff = abs(lagI-lagJ);

        if lagI>=1 && lagI<=maxOrder && lagJ>=1 && lagJ<=maxOrder

            kappa(lagI+1,lagJ+1) = gammaV(maxLag+1+lagDiff);
            
        elseif minLocal<=maxOrder && maxOrder<maxLocal && maxLocal<=2*maxOrder
            
            kappa(lagI+1,lagJ+1) = gammaV(maxLag+1+lagDiff);
            if arOrder ~= 0
                kappa(lagI+1,lagJ+1) = kappa(lagI+1,lagJ+1)-(arParam*...
                    gammaV(maxLag+2-lagDiff:maxLag+1+arOrder-lagDiff));
            end
           
        elseif min(lagI,lagJ) > maxOrder && maOrder ~= 0

            maPoly1 = [1 maParam];
            maPoly2 = [maPoly1(1+lagDiff:end) zeros(1,min(maOrder+1,lagDiff))]';
            kappa(lagI+1,lagJ+1) = maPoly1*maPoly2;
            
        end
        
    end
end
