function [kappa,errFlag] = transAutoCov(arOrder,maOrder,arParam,maParam,...
    maxLagI,maxLagJ,checkRoots)
%TRANSAUTOCOV calculates the autocovariance function for the innovations algorithm (see Remark)
%
%SYNOPSIS [kappa,errFlag] = transAutoCov(arOrder,maOrder,arParam,maParam,...
%    maxLagI,maxLagJ,checkRoots)
%
%INPUT  arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       maxLag      : Maximum lag at which autocovariance function is
%                     calculated.
%       checkRoots  : 1 if AR polynomial roots are to be checked, 0 otherwise.
%
%OUTPUT kappa       : Autocovariance matrix for lags 0 to maxLag.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARK This function computes the autocovariance function of the transformed process 
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
if arOrder < 1
    disp('--transAutoCov: "arOrder" should be >= 1!');
    errFlag = 1;
end
if maOrder < 1
    disp('--transAutoCov: "maOrder" should be >= 1!');
    errFlag = 1;
end
if errFlag
    disp('--transAutoCov: Please fix input data!');
    kappa = [];
    return
end
[nRow,nCol] = size(arParam);
if nRow ~= 1
    disp('--transAutoCov: "arParam" should be a row vector!');
    errFlag = 1;
end
if nCol ~= arOrder
    disp('--transAutoCov: Wrong length of "arParam"!');
    errFlag = 1;
end
[nRow,nCol] = size(maParam);
if nRow ~= 1
    disp('--transAutoCov: "maParam" should be a row vector!');
    errFlag = 1;
end
if nCol ~= maOrder
    disp('--transAutoCov: Wrong length of "maParam"!');
    errFlag = 1;
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


%compute (autocovariance function)/(noise varaince) of an ARMA(p,q) process
maxLag = max(maxLagI,maxLagJ);
[gammaV,errFlag] = relAutoCov(arOrder,maOrder,arParam,maParam,maxLag,checkRoots);
gammaV = [gammaV(end:-1:2); gammaV(1); gammaV(2:end)];

maxOrder = max(arOrder,maOrder);  %needed for calculation of kappa

%initialize kappa (autocovariance matrix of transformed process)
kappa = zeros(maxLagI+1,maxLagJ+1);

for lagI = 0:maxLagI
    for lagJ = 0:maxLagJ
        
        minLocal = min(lagI,lagJ); %minimum of pair
        maxLocal = max(lagI,lagJ); %maximum of pair
        lagDiff = abs(lagI-lagJ);  %absolute difference between the two lags

        %implementation of Eq. 3.3.3
        if lagI>=1 && lagI<=maxOrder && lagJ>=1 && lagJ<=maxOrder

            kappa(lagI+1,lagJ+1) = gammaV(maxLag+1+lagDiff);
            
        elseif minLocal<=maxOrder && maxOrder<maxLocal && maxLocal<=2*maxOrder
            
            kappa(lagI+1,lagJ+1) = gammaV(maxLag+1+lagDiff);
            if arOrder ~= 0
                kappa(lagI+1,lagJ+1) = kappa(lagI+1,lagJ+1)-(arParam*...
                    gammaV(maxLag+2-lagDiff:maxLag+1+arOrder-lagDiff));
            end
           
        elseif minLocal > maxOrder && maOrder ~= 0

            maPoly1 = [1 maParam];
            maPoly2 = [maPoly1(1+lagDiff:end) zeros(1,min(maOrder+1,lagDiff))]';
            kappa(lagI+1,lagJ+1) = maPoly1*maPoly2;
            
        end
        
    end
end
