function [gamma,gammaV,errFlag] = autoCov(arOrder,maOrder,arParam,maParam,variance,maxLag)
%AUTOCOV calculates the autocovariance function of an ARMA(p,q) model. 
%
%SYNOPSIS [gamma,gammaV,errFlag] = autoCov(arOrder,maOrder,arParam,maParam,variance,maxLag)
%
%INPUT  arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       variance    : Variance of white noise in model. 
%       maxLag      : Maximum lag at which autocovariance function is
%                     calculated.
%
%OUTPUT gamma       : Autocovariance function for lags 0 to maxLag.
%       gammaV      : Autocovariance/noise variance for lags 0 to maxLag
%                     (useful for noise free simulations).
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('autoCov')
    disp('--autoCov: Incorrect number of input arguments!');
    errFlag  = 1;
    gamma = [];
    gammaV = [];
    return
end

%check input data
if arOrder < 0
    disp('--autoCov: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--autoCov: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--autoCov: Please fix input data!');
    gamma = [];
    gammaV = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam);
    if nRow ~= 1
        disp('--autoCov: "arParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= arOrder
            disp('--autoCov: Wrong length of "arParam"!');
            errFlag = 1;
        end
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--autoCov: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam);
    if nRow ~= 1
        disp('--autoCov: "maParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= maOrder
            disp('--autoCov: Wrong length of "maParam"!');
            errFlag = 1;
        end
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--autoCov: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if variance < 0
    disp('--autoCov: White noise variance should be nonnegative!');
    errFlag;
end
if maxLag < 0
    disp('--autoCov: "maxLag" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--autoCov: Please fix input data!');
    gamma = [];
    gammaV = [];
    return
end


[maInfParam,errFlag] = arma2ma(arOrder,maOrder,arParam,maParam,2*maxLag);
if errFlag
    gamma = [];
    gammaV = [];
    return
end

gammaV = zeros(maxLag+1,1);

for lag = 0:maxLag
    gammaV(lag+1) = sum(maInfParam(1:end-lag).*maInfParam(1+lag:end));
end

gamma = variance*gammaV;
