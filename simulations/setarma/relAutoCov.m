function [gammaV,errFlag] = relAutoCov(arOrder,maOrder,arParam,maParam,maxLag)
%RELAUTOCOV calculates the autocovariance function of an ARMA model. 
%
%SYNOPSIS [gammaV,errFlag] = relAutoCov(arOrder,maOrder,arParam,maParam,maxLag)
%
%INPUT  arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       maxLag      : Maximum lag at which autocovariance function is
%                     calculated.
%
%OUTPUT gammaV      : Autocovariance/(noise variance) for lags 0 to maxLag
%                     (useful for noise free simulations).
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('relAutoCov')
    disp('--relAutoCov: Incorrect number of input arguments!');
    errFlag  = 1;
    gammaV = [];
    return
end

%check input data
if arOrder < 0
    disp('--relAutoCov: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--relAutoCov: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--relAutoCov: Please fix input data!');
    gammaV = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam);
    if nRow ~= 1
        disp('--relAutoCov: "arParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= arOrder
            disp('--relAutoCov: Wrong length of "arParam"!');
            errFlag = 1;
        end
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--relAutoCov: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam);
    if nRow ~= 1
        disp('--relAutoCov: "maParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= maOrder
            disp('--relAutoCov: Wrong length of "maParam"!');
            errFlag = 1;
        end
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--relAutoCov: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maxLag < 0
    disp('--relAutoCov: "maxLag" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--relAutoCov: Please fix input data!');
    gammaV = [];
    return
end


[maInfParam,errFlag] = arma2ma(arOrder,maOrder,arParam,maParam,2*maxLag);
if errFlag
    gammaV = [];
    return
end

gammaV = zeros(maxLag+1,1);

for lag = 0:maxLag
    gammaV(lag+1) = maInfParam(1:end-lag)'*maInfParam(1+lag:end);
end
