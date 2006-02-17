function [gamma,errFlag] = autoCorrARMA(arParam,maParam,maxLag)
%autoCorrARMA calculates the autocorrelaton function of a causal ARMA model, normalized by its white noise variance. 
%
%SYNOPSIS [gamma,errFlag] = autoCorrARMA(arParam,maParam,maxLag)
%
%INPUT  arParam : Row vector of autoregression coefficients. Must satisfy causality condition.
%       maPAram : Row vector of moving average coefficients.
%       maxLag  : Maximum lag at which autocorrelation function is
%                 calculated.
%
%OUTPUT gamma   : Autocorrelation/(white noise variance).
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used
if nargin ~= nargin('autoCorrARMA')
    disp('--autoCorrARMA: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(arParam)
    [nRow,arOrder] = size(arParam);
    if nRow ~= 1
        disp('--autoCorrARMA: "arParam" should be a row vector!');
        errFlag = 1;
    end
else
    arOrder = 0;
end
if ~isempty(maParam)
    [nRow,maOrder] = size(maParam);
    if nRow ~= 1
        disp('--autoCorrARMA: "maParam" should be a row vector!');
        errFlag = 1;
    end
else
    maOrder = 0;
end
if maxLag < 0
    disp('--autoCorrARMA: "maxLag" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--autoCorrARMA: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of autocorrelation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%write AMRA(p,q) process as an MA(infinity process). 
%(coefficients upto 3*maxLag are computed so that number of MA(infinity)
%coefficients used to calcaulte gamma(maxLag) = 2/3 the number of
%MA(infinity) coefficients used to calculate gamma(0)).
[maInfParam,errFlag] = arma2ma(arParam,maParam,3*maxLag);
if errFlag
    return
end

%initialize autocorrelation
gamma = zeros(maxLag+1,1);

%calculate autocorrelation using Eq. 3.2.3
for lag = 0:maxLag
    gamma(lag+1) = maInfParam(1:end-lag)*maInfParam(1+lag:end)';
end


%%%%% ~~ the end ~~ %%%%%
