function [maInfParam,errFlag] = arma2ma(arOrder,maOrder,arParam,maParam,upLim,...
    checkRoots)
%MAINFPARAM converts a causal ARMA(p,q) model into a MA(infinity) model. 
%
%SYNOPSIS [maInfParam,errFlag] = arma2ma(arOrder,maOrder,arParam,maParam,upLim,...
%    checkRoots)
%
%INPUT  arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       upLim       : Number of coefficients in MA(infinity) model to be
%                     calculated.
%       checkRoots  : 1 if AR polynomial roots are to be checked, 0 otherwise.
%
%OUTPUT maInfParam  : Parameters in MA(infinity) model
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arma2ma')
    disp('--arma2ma: Incorrect number of input arguments!');
    errFlag  = 1;
    maInfParam = [];
    return
end

%check input data
if arOrder < 0
    disp('--arma2ma: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--arma2ma: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if errFlag
    disp('--arma2ma: Please fix input data!');
    maInfParam = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam);
    if nRow ~= 1
        disp('--arma2ma: "arParam" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= arOrder
            disp('--arma2ma: Wrong length of "arParam"!');
            errFlag = 1;
        end
        if checkRoots
            r = abs(roots([-arParam(end:-1:1) 1]));
            if ~isempty(find(r<=1))
                disp('--arma2ma: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
                errFlag = 1;
            end
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam);
    if nRow ~= 1
        disp('--arma2ma: "maParam" should be a row vector!');
        errFlag = 1;
    end
    if nCol ~= maOrder
        disp('--arma2ma: Wrong length of "maParam"!');
        errFlag = 1;
    end
end
if upLim <= 0
    disp('--arma2ma: "upLim" should be a positive integer!');
    errFlag = 1;
end
if errFlag
    disp('--arma2ma: Please fix input data!');
    maInfParam = [];
    return
end

maxOrder = max(arOrder,maOrder);

%rewrite MA part of ARMA process
maPoly = zeros(upLim,1);
maPoly(1) = 1;
maPoly(2:maOrder+1) = maParam;

%initialize MA(infinity) coefficients
maInfParam = zeros(upLim,1);

%calculate coefficients using Eq. 3.1.7
if arOrder ~= 0 %if there is an AR part
    maInfParam(1) = 1; %zeroth order term
    for i=2:upLim
        j1 = min(i-1,arOrder);
        j2 = max(1,i-arOrder);
        maInfParam(i) = maPoly(i) + arParam(1:j1)*maInfParam(i-1:-1:j2);
    end
else %if there is no AR part
    maInfParam(1:maOrder+1) = [1 maParam];
end
