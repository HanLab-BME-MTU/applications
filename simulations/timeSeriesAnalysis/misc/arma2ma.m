function [maInfParam,errFlag] = arma2ma(arParam,maParam,numMACoef)
%ARMA2MA converts a causal ARMA(p,q) model into a MA(infinity) model. 
%
%SYNOPSIS [maInfParam,errFlag] = arma2ma(arParam,maParam,numMACoef)
%
%INPUT  arParam    : Row vector of autoregression coefficients.
%       maPAram    : Row vector of moving average coefficients.
%       numMACoef  : Number of coefficients in MA(infinity) model to be
%                     calculated. Optional. Default: 50.
%
%OUTPUT maInfParam : Parameters in MA(infinity) model.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm used here is taken from Brockwell and Davis,
%"Introduction to Time Series and Forecasting", 2nd ed., pp. 84-85.
%
%Khuloud Jaqaman, February 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maInfParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used
if nargin < 2 || nargin > 3
    disp('--arma2ma: You should input AR and MA coefficients!')
    disp('           Number of MA coefficients to calculate, however, is optional!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(arParam)
    [nRow,arOrder] = size(arParam);
    if nRow ~= 1
        disp('--arma2ma: "arParam" should be a row vector!');
        errFlag = 1;
    else
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--arma2ma: Causality requires the polynomial defining the autoregressive part of the model to have no zeros for z <= 1!');
            errFlag = 1;
        end
    end
else
    arOrder = 0;
end
if ~isempty(maParam)
    [nRow,maOrder] = size(maParam);
    if nRow ~= 1
        disp('--arma2ma: "maParam" should be a row vector!');
        errFlag = 1;
    end
else
    maOrder = 0;
end

%assign default number of MA coefficients if user did not specify a number
if nargin == 2
    numMACoef = 50;
else
    if numMACoef <= 0
        disp('--arma2ma: "numMACoef" should be a positive integer!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--arma2ma: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conversion from ARMA(p,q) to MA(infinity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxOrder = max(arOrder,maOrder);

%rewrite MA part of ARMA process
maPoly = zeros(1,numMACoef);
maPoly(1) = 1;
maPoly(2:maOrder+1) = maParam;

%initialize MA(infinity) coefficients
maInfParam = zeros(1,numMACoef);

%calculate coefficients using Eq. 3.1.7
if arOrder ~= 0 %if there is an AR part
    maInfParam(1) = 1; %zeroth order term
    for i=2:numMACoef %index j in Eq. 3.1.7
        j1 = min(i-1,arOrder); %index of AR coef in Eq. 3.1.7 goes up from 1 to j1
        j2 = max(1,i-arOrder); %index of maInfParam in Eq. 3.1.7 goes down from j-1 (i-1 here) to j2
        maInfParam(i) = maPoly(i) + arParam(1:j1)*maInfParam(i-1:-1:j2)';
    end
else %if there is no AR part
    maInfParam(1:maOrder+1) = [1 maParam]; %trivial case
end


%%%%% ~~ the end ~~ %%%%%

