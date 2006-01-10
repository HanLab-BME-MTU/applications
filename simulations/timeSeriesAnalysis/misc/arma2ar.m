function [arInfParam,errFlag] = arma2ar(arParam,maParam,numARCoef)
%ARMA2AR converts an invertible ARMA(p,q) model into an AR(infinity) model. 
%
%SYNOPSIS [arInfParam,errFlag] = arma2ar(arParam,maParam,numARCoef)
%
%INPUT  arParam    : Row vector of autoregression coefficients.
%       maPAram    : Row vector of moving average coefficients.
%       numARCoef  : Number of coefficients in AR(infinity) model to be
%                     calculated. Optional. Default: 50.
%
%OUTPUT arInfParam : Parameters in AR(infinity) model.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm used here is taken from Brockwell and Davis,
%"Introduction to Time Series and Forecasting", 2nd ed., p. 86.
%
%Khuloud Jaqaman, September 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arInfParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used
if nargin < 2 || nargin > 3
    disp('--arma2ar: You should input AR and MA coefficients!')
    disp('           Number of AR coefficients to calculate, however, is optional!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(arParam)
    [nRow,arOrder] = size(arParam);
    if nRow ~= 1
        disp('--arma2ar: "arParam" should be a row vector!');
        errFlag = 1;
    end
else
    arOrder = 0;
end
if ~isempty(maParam)
    [nRow,maOrder] = size(maParam);
    if nRow ~= 1
        disp('--arma2ar: "maParam" should be a row vector!');
        errFlag = 1;
    else
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--arma2ar: Invertibility requires the polynomial defining the moving average part of the model to have no zeros for z <= 1!');
            errFlag = 1;
        end
    end
else
    maOrder = 0;
end

%assign default number of AR coefficients if user did not specify a number
if nargin == 2
    numARCoef = 50;
else
    if numARCoef <= 0
        disp('--arma2ar: "numARCoef" should be a positive integer!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--arma2ar: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conversion from ARMA(p,q) to AR(infinity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxOrder = max(arOrder,maOrder);

%rewrite MA part of ARMA process
arPoly = zeros(1,numARCoef);
arPoly(1) = -1;
arPoly(2:arOrder+1) = arParam;

%initialize AR(infinity) coefficients
arInfParam = zeros(1,numARCoef);

%calculate coefficients using Eq. 3.1.8
if maOrder ~= 0 %if there is an MA part
    arInfParam(1) = 1; %zeroth order term
    for i=2:numARCoef %index j in Eq. 3.1.8
        j1 = min(i-1,maOrder); %index of MA coef in Eq. 3.1.8 goes up from 1 to j1
        j2 = max(1,i-maOrder); %index of arInfParam in Eq. 3.1.8 goes down from j-1 (i-1 here) to j2
        arInfParam(i) = -arPoly(i) - maParam(1:j1)*arInfParam(i-1:-1:j2)';
    end
else %if there is no MA part
    arInfParam(1:maOrder+1) = [1 -arParam]; %trivial case
end


%%%%% ~~ the end ~~ %%%%%
