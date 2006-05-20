function [maParam,errFlag] = levinsonDurbinExpoMA(maParamP)
%LEVINSONDURBINEXPOMA determines MA coefficients from "partial MA coefficients" using Levinson-Durbin recursions.
%
%SYNOPSIS [maParam,errFlag] = levinsonDurbinExpoMA(maParamP)
%
%INPUT  maParamP: Parameters (row vector) from which partial MA 
%                 coefficients are obtained via the equation
%                 partial MA coef. = (1-exp(maParamP))/(1+exp(maParamP))
%
%OUTPUT maParam : moving average coefficients (row vector).
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARKS The recursion used is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980), Eq. 6.4. 
%
%Khuloud Jaqaman, July 2004

%initialize output
maParam = [];
errFlag = [];

%get MA order
maOrder = length(maParamP);

%convert parameters to partial MA coefficients
maParamP = (1-exp(maParamP))./(1+exp(maParamP));

%compute MA coefficients from partial MA coefficients
temp = [];
for i=1:maOrder
    temp = [temp zeros(2,1)];
    temp(2,i) = maParamP(i);
    for j=1:i-1
        temp(2,j) = temp(1,j) + maParamP(i)*temp(1,i-j);
    end
    temp(1,:) = temp(2,:);
end

maParam = temp(2,:);


%%%%% ~~ the end ~~ %%%%%
