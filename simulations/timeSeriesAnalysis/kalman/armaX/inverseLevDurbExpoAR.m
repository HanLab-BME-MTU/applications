function [arParamP,errFlag] = inverseLevDurbExpoAR(arParam)
%INVERSELEVDURBEXPOAR determines "partial AR coefficients" from AR coefficients by inverting the Levinson-Durbin recursions.
%
%SYNOPSIS [arParamP,errFlag] = inverseLevDurbExpoAR(arParam)
%
%INPUT  arParam : Autoregressive coefficients (row vector).
%
%OUTPUT 
%       arParamP: Parameters (row vector) from which partial AR 
%                 coefficients are obtained via the equation
%                 partial AR coef. = (1-exp(arParamP))/(1+exp(arParamP))
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARKS The recursions used are the inverse of the recursions presented in
%        R. H. Jones, "Maximum Likelihood Fitting of ARMA Models to Time 
%        Series with Missing Observations", Technometrics 22: 389-395 (1980), 
%        Eq. 6.2. 
%
%Khuloud Jaqaman, April 2006

%initialize output
arParamP = [];
errFlag  = [];

%get AR order
arOrder = length(arParam);

%do the recursion
temp = arParam;
for i=arOrder:-1:1
    arParamP(i) = temp(end);
    multFact = 1/(1-arParamP(i)^2);
    sol = [];
    for j=1:i-1
        sol(j) = multFact*(arParamP(i)*temp(i-j) + temp(j));
    end
    temp = sol;
end

%convert partial AR coefficients to output
arParamP = log((1-arParamP)./(1+arParamP));


%%%%% ~~ the end ~~ %%%%%
