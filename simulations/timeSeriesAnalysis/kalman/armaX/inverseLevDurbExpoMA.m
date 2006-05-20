function [maParamP,errFlag] = inverseLevDurbExpoMA(maParam)
%INVERSELEVDURBEXPOAR determines "partial AR coefficients" from AR coefficients by inverting the Levinson-Durbin recursions.
%
%SYNOPSIS [maParamP,errFlag] = inverseLevDurbExpoMA(maParam)
%
%INPUT  maParam : Moving average coefficients (row vector).
%
%OUTPUT 
%       maParamP: Parameters (row vector) from which partial MA 
%                 coefficients are obtained via the equation
%                 partial MA coef. = (1-exp(maParamP))/(1+exp(maParamP))
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARKS The recursions used are the inverse of the recursions presented in
%        R. H. Jones, "Maximum Likelihood Fitting of ARMA Models to Time 
%        Series with Missing Observations", Technometrics 22: 389-395 (1980), 
%        Eq. 6.4. 
%
%Khuloud Jaqaman, April 2006

%initialize output
maParamP = [];
errFlag  = [];

%get MA order
maOrder = length(maParam);

%do the recursion
temp = maParam;
for i=maOrder:-1:1
    maParamP(i) = temp(end);
    multFact = 1/(maParamP(i)^2-1);
    sol = [];
    for j=1:i-1
        sol(j) = multFact*(maParamP(i)*temp(i-j) - temp(j));
    end
    temp = sol;
end

%convert partial AR coefficients to output
maParamP = log((1-maParamP)./(1+maParamP));


%%%%% ~~ the end ~~ %%%%%
