function [innovErr,innovCoef,errFlag] = innovationsAlg(kappa,maxIndex)
%INNOVATIONSALG calculates the parameters used in the innovations algorithm
%
%SYNOPSIS [innovErr,innovCoef,errFlag] = innovationsAlg(kappa,maxIndex)
%
%INPUT  kappa       : Autocovariance matrix of transformed series,
%                     calculated via "tranAutoCov"
%       maxIndex    : maximum index of series whose values are to be
%                     predicted.
%
%OUTPUT innovErr    : Mean squared errors of prediction.
%       unnovCoef   : Prediction coefficients in algorithm.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('innovationsAlg')
    disp('--innovationsAlg: Incorrect number of input arguments!');
    errFlag  = 1;
    kappa = [];
    return
end

%check input data
if maxIndex <= 0
    disp('--innovationsAlg: "maxIndex" should be a positive integer!');
    errFlag = 1;
end
[nRow,nCol] = size(kappa);
if nRow < maxIndex + 2
    disp('--innovationsAlg: "kappa" should have more than "maxIndex+2" rows!');
    errFlag = 1;
end
if nCol < maxIndex + 2
    disp('--innovationsAlg: "kappa" should have more than "maxIndex+2" columns!');
    errFlag = 1;
end

innovErr = zeros(maxIndex+1,1);
innovCoef = zeros(maxIndex+1);

innovErr(1) = kappa(2,2);

for n = 1:maxIndex
    
    innovCoef(n,n) = kappa(n+2,2)/innovErr(1);
    
    for k=1:n-1        
        innovCoef(n,n-k) = (kappa(n+2,k+2)-sum(innovCoef(k,k:-1:1).*...
            innovCoef(n,n:-1:n-k+1).*innovErr(1:k)'))/innovErr(k+1);
    end
    
    innovErr(n+1) = kappa(n+2,n+2)-(innovCoef(n,n:-1:1).^2)*innovErr(1:n);
    
end
