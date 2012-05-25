
%    prmVect : Parameter vector with components [A1 lambda1 k1 ... ], where
%              lambda : scale parameter
%                   k : shape parameter
%                   A : contribution of each component
%
% Mean: lambda*Gamma(1+1/k)

function [w W] = weibullMixture(x, prmVect, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('prmVect', @(x) mod(numel(x),3)==0);
ip.addParamValue('Mode', 'CDF', @(x) any(strcmpi(x, {'PDF','CDF'})));
ip.parse(x, prmVect, varargin{:});

N = numel(prmVect)/3;
prmVect = abs(prmVect);
lambda = prmVect(1:3:end);
k = prmVect(2:3:end);
A = prmVect(3:3:end);


W = zeros(N,numel(x));
switch ip.Results.Mode
    case 'PDF'
        for n = 1:N
            xl = x/lambda(n);
            W(n,:) = A(n) * k(n)/lambda(n) * xl.^(k(n)-1).*exp(-xl.^k(n));
        end
    case 'CDF'
        for n = 1:N
            W(n,:) = A(n) * (1 - exp(-(x/lambda(n)).^k(n)));
        end
end

w = sum(W, 1);