
%    prmVect : Parameter vector with components [A1 lambda1 k1 ... ], where
%              lambda : scale parameter
%                   k : shape parameter
%                   A : contribution of each component
%
% Mean: lambda*Gamma(1+1/k)

function [w W J] = weibullMixture(x, prmVect, varargin)

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
J = zeros(numel(x), numel(prmVect));
switch ip.Results.Mode
    case 'PDF'
        for i = 1:N
            xl = x/lambda(i);
            W(i,:) = k(i)/lambda(i) * xl.^(k(i)-1).*exp(-xl.^k(i)); % *A(i)
            J(:,3*(i-1)+1) = A(i) * exp(-xl.^k(i))*k(i)^2.*(-1+xl.^k(i)).*xl.^(k(i)-1)/lambda(i)^2;
            J(:,3*(i-1)+2) = A(i) * exp(-xl.^k(i)).*xl.^k(i).*(1-k(i)*(-1+xl.^k(i)).*log(xl))./x;
            J(:,3*(i-1)+3) = k(i)/lambda(i) * xl.^(k(i)-1).*exp(-xl.^k(i));
        end
    case 'CDF'
        for i = 1:N
            xl = x/lambda(i);
            W(i,:) = (1 - exp(-xl.^k(i))); % * A(i)
            J(:,3*(i-1)+1) = A(i) * exp(-xl.^k(i))*k(i).*xl.^k(i)/lambda(i);
            J(:,3*(i-1)+2) = A(i) * exp(-xl.^k(i)).*xl.^k(i).*log(x./lambda(i));
            J(:,3*(i-1)+3) = (1 - exp(-(x/lambda(i)).^k(i)));
        end
end

w = A*W;