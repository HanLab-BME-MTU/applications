
%    prmVect : Parameter vector with components [k1 n1 A1 ... ], where
%                   k : rate parameter
%                   n : shape parameter
%                   A : contribution of each component
%
% Mean: n/k

function [w, W, J] = gammaMixture(x, prmVect, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('prmVect', @(x) mod(numel(x),3)==0);
ip.addParamValue('Mode', 'CDF', @(x) any(strcmpi(x, {'PDF','CDF'})));
ip.parse(x, prmVect, varargin{:});

N = numel(prmVect)/3;
prmVect = abs(prmVect);
k = prmVect(1:3:end);
n = prmVect(2:3:end);
A = prmVect(3:3:end);


W = zeros(N,numel(x));
switch ip.Results.Mode
    case 'PDF'
        J = zeros(numel(x), numel(prmVect));
        for i = 1:N
            W(i,:) = k(i)^n(i)*x.^(n(i)-1).*exp(-k(i)*x)/gamma(n(i));
            J(:,3*(i-1)+1) = A(i) * k(i)^(n(i)-1)*x.^(n(i)-1).*exp(-k(i)*x)/gamma(n(i)).*(n(i)-k(i)*x);
            J(:,3*(i-1)+2) = W(i,:)*(log(k(i))+log(n(i))-psi(n(i)));
            J(:,3*(i-1)+3) = k(i)^n(i)*x.^(n(i)-1).*exp(-k(i)*x)/gamma(n(i));
        end
    case 'CDF'
        J = [];
        for i = 1:N
            W(i,:) = gammainc(k(i)*x, n(i), 'lower');
            %J(:,3*(i-1)+1) = 
            %J(:,3*(i-1)+2) = 
            %J(:,3*(i-1)+3) = gammainc(k(i)*x, n(i), 'lower');
        end
end

w = A*W;
