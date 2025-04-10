% [w W] = nWeibull(t, prmVect, distribution)
% Computes the sum of N Weibull distributions.
%
% Inputs:
%           t              : time vector
%           prmVect        : parameter vector of size 3*N+1
%           (distribution) : 'PDF' (default), 'CDF', or 'SVF'
%
% Outputs:  w : Summed distributions 
%           W : Individual underlying distributions (without offset)
%
% Structure of the parameter vector: [offset A(1) lambda(1) k(1) ... A(N) lambda(N) k(N)]

% Francois Aguet, Jan 2010


function [w W] = nWeibull(t, prmVect, distribution)

if nargin<3
    distribution = 'PDF';
end

if ~strcmp(distribution, 'PDF') && ~strcmp(distribution, 'CDF') && ~strcmp(distribution, 'SVF')
    error('err:thise', 'Valid ''functionType'':\n   ''PDF'' (probability density function)\n   ''CDF'' (cumulative distribution function)\n   ''SVF'' (survival function)');
end;

N = (length(prmVect)-1)/3;
offset = prmVect(1);
A = abs(prmVect(2:3:end));

lambda = abs(prmVect(3:3:end));
k = abs(prmVect(4:3:end));

W = zeros(N, length(t));
switch distribution
    case 'PDF'
        for n = 1:N
            tl = t/lambda(n);
            W(n,:) = A(n) * k(n)/lambda(n) * tl.^(k(n)-1).*exp(-tl.^k(n));
        end
    case 'CDF'
        for n = 1:N
            W(n,:) = A(n) * (1 - exp(-(t/lambda(n)).^k(n)));
        end
    case 'SVF'
        for n = 1:N
            W(n,:) = A(n) * exp(-(t/lambda(n)).^k(n));
        end
end
w = sum(W, 1) + offset;