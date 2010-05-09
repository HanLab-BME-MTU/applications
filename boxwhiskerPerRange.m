function [delta] = boxwhiskerPerRange(tau, k, pct)
% Returns the range of for a box-and-whisker plot for Weibull distribution.
% In the interval [tau-delta, tau+delta], the integral of the distribution is equal to 'pct'.

% Francois Aguet, April 2010
% For k=1 the exact solution is: delta = tau*asinh(exp(1)*pct(n)/2).

if tau <= 0
    error('''tau'' must be greater than zero.');
end
if k <= 0
    error('''k'' must be greater than zero.');
end

f = @(a) exp(-(1-a/tau)^k) - exp(-(1+a/tau)^k) - pct;
delta = fzero(f, tau);