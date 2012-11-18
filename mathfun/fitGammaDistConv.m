%
% Implements the form corresponding to convolution of 'n' steps with rate 'k'

% Francois Aguet, 11/17/2012

function [k, n, a, BIC] = fitGammaDistConv(x, f1, f2, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x', @isvector);
ip.addRequired('f1', @isvector);
ip.addRequired('f2', @isvector);
ip.addParamValue('Mode', 'kna', @ischar);
ip.addParamValue('k', NaN, @isscalar);
ip.addParamValue('n', NaN, @isscalar);
ip.addParamValue('a', NaN, @isscalar);
ip.parse(x, f1, f2, varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

p0 = [ip.Results.k ip.Results.n];% ip.Results.a];

estIdx = false(1,2); % [k n a]
estIdx(regexp('kn', ['[' ip.Results.Mode ']'])) = true;

mu = sum(f1/sum(f1).*x);
prmVect = [1/mu 2];
idx = ~isnan(p0);
prmVect(idx) = p0(idx);

F1 = fft(f1);
N = sum(estIdx);
[p,RSS] = lsqnonlin(@cost, prmVect(estIdx), zeros(N,1), inf(N,1), opts, x, F1, f2, prmVect, estIdx);
prmVect(estIdx) = p;
k = prmVect(1);
n = prmVect(2);
% a = prmVect(3);

d = numel(x);
BIC = d*log(RSS/d) + N*log(d); 

fprintf('k = %.3f, n = %.3f, BIC = %.3f\n', k, n, BIC);

dx = x(2)-x(1);

g = exp(-k*x)*k^n.*x.^(n-1)/gamma(n);
g = g/sum(g)/dx;
G = fft(g);
fx = ifft(G.*F1)*dx;

figure;
hold on;
plot(x, f1, 'k');
plot(x, f2, 'b');
plot(x, fx, 'r--');


function v = cost(p, x, F1, f2, prmVect, estIdx)
dx = x(2)-x(1);
prmVect(estIdx) = p;
k = prmVect(1);
n = prmVect(2);

g = exp(-k*x)*k^n.*x.^(n-1)/gamma(n);
g = g/sum(g)/dx;

G = fft(g);

v = ifft(G.*F1)*dx - f2;
