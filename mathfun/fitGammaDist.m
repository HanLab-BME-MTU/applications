% 
% Implements the form corresponding to convolution of 'n' steps with rate 'k'

function [k, n, x, f, a, K] = fitGammaDist(samples)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);


% EDF of the samples
[f_ecdf, x_ecdf] = ecdf(samples);

[p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, [2/mean(samples) 2], [], [], opts, x_ecdf, f_ecdf);
k = p(1);
n = p(2);

% error propagation, parameter correlations
J = full(J);


C = resnorm/(numel(samples)-numel(p)-1)*inv(J'*J); %#ok<MINV>
param_pstd = sqrt(diag(C))';
K = corrFromC(C)';
K = K(2,1);

% area corresponding to missing data
a = 1-gammainc(k*x_ecdf(1),n, 'lower');

x = linspace(0, x_ecdf(end), 1000);
f = k^n*x.^(n-1).*exp(-k*x)/gamma(n);


function v = cost(p, x, f_ecdf)
k = p(1);
n = p(2);
cdf = gammainc(k*x,n, 'lower');

% normalization for missing data
T = cdf(1);
cdf = (cdf-T)/(1-T);
v = cdf - f_ecdf;



function K = corrFromC(C)
n = size(C,1);
K = zeros(n,n);

idx = pcombs(1:n);
i = idx(:,1);
j = idx(:,2);
ij = i+n*(j-1);
ii = i+n*(i-1);
jj = j+n*(j-1);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));