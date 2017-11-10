% 
% Implements the form corresponding to convolution of 'n' steps with rate 'k'

function [k, n, x, f, a, K] = fitGammaDist(in1, in2)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);


if nargin==1
    samples = in1;

    % EDF of the samples
    [f_ecdf, x_ecdf] = ecdf(samples);
    
    [p,resnorm,~,~,~,~,J] = lsqnonlin(@costEDF, [2/mean(samples) 2], [], [], opts, x_ecdf, f_ecdf);
    k = p(1);
    n = p(2);
    
    % error propagation, parameter correlations
    J = full(J);
    
    
    C = resnorm/(numel(samples)-numel(p)-1)*inv(J'*J); %#ok<MINV>
    param_pstd = sqrt(diag(C))';
    K = corrMatFromCov(C)';
    K = K(2,1);
    
    % area corresponding to missing data
    a = 1-gammainc(k*x_ecdf(1),n, 'lower');
    
    x = linspace(0, x_ecdf(end), 1000);
    f = k^n*x.^(n-1).*exp(-k*x)/gamma(n);
    
else
    x = in1;
    f = in2;
    dx = x(2)-x(1);
    f = f/sum(f)/dx;

    lb = zeros(2,1);
    ub = Inf(2,1);
    mu = sum(f/sum(f).*x);
    prmVect = [1/mu 1];
      
    [p,RSS] = lsqnonlin(@cost, prmVect, lb, ub, opts, x, f);
    k = p(1);
    n = p(2);
    f = k^n*x.^(n-1).*exp(-k*x)/gamma(n);
    a = sum(f)*dx;
    fprintf('k = %.3f, n = %.3f, a = %.3f\n', k, n, a);
    

end
    

function v = cost(p, x, f)
k = p(1);
n = p(2);
dx = x(2)-x(1);
g = exp(-k*x)*k^n.*x.^(n-1)/gamma(n);

v = g/sum(g)/dx - f;



function v = costEDF(p, x, f_ecdf)
k = p(1);
n = p(2);
cdf = gammainc(k*x,n, 'lower');

% normalization for missing data
T = cdf(1);
cdf = (cdf-T)/(1-T);
v = cdf - f_ecdf;
