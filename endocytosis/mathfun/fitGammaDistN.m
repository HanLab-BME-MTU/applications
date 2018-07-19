% 
% Implements the form corresponding to convolution of 'n' steps with rate 'k'

% samples is a cell array
function [k, n, x, f, F, a] = fitGammaDistN(samples)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

nd = numel(samples);

% EDF of the samples
x_ecdf = cell(1,nd);
f_ecdf = cell(1,nd);
for i = 1:nd
    [f_ecdf{i}, x_ecdf{i}] = ecdf(samples{i});
end
    
[p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, [2/mean(samples{1}) 2*ones(1,nd)], [], [], opts, x_ecdf, f_ecdf);
k = p(1);
n = p(2:end);

% error propagation, parameter correlations
% J = full(J);


% C = resnorm/(numel(samples)-numel(p)-1)*inv(J'*J); %#ok<MINV>
% param_pstd = sqrt(diag(C))';
% K = corrMatFromCov(C)';
% K = K(2,1);


x = linspace(0, x_ecdf{nd}(end), 1000);
a = ones(1,nd);
f = cell(1,nd);
F = cell(1,nd);
for i = 1:nd
    % area corresponding to data
    a(i) = 1-gammainc(k*x_ecdf{i}(1),n(i), 'lower');
    % pdf
    f{i} = k^n(i)*x.^(n(i)-1).*exp(-k*x)/gamma(n(i));
    F{i} = gammainc(k*x,n(i), 'lower');
end



function v = cost(p, x_ecdf, f_ecdf)
k = p(1);
n = p(2:end);
nd = numel(p)-1;

v = cell(1,nd);
for i = 1:nd
    cdf = gammainc(k*x_ecdf{i},n(i), 'lower');
    % normalization for missing data
%     T = cdf(1);
%     cdf = (cdf-T)/(1-T);
    v{i} = cdf - f_ecdf{i};
end
v = vertcat(v{:});
