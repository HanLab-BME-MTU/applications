% Francois Aguet, 09/06/2011

function [k, k_pstd, t_edf, f_edf, cdf] = fitExpCDF(samples)

opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

[f_edf, t_edf] = ecdf(samples);
k_init = 1/mean(samples);

[k,resnorm,~,~,~,~,J] = lsqnonlin(@cost, k_init, 0, Inf, opts, t_edf, f_edf);

C = resnorm*full(inv(J'*J));
k_pstd = sqrt(diag(C)/(numel(samples)-2));

cdf = 1 - exp(-k*t_edf);


function [v J] = cost(k, t_edf, f_edf)
v = 1 - exp(-k*t_edf) - f_edf;
J = t_edf.*exp(-k*t_edf);
