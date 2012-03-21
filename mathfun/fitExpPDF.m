function [mu, mu_std, a] = fitExpPDF(samples)


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

[f_edf, t_edf] = ecdf(samples);

% figure;
% hold on;
% plot(t_edf, f_edf);

[mu,resnorm,~,~,~,~,J] = lsqnonlin(@costFct, mean(samples), [], [], opts, t_edf, f_edf);

J = full(J);
mu_std = sqrt( resnorm/(numel(samples)-2) * inv(J'*J) ); %#ok<MINV>
a = exp(-1./mu.*t_edf(1)); % area of data

% figure;
% hold on;
% plot(t_edf, f_edf);
% cdf = 1 - exp(-1./mu.*t_edf);
% T = cdf(1);
% cdf = (cdf - T) / (1-T);
% 
% plot(t_edf, cdf, 'r--');


function v = costFct(mu, t_edf, f_edf)
cdf = 1 - exp(-1./mu.*t_edf);
T = cdf(1);
cdf = (cdf - T) / (1-T);
% plot(t_edf, cdf, 'r--');
v = cdf - f_edf;