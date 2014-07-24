function [mu, mu_std, A, pdf] = fitExpToHist(ti, ni)


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

mu0 = nansum(ni.*ti)/nansum(ni);

[p,resnorm,~,~,~,~,J] = lsqnonlin(@costFct, [mu0 1], [], [], opts, ti, ni);
mu = p(1);
A = p(2);

J = full(J);
mu_std = sqrt( resnorm/(numel(ni)-2) * inv(J'*J) ); %#ok<MINV>

pdf = A * 1/mu * exp(-1/mu*ti);
pdf(isnan(ni)) = NaN;


function v = costFct(p, ti, ni)
mu = p(1);
A = p(2);

pdf = A * 1/mu * exp(-1/mu*ti);
v = pdf-ni;
v(isnan(v)) = [];