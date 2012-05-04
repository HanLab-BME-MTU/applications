function [mu sigma x g] = fitGaussianModeToHist(xi, ni, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

dxi = xi(2)-xi(1);
A0 = max(ni);
mu0 = sum(ni.*xi)*dxi;
sigma0 = sqrt(sum(ni.*xi.^2)*dxi - mu0^2);

[p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, [A0 mu0 sigma0], [0 0], [Inf Inf], opts, xi, ni);



A = p(1);
mu = p(2);
sigma = p(3);

% A = max(ni);
% mu = p(1);
% sigma = p(2);

x = 0:dxi/10:xi(end);
g = A * exp(-(x-mu).^2/(2*sigma^2));

if ip.Results.Display
    figure;
    hold on;
    plot(xi, ni, 'k.-');
    plot(x, g, 'r');
end



function v = cost(p, xi, ni)
A = p(1);
mu = p(2);
sigma = p(3);
% A = max(ni);
% mu = p(1);
% sigma = p(2);

g = A * exp(-(xi-mu).^2/(2*sigma^2));
v = g-ni;
v(xi>mu+sigma) = 0;

