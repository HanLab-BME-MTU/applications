%p = getGaussianConvPrms(x, pdf1, pdf2, varargin) determines the parameters of the Gaussian which convolved with pdf1 gives pdf2
%
% Example:
%
% sigma = 1;
% mu = 6;
% dx = 0.1;
% x = 0:dx:20;
% g = exp(-(x-mu).^2/(2*sigma^2));
% g = g/sum(g)/dx;
% g2 = conv(g,g, 'full')*dx;
% g2 = g2(1:numel(g));
% getGaussianConvPrms(x, g, g2, 'Display', true);

% Francois Aguet, 03/06/2012

function [p, y] = getGaussianConvPrms(x, pdf1, pdf2, varargin)

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

dx = x(2)-x(1);
pdf1 = pdf1/sum(pdf1);
pdf2 = pdf2/sum(pdf2);

% estimate shift
mu0 = sum((pdf2-pdf1).*x);

lb = [0 0.5]; % lower limits on mu, sigma
ub = [Inf Inf];
[p,resnorm,~,~,~,~,J] = lsqnonlin(@costX, [mu0 2], lb, ub, opts, x, pdf1, pdf2);

y = filterX(p(1), p(2), pdf1);
p = p*dx;
y = y/sum(y)/dx;

if ip.Results.Display
    
    figure;
    hold on;
    plot(x, pdf1, 'k');
    plot(x, pdf2, 'b-');
    plot(x, y, 'r--');
    legend('f_0(x)', 'f_1(x)', '(f_0 \ast g)(x)');
end


function v = costX(p, x, pdf1, pdf2)
y = filterX(p(1), p(2), pdf1);
v = y - pdf2;


function y = filterX(mu, sigma, f)
mu0 = floor(mu);

% first shift
f = zeroshift(f, mu0);

% then convolve with quasi-centered PSF
dm = mu-mu0;
w = ceil(4*sigma);
x = -w:w;
g = exp(-(x-dm).^2/(2*sigma^2));
g = g/sum(g);

y = conv(padarrayXT(f, [0 (numel(x)-1)/2]), g, 'valid');


% s: integer
function y = zeroshift(x, s)
if s>=0
    y = [zeros(1,s) x(1:end-s)];
else
    y = [x(1+s:end) zeros(1,s)];
end
