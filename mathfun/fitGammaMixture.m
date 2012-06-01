

function [prmVect, a, BIC] = fitGammaMixture(x, f, varargin)

opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', 1e5, ...
    'Display', 'off', ...
    'TolX', 1e-12, ...
    'Tolfun', 1e-12,...
    'Algorithm', 'active-set');

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addOptional('mode', 'kna', @ischar);
ip.addParamValue('FitMode', 'CDF', @(x) any(strcmpi(x, {'PDF','CDF'})));
ip.addParamValue('N', 3);
ip.parse(x, f, varargin{:});

N = ip.Results.N;
estIdx = regexpi('kna', ['[' ip.Results.mode ']']);
estIdx = arrayfun(@(i) 3*i+estIdx, 0:N-1, 'UniformOutput', false);
estIdx = [estIdx{:}];
% estIdx(end) = 0; % contribution of last component is constrained
% estIdx = estIdx(1:end-1);


% transform into EDF format for interpolation
if strcmpi(ip.Results.FitMode, 'CDF')
    [fu, idx] = unique(f, 'first');
    opts = optimset(opts, 'Jacobian', 'off');
else
    [fu, idx] = unique(cumsum(f)*(x(2)-x(1)), 'first');
end
% Mean of the Gamma distribution: n/k
mu0 = interp1(fu, x(idx), (1:N)/(N+1));


% initialize with shape parameters 1,2,..,2
n0 = 1:N;
% n0 = 2*n0;
k0 = n0./mu0;
A0 = ones(1,N)/N;

prmVect = reshape([k0; n0; A0], [3*N 1])';
lb = repmat([0 1 0], [1 N]);
ub = repmat([Inf Inf 1], [1 N]);

% [p,RSS,~,~,~,~,J] = lsqnonlin(@cost, prmVect(estIdx), lb(estIdx), ub(estIdx), opts, x, f, prmVect, estIdx, ip.Results.FitMode);
A = repmat([0 0 1], [1 N]); % equality constraint

switch ip.Results.FitMode
    case 'PDF'
        [p RSS] = fmincon(@(i) costPDF(i, x, f, prmVect, estIdx), prmVect(estIdx), [], [], A(estIdx), 1, lb(estIdx), ub(estIdx), [], opts);
        
    case 'CDF'
        [p RSS] = fmincon(@(i) costCDF(i, x, f, prmVect, estIdx), prmVect(estIdx), [], [], A(estIdx), 1, lb(estIdx), ub(estIdx), [], opts);
end
n = numel(f);
BIC = n*log(RSS/n) + numel(p)*log(n);

prmVect(estIdx) = p;
a = gammaMixture(x(1), prmVect, 'Mode', 'CDF');

% [w W] = gammaMixture(x, prmVect, 'Mode', 'CDF');
% figure;
% hold on;
% plot(x, f*(1-a)+a, 'k');
% plot(x, W, 'b');
% plot(x, w, 'r');

% [w W] = gammaMixture(x, prmVect, 'Mode', 'PDF');
% figure;
% hold on;
% plot(x, f, 'k');
% plot(x, W, 'b');
% plot(x, w, 'r');



function [v J] = costPDF(p, x, f, prmVect, estIdx)
prmVect(estIdx) = p;
% round shape parameters
% n = prmVect(2:3:end);
% n = round(n);
% n(n<1) = 1;
% prmVect(2:3:end) = n;

% missing data:
a = gammaMixture(x(1), prmVect, 'Mode', 'CDF');


[w,~,J] = gammaMixture(x, prmVect, 'Mode', 'PDF');
v = sum((w - (1-a)*f).^2);


function v = costCDF(p, x, f, prmVect, estIdx)
prmVect(estIdx) = p;
% missing data:
a = gammaMixture(x(1), prmVect, 'Mode', 'CDF');
a = 0;

w = gammaMixture(x, prmVect, 'Mode', 'CDF');
v = sum((w - ((1-a)*f+a)).^2);
