%[k, y, BIC] = fitMultiStep(x, f, varargin) implements a multi-step process of the form S1 --> S2 --> ... --> Sn
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
%
% FFT-based implementation


% Francois Aguet, 03/12/2013

function [k, k_pstd, y, i, BIC, C] = fitMultiStep(x, f, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addOptional('StepRange', 1:5);
% ip.addParamValue('MaxSteps', 6);
ip.addParamValue('Display', false, @islogical);
ip.parse(x, f, varargin{:});
% K = ip.Results.MaxSteps;
stepRange = ip.Results.StepRange;

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-10, ...
    'Tolfun', 1e-10);

dx = x(2)-x(1);
% initial value
mu = sum(f.*x*dx)/sum(f*dx);
k0 = 1/mu;

K = numel(stepRange);

n = numel(f);
RSS = NaN(1,K);
BIC = NaN(1,K);
kvec = cell(1,K);
J = cell(1,K);

% kref = [ 0.0347    0.1743    0.2666    0.0376    0.1102 0.1]; % temp fix
% k0 = kref(1:stepRange(1));
for i = 1:K
    %lb = zeros(1,i);
    lb = zeros(1,numel(k0));
    [k,RSS(i),~,~,~,~,J{i}] = lsqnonlin(@multistepCost, k0, lb, [], opts, x, f);
    BIC(i) = n*log(RSS(i)/n) + numel(k)*log(n);
    kvec{i} = k;
    k0 = [k max(k)/2];
%     k0 = [k min(k)*2];
    %k0 = rand(1,numel(k)+1);
    %k0 = 0.01*ones(1,numel(k)+1);
%     k0 = kref(1:numel(k)+1);
end

% best fit
[~,i] = min(BIC);
k = kvec{i};

% correlation btw. parameters of best fit
J = full(J{i});
C = RSS(i)/(numel(f)-numel(k)-1)*inv(J'*J); %#ok<MINV>
k_pstd = sqrt(diag(C))';
K = corrMatFromCov(C)';

% m = stepModelFFT(k, x);
y = stepModelODE(k, x);
% figure; hold on; plot(x, y, 'r'), plot(x, m, 'g--')


if ip.Results.Display
    plotCorrelationMatrix(K, 'TickLabels', arrayfun(@(i) ['k_' num2str(i)], 1:i, 'unif', 0));
    
    setupFigure();
    plot(k, 'k.');
    errorbar(k, k_pstd)
    
    figure;
    plot(1:i, BIC);
    xlabel('# steps');
    ylabel('BIC');
end


function v = multistepCost(k, x, f)
dx = x(2)-x(1);
xi = 0:dx:x(end);
mi = stepModelFFT(k,xi);
% mi = stepModelODE(k,xi);
m = interp1(xi,mi,x);

v = (m/sum(m) - f/sum(f))/dx;


% Models of the form:
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn

% FFT-based implementation
function m = stepModelFFT(k, x)
N = numel(x);
dx = x(2)-x(1);
w = ((0:N-1)-floor(N/2))/dx/N*2*pi;
F = ones(size(x))/dx;
for i = 1:numel(k)
    F = F.*k(i)./(k(i)+1j*w);
end
m = abs(ifft(ifftshift(F)));

% ODE solver-based implementation
function m = stepModelODE(k, x)
S0 = [1 zeros(1,numel(k))];
sol = ode45(@(t,y) ksteps(t, y, k), [0 x(end)], S0);
Y = deval(sol, x);
m = Y(end-1,:);
m = m/sum(m)/(x(2)-x(1));

function dy = ksteps(~, y, k)
S = -diag([k 0]) + diag(k,-1);
dy = S*y;
