%[k, y, BIC] = fitMultiStep(x, f, varargin) implements a multi-step process of the form S1 --> S2 --> ... --> Sn
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn

% Francois Aguet, 03/12/2013

function [k, y, BIC] = fitMultiStep(x, f, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('MaxSteps', 6);
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});
K = ip.Results.MaxSteps;


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

dx = x(2)-x(1);
% initial value
mu = sum(f.*x*dx)/sum(f*dx);
k0 = 1/mu;

n = numel(f);
BIC = zeros(1,K);
k_est = cell(1,K);
for i = 1:K
    
    [k,RSS] = lsqnonlin(@multistepCost, k0, zeros(1,i), [], opts, x, f);
    BIC(i) = n*log(RSS/n) + numel(k)*log(n);
    k_est{i} = k;   
    k0 = [k min(k)/2];
    if i>1 && BIC(i)>BIC(i-1);
        break
    end        
end
BIC = BIC(1:i);
k = k_est{i-1};

S0 = [1 zeros(1,numel(k))];
sol = ode45(@(t,y) ksteps(t, y, k), [0 x(end)], S0);
Y = deval(sol, x);
y = Y(end-1,:);
y = y/sum(y)/dx;

if ip.Results.Display
    figure;
    plot(1:i, BIC);
    xlabel('# steps');
    ylabel('BIC');
end

function v = multistepCost(k, x, f)
N = numel(x);
dx = x(2)-x(1);
w = ((0:N-1)-floor(N/2))/dx/N*2*pi;
F = ones(size(x))/dx;
for i = 1:numel(k)
    F = F.*k(i)./(k(i)+1j*w);
end
g = abs(ifft(ifftshift(F)));
v = g - f;

%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
function dy = ksteps(~, y, k)
S = -diag([k 0]) + diag(k,-1);
dy = S*y;
