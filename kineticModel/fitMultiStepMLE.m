%[k, y, BIC] = fitMultiStep(x, f, varargin) implements a multi-step process 
% of the form 
%
%    k1     k2      kn
% S0 --> S1 --> ... --> Sn
%
% FFT-based implementation


% Francois Aguet, 03/12/2013

function [k, k_pstd, y, i, BIC, C] = fitMultiStepMLE(x, data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addOptional('StepRange', 1:5);
% ip.addParamValue('MaxSteps', 6);
ip.addParamValue('Display', false, @islogical);
ip.parse(x, data, varargin{:});
% K = ip.Results.MaxSteps;
stepRange = ip.Results.StepRange;

dx = x(2)-x(1);
% initial value
mu = sum(data.*x*dx)/sum(data*dx);
k0 = 1/mu;

K = numel(stepRange);

n = numel(data);
RSS = NaN(1,K);
BIC = NaN(1,K);
kvec = cell(1,K);
J = cell(1,K);

kref = [ 0.0347    0.1743    0.2666    0.0376    0.1102 0.1]; % temp fix
% k0 = kref(1:stepRange(1));
% res = cell(1,K);
% res = struct([]);
% res = [];
p = cell(1,K);
for i = 1:K

    %     k0 = 0.15+zeros(i,1);
    %k0 = kref(1:i);
    k0 = rand(1,i)*0.05;
    [p{i}, res(i)] = fitHistMLE(x, data, @stepModelFFT, k0);
    k = p{i}(:,1);
   
    % update initial values
    k0 = [k; max(k)/2];
end

% best fit
BIC = vertcat(res.BIC);

% pairwise t-test between successive BIC values
hval = zeros(1,K-1);
pval = zeros(1,K-1);
for i = 1:K-1
    [hval(i), pval(i)] = ttest2Parameters(BIC(i:i+1,1), BIC(i:i+1,2), numel(data)*[1 1]);
end

% last significant difference:
% [~,i] = min(BIC(:,1));
%idx = find(hval==1, 1, 'last')+1;
idx = find(hval==0, 1, 'first');
k = p{idx}(:,1);
k_pstd = p{idx}(:,2);



% m = stepModelFFT(k, x);
% y = stepModelODE(k, x);
% figure; hold on; plot(x, y, 'r'), plot(x, m, 'g--')

%%
if ip.Results.Display
    ha = setupFigure(1,3, 'XSpace', [2 2 0.5], 'YSpace', [1.5 0.5 1]);
    plot(ha(1), x, data, 'k', 'LineWidth', 1);
    plot(ha(1), x, res(idx).model(:,1), 'r', 'LineWidth', 1);
    xlabel(ha(1), 'Lifetime (s)');
    ylabel(ha(1), '# Events');
    
    % rates
    he = errorbar(ha(2), k, k_pstd, 'LineStyle', 'none',...
        'Color', 0.6*[1 1 1], 'LineWidth', 2);
    setErrorbarStyle(he, 0);
    plot(ha(2), 1:idx, k, 'k.', 'MarkerSize', 12);
    set(ha(2), 'XLim', [0.5 idx+0.5], 'XTick', 1:3);
    YLim = get(ha(2), 'YLim');
    set(ha(2), 'YLim', [0 YLim(2)]);
    xlabel(ha(2), 'Step #');
    ylabel(ha(2), 'Rate (s^{-1})');
    
    % BIC
    he = errorbar(ha(3), BIC(:,1), BIC(:,2), 'LineStyle', 'none',...
        'Color', 0.6*[1 1 1], 'LineWidth', 2);
    setErrorbarStyle(he, 0);
    plot(ha(3), BIC(:,1), 'k.', 'MarkerSize', 12);
    set(ha(3), 'XLim', [0.5 K+0.5], 'XTick', 1:K);
    % indicate significant difference
    YLim = get(ha(3), 'YLim');
    i = idx-1;
    plot(ha(3), i+[0.05 0.05 0.95 0.95],...
        max(BIC(i:i+1,1)+BIC(i:i+1,2))+diff(YLim)*[0.04 0.08 0.08 0.04],...
        'k', 'LineWidth', 1);
    text(i+0.5, BIC(i,1)+BIC(i,2)+diff(YLim)*0.1, ['p < 10^{' num2str(ceil(log10(pval(2)))) '}'], 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center', 'Parent', ha(3));
    % indicate n.s. differences
    for i = idx:K-1
        if hval(i)==0
            plot(ha(3), i+[0.05 0.05 0.95 0.95],...
                max(BIC(i:i+1,1)+BIC(i:i+1,2))+diff(YLim)*[0.04 0.08 0.08 0.04],...
                'k', 'LineWidth', 1);
            text(i+0.5, BIC(i,1)+BIC(i,2)+diff(YLim)*0.1, 'n.s.', 'VerticalAlignment', 'bottom',...
                'HorizontalAlignment', 'center', 'Parent', ha(3));
        end
    end
    xlabel(ha(3), 'Model #');
    ylabel(ha(3), 'BIC');
    
    %K = corrMatFromCov(res(idx).C);
    %plotCorrelationMatrix(K, 'TickLabels', arrayfun(@(i) ['k_' num2str(i)], 1:i, 'unif', 0));

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
F = ones(size(x));%/dx;
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
m = m/sum(m);%/(x(2)-x(1));

function dy = ksteps(~, y, k)
S = -diag([k 0]) + diag(k,-1);
dy = S*y;
