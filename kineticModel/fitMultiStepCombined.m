%[k, y, BIC] = fitMultiStep(x, f, varargin) implements a multi-step process of the form S1 --> S2 --> ... --> Sn
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
%
% FFT-based implementation


% Francois Aguet, 03/12/2013

function [k, y, BIC, C] = fitMultiStepCombined(x, f, varargin)

ip = inputParser;
ip.CaseSensitive = false;
% ip.addParamValue('MaxSteps', 6);
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});
% K = ip.Results.MaxSteps;


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-10, ...
    'Tolfun', 1e-10);

dx = x(2)-x(1);
% initial value
mu = sum(f{1}.*x*dx)/sum(f{1}*dx);

% % config = [2 3 4];
% % config = [2 3 3];
% % config = [1 3 4];
% config = [1 2 3];

% # data points
n = sum(cellfun(@numel, f));

% config structure: #common steps, total #step model 1, total #steps model 2

% config = {[1 2 3], [1 2 4], [1 3 4], [1 3 3],...
%     [2 3 3], [2 3 4], [3 4 4], [3 4 5]};

config = {[0 1 1], [1 2 2], [1 2 3], [1 2 4], [1 3 4], [1 3 3],...
    [2 3 3], [2 3 4], [3 4 4], [3 4 5]};

nc = numel(config);
RSS = NaN(1,nc);
BIC = NaN(1,nc);
kvec = cell(1,nc);
J = cell(1,nc);

% loop through configurations
for i = 1:nc
    
    % # rate params to estimate:
    nk = config{i}(2) + config{i}(3) - config{i}(1);
    k0 = nk/mu*ones(1,nk);

    [k, RSS(i), ~, ~, ~, ~, J{i}] = lsqnonlin(@multistepCost, k0, zeros(1,nk), [], opts, config{i}, x, f);
    BIC(i) = n*log(RSS(i)/n) + numel(k)*log(n);
    kvec{i} = k;
end

%%
% best fit
[~,i] = min(BIC);
k = kvec{i};
k1 = k(1:config{i}(2));
k2 = [k(1:config{i}(1)) k(config{i}(2)+1:end)];

y1 = stepModelODE(k1, x);
y2 = stepModelODE(k2, x);
y = {y1,y2}; % fct output


% correlation btw. parameters of best fit
J = full(J{i});
C = RSS(i)/(n-numel(k)-1)*inv(J'*J); %#ok<MINV>
k_pstd = sqrt(diag(C))';
K = corrMatFromCov(C)';


% J1 = J(1:120,1:3);
% C1 = RSS(i)/(n-numel(k)-1)*inv(J1'*J1); %#ok<MINV>
% k1_pstd = sqrt(diag(C1))';
% 
% J2 = J(121:end,[1 2 4 5]);
% C2 = RSS(i)/(n-numel(k)-1)*inv(J2'*J2); %#ok<MINV>
% k2_pstd = sqrt(diag(C2))';


k
k_pstd

%%
if ip.Results.Display
    setupFigure();%'AxesHeight', 5.25, 'AxesWidth', 9);
    fset = loadFigureSettings('print');
    plot(x, f{1}, 'LineWidth', 1, 'Color', hsv2rgb([0.6 0.3 0.9]));
    plot(x, y1, 'LineWidth', 1, 'Color', hsv2rgb([0.6 1 0.9]));
    plot(x, f{2}, 'LineWidth', 1, 'Color', hsv2rgb([1/3 0.3 0.9]));
    plot(x, y2, 'LineWidth', 1, 'Color', hsv2rgb([1/3 1 0.9]));
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    axis([0 120 0 0.05]);
    
    % BIC inset
    axes(fset.axOpts{:}, 'Units', 'Normalized', 'Position', [5.5/8 3/5.5 2/8 2/5.5], 'TickLength', fset.TickLength/2*6);
    hold on;
    plot(1:nc, BIC, '.', 'Color', hsv2rgb([0 0.3 0.9]), 'MarkerSize', 9);
    plot(i, BIC(i), '.', 'Color', hsv2rgb([0 1 0.9]), 'MarkerSize', 10);
    set(gca, 'XLim', [0.5 nc+.5], 'XTick', 1:nc);
    xlabel('Model #', fset.sfont{:});
    ylabel('BIC', fset.sfont{:});
   
    
    plotCorrelationMatrix(K, 'TickLabels', arrayfun(@(i) ['k_' num2str(i)], 1:numel(k), 'unif', 0));
    
    
    setupFigure();
    plot(k, 'k.');
    errorbar(k, k_pstd)

end


% k structure: [model 1, rem. of model 2]
function v = multistepCost(k, config, x, f)
dx = x(2)-x(1);
xi = 0:dx:x(end);

k1 = k(1:config(2));
k2 = [k(1:config(1)) k(config(2)+1:end)];

m1 = interp1(xi, stepModelFFT(k1,xi), x);
m2 = interp1(xi, stepModelFFT(k2,xi), x);

v = [(m1/sum(m1) - f{1}/sum(f{1}))/dx (m2/sum(m2) - f{2}/sum(f{2}))/dx];


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
