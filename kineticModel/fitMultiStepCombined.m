%[k, y, BIC] = fitMultiStepCombined(x, f, varargin) implements a multi-step process
% of the form S1 --> S2 --> ... --> Sn to two data sets representing distinct exit states
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
%
% FFT-based implementation

% Francois Aguet, 08/04/2013

function [k, k_pstd, y, BIC, C] = fitMultiStepCombined(x, f, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-10, ...
    'Tolfun', 1e-10);

dx = x(2)-x(1);
% initial value
mu = sum(f{1}.*x*dx)/sum(f{1}*dx);

% # data points
n = sum(cellfun(@numel, f));

% config structure: #common steps, total #step model 1, total #steps model 2
% config = {[0 1 1], [1 2 2], [1 2 3], [1 2 4], [1 3 3], [1 3 4],...
%     [2 3 3], [2 3 4], [3 4 4], [3 4 5]};%, [2 4 5]};

config = {[2 3 4]};

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


% best fit
[~,i] = min(BIC);
k = kvec{i};
nk = numel(k);
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


if ip.Results.Display
    % 1) Plot data w/ best model
    ah = 3.5*1.5;
    aw = 6*1.5;
    setupFigure('AxesHeight', ah, 'AxesWidth', aw);
    fset = loadFigureSettings('print');
    %figure(fset.fOpts{:});
    %axes(fset.axOpts{:});
    plot(x, f{1}, 'LineWidth', 1, 'Color', hsv2rgb([0.6 0.3 0.9]), 'LineWidth', 1);
    plot(x, y1, 'LineWidth', 1, 'Color', hsv2rgb([0.6 1 0.9]), 'LineWidth', 1.25);
    plot(x, f{2}, 'LineWidth', 1, 'Color', hsv2rgb([1/3 0.3 0.9]), 'LineWidth', 1);
    plot(x, y2, 'LineWidth', 1, 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 1.25);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    axis([0 160 0 0.05]);
    set(gca, 'XTick', 0:20:160, 'TickLength', fset.TickLength/6*6);
    
    % BIC inset
    axes(fset.axOpts{:}, 'Units', 'Normalized', 'Position', [4.5/8 3/5.5 3/8 2/5.5],...
        'TickLength', fset.TickLength/3*6);
    hold on;
    plot(1:nc, BIC, '.', 'Color', hsv2rgb([0 0.3 0.9]), 'MarkerSize', 11);
    plot(i, BIC(i), '.', 'Color', hsv2rgb([0 1 0.9]), 'MarkerSize', 12);
    set(gca, 'XLim', [0.5 nc+.5], 'XTick', 1:nc);
    xlabel('Model #', fset.sfont{:});
    ylabel('BIC', fset.sfont{:});
    
    % residuals
%     %%
%     setupFigure('AxesHeight', ah/3, 'AxesWidth', aw);
%     set(gca, 'TickLength', fset.TickLength/2*2.5);
%     plot([0 160], [0 0], 'k--')
%     plot(x, f{1}-y1, 'Color', hsv2rgb([0.6 0.3 0.9]), 'LineWidth', 1);
%     plot(x, f{2}-y2, 'Color', hsv2rgb([1/3 0.3 0.9]), 'LineWidth', 1);
%     axis([0 160 -0.008 0.008]);
%     %%
    
    % 2) Correlation matrix
    klabel = arrayfun(@(i) ['k_' num2str(i)], 1:nk, 'unif', 0);
    plotCorrelationMatrix(K, 'TickLabels', klabel);
    
    
    % 3) Rates +/- s.d.
    xa = 1:nk;
    
    setupFigure();
    
    he = errorbar(k, k_pstd, 'Color', 0.4*[1 1 1], 'LineStyle', 'none', 'LineWidth', 1);
    setErrorbarStyle(he, 0.1);
    plot(k, 'k.', 'MarkerSize', 12);
    set(gca, 'XTick', xa, 'XTickLabel', ' ', 'YTick', 0:0.1:1, fset.sfont{:});
    axis([0.5 nk+0.5 0 0.5]);
    YLim = get(gca, 'YLim');
    
    arrayfun(@(i) text(xa(i), YLim(1)-0.05*diff(YLim), klabel(i), fset.sfont{:},...
        'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'Interpreter', 'TeX'), 1:nk, 'unif', 0);
    xlabel('Rate constants (s^{-1})', fset.lfont{:});
    %ylabel('k (s^{-1})', fset.lfont{:});
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
