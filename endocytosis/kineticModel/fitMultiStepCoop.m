% Francois Aguet, 06/02/2013

function [k, m, t_fine, model, BIC, K] = fitMultiStepCoop(x, f, ns, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addRequired('ns');
ip.addOptional('SelV', []);
ip.addParamValue('Display', false, @islogical);
ip.parse(x, f, ns, varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

% initial values
dx = x(2)-x(1);
mu = sum(f.*x*dx)/sum(f*dx);
k0 = 1/mu;

dxi = dx/10;

f = f/sum(f)/dx;

% configurations to test
selV = ip.Results.SelV;
if isempty(selV)
    switch ns
        case 2
            selV = {[0 0], [1 0], [0 1], [1 1]};
            %selV = {[0 0]};
        case 3
            selV = {[0 0 0], [1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1], [1 1 1]};
            %selV = {[1 1 0]};
    end
end

nc = numel(selV);

k0 = [0.03 0.03 0.11];

BIC = zeros(1,nc);
RSS = zeros(1,nc);
kvec = cell(1,nc);
mvec = cell(1,nc);
for i = 1:nc
    sel = selV{i};
    lb = [zeros(1,ns) zeros(1,sum(sel))];
    
    m0 = 1.5*ones(1,ns);
    %init = [linspace(k0/2,k0*2,ns) m0(sel==1)];
    init = [k0(1:ns) m0(sel==1)];

    [p,RSS(i),~,~,~,~,J] = lsqnonlin(@multistepCost, init, lb, [], opts, x, f, sel);
    p = abs(p);
    kvec{i} = p(1:ns);
    m = ones(1,ns);
    m(sel==1) = p(ns+1:end);
    mvec{i} = m;
    
    n = numel(f);
    BIC(i) = n*log(RSS(i)/n) + numel(p)*log(n);
end

% best fit
[~,i] = min(BIC);
k = kvec{i};
m = mvec{i};

% sol = ode45(@(t,y) multiStepCoop(t, y, p(1:ns), p(ns+1:end)), [0 x(end)], S0);
% Y = deval(sol, t_fine);
% model = Y(end-1,:);
% model = model/sum(model)/dxi*2;

% correlation btw. parameters of best fit
J = full(J);
C = RSS(i)/(numel(f)-numel(k)-1)*inv(J'*J); %#ok<MINV>
prmStd = sqrt(diag(C))';
k_pstd = prmStd(1:ns);
m_pstd = prmStd(ns+1:end);

K = corrMatFromCov(C)';


t_fine = x(1):dxi:x(end);
model = multiStepModelCoop(t_fine, k, m);


%
% figure(fset.fOpts{:});
% axes(fset.axOpts{:});
% hold on;
% 
% plot(x, f, 'k', 'LineWidth', 1);
% plot(t_fine, model, 'r--', 'LineWidth', 1);
% 
% axis([0 120 0 0.05]);

% figure(fset.fOpts{:});
% axes(fset.axOpts{:});
% hold on;
% plot(t, mean(lftResDC3.lftHist_Aneg(:,idx),1), 'LineWidth', 1, 'Color', hsv2rgb([0.6 0.3 0.9]));
% plot(t, y1, 'LineWidth', 1, 'Color', hsv2rgb([0.6 1 0.9]));
% plot(t, mean(lftResDC3.lftHist_Apos(:,idx),1), 'LineWidth', 1, 'Color', hsv2rgb([1/3 0.3 0.9]));
% plot(t, y2, 'LineWidth', 1, 'Color', hsv2rgb([1/3 1 0.9]));
% xlabel('Lifetime (s)', fset.lfont{:});
% ylabel('Frequency', fset.lfont{:});
% axis([0 120 0 0.05]);


% k = 0.5;
% p = 3;
% S0 = [1 0]; % two states
% sol = ode45(@(t,y) multiStepCoop(t, y, k, p), [0 x(end)], S0);
% Y = deval(sol, t_fine);
% plot(t_fine, Y);


% fset = loadFigureSettings('');
% figure(fset.fOpts{:});
% axes(fset.axOpts{:});
% hold on;
% % h = linspace(0.55, 0.6, 3);
% h = linspace(0., 0.6, 3);
% v = ones(1,3);
% s = 0.6:0.2:1;
% cmap = hsv2rgb([h' s' v']);
% set(gca, 'ColorOrder', cmap);
% 
% dx = x(2)-x(1);
% dxi = dx/10;
% t_fine = 0:dxi:x(end);
% 
% 
% k = [0.1 0.1];
% p = [20 1];
% S0 = [1 0 0]; % two states
% sol = ode45(@(t,y) multiStepCoop(t, y, k, p), [0 x(end)], S0);
% % k = [0.03 0.95 2];
% % sol = ode45(@(t,y) multiStepCoop(t, y, k, [5 4]), [0 t(end)], S0);
% Y = deval(sol, t_fine);
% Y = (diag(1./sum(Y,2)))*Y; % normalize
% plot(t_fine, Y);
% set(gca, 'XLim', [0 120], 'YLim', [0 0.02]);
% % - VS - 
% k = [0.05 0.5];
% p = [1 5];
% S0 = [1 0 0]; % two states
% sol = ode45(@(t,y) multiStepCoop(t, y, k, p), [0 x(end)], S0);
% Y = deval(sol, t_fine);
% Y = (diag(1./sum(Y,2)))*Y; % normalize
% plot(t_fine, Y, '--');


% parameters: [rates, exponents]
function v = multistepCost(p, x, f, sel)
p = abs(p);
% split p into k and m
ns = numel(sel);
k = p(1:ns);
m = ones(1,ns);
m(sel==1) = p(ns+1:end);
v = multiStepModelCoop(x, k, m) - f;


function m = multiStepModelCoop(x, k, m)
S0 = [1 zeros(1,numel(k))];
sol = ode45(@(t,y) multiStepCoop(t, y, k, m), [0 x(end)], S0);
Y = deval(sol, x);
m = abs(Y(end-1,:));
m = m/sum(m)/(x(2)-x(1));

function dy = multiStepCoop(~, y, k, m)
M = getTransitionMatrix(k);
dy = M*y.^([m 1]');

% state transition matrix for multiple sequential steps
function M = getTransitionMatrix(k)
M = diag([-k 0])+diag(k,-1);
