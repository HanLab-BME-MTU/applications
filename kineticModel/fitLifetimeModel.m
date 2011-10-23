function res = fitLifetimeModel(t, lftDist)

% figure; plot(tvec, cumsum(hvec.*[0 diff(tvec)]))

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

% cf1 = [184 238 255]/255; % blue
% ce1 = [0 186 255]/255;
% 
% cf2 = [255 194 213]/255; % red
% ce2 = [223 8 0]/255;
% 
% cf3 = [1 1 1]*0.6;
% ce3 = [1 1 1]*0.3;
% 
% tfont = {'FontName', 'Helvetica', 'FontSize', 12};
sfont = {'FontName', 'Helvetica', 'FontSize', 20};
lfont = {'FontName', 'Helvetica', 'FontSize', 24};


dt = t(2)-t(1);
dti = dt/10;
t_fine = 0:dti:t(end);

n = numel(lftDist);

res = struct([]);
for i = 1:4

    % # states
    ns = i*2;
    
    % Intializations & bounds
    S0 = [1 zeros(1,ns-1)];
    k0 = 0.05 * ones(1,ns-1);
    lb = zeros(1,ns-1);
    ub = Inf(1,ns-1);

    [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(str2func(['@costPop' num2str(i) 'Model']), k0, lb, ub, opts, t, lftDist, S0);
    
    % BIC and correlation matrix
    J = full(J);
    C = resnorm/(n-numel(k)-1)*inv(J'*J);
    k_std = sqrt(diag(C));
    K = corrFromC(C);

    res(i).k = k;
    res(i).k_std = k_std;
    res(i).corr = K';
    res(i).BIC = n*log(resnorm/n) + numel(k)*log(n);
end

bicV = [res.BIC];
minIdx = find(bicV==min(bicV), 1, 'first');


% plot result of best fit
plotKineticModelRates(res(minIdx).k, res(minIdx).k_std, res(minIdx).corr);

% compute scaling factor/renormalize (or: normalize input to model??)
hf = str2func(['pop' num2str(minIdx) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, res(minIdx).k), [0 t(end)], S0);

modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);
nf = sum(dt*interp1(t_ode, modelPDF, t));

% interpolate at fine scale for display
Y = interp1(t_ode, Y, t_fine);
modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);


%------------------------------------
% Display result of best fit
%------------------------------------
figure('Position', [440 378 560 360], 'PaperPositionMode', 'auto');
hold on;
hp(1) = plot(t, lftDist, '.', 'MarkerSize', 20, 'Color', [0 0 0]);
hp(2) = plot(t_fine, Y(end,2)*Y(:,1)/nf, 'Color', [0 0.8 0], 'LineWidth', 2);
plot(t_fine, Y(end,4)*Y(:,3)/nf, 'Color', [0 0.8 0], 'LineWidth', 2);
hp(3) = plot(t_fine, modelPDF/nf, 'r--', 'LineWidth', 4);

axis([0 100 0 0.12]);
set(gca, 'LineWidth', 2, 'Layer', 'top', sfont{:});
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});

hl = legend(hp, 'Meas. lifetime', 'Pop. lifetimes', 'Model');
set(hl, 'Box', 'off');


%%
% Plot BIC

figure;
hold on;
plot(bicV, 'r.', 'MarkerSize', 20);
set(gca, 'LineWidth', 2, 'Layer', 'top', sfont{:}, 'XLim', [0.5 numel(bicV)+0.5], 'XTick', 1:numel(bicV));
xlabel('# populations', lfont{:});
ylabel('BIC', lfont{:});






function v = costPop1Model(kVect, t, lftDist, S_init)
[ti, Y] = ode45(@(t,y) pop1Model(t, y, kVect), [0 t(end)], S_init);
modelPDF = Y(end,2)*Y(:,1);

% interpolate ODE output to input grid
modelPDF = interp1(ti, modelPDF, t); % same length as input
% renormalize (or: normalize input to model??)
dt = t(2)-t(1);
modelPDF = modelPDF/sum(dt*modelPDF);
v = modelPDF - lftDist;


function v = costPop2Model(kVect, t, lftDist, S_init)
[ti, Y] = ode45(@(t,y) pop2Model(t, y, kVect), [0 t(end)], S_init);
modelPDF = Y(end,2)*Y(:,1) + Y(end,4)*Y(:,3);

% interpolate ODE output to input grid
modelPDF = interp1(ti, modelPDF, t); % same length as input
% renormalize (or: normalize input to model??)
dt = t(2)-t(1);
modelPDF = modelPDF/sum(dt*modelPDF);
v = modelPDF - lftDist;


function v = costPop3Model(kVect, t, lftDist, S_init)
[ti, Y] = ode45(@(t,y) pop3Model(t, y, kVect), [0 t(end)], S_init);
modelPDF = Y(end,2)*Y(:,1) + Y(end,4)*Y(:,3) + Y(end,6)*Y(:,5);

% interpolate ODE output to input grid
modelPDF = interp1(ti, modelPDF, t); % same length as input
% renormalize (or: normalize input to model??)
dt = t(2)-t(1);
modelPDF = modelPDF/sum(dt*modelPDF);
v = modelPDF - lftDist;


function v = costPop4Model(kVect, t, lftDist, S_init)
[ti, Y] = ode45(@(t,y) pop4Model(t, y, kVect), [0 t(end)], S_init);
modelPDF = Y(end,2)*Y(:,1) + Y(end,4)*Y(:,3) + Y(end,6)*Y(:,5) + Y(end,8)*Y(:,7);

% interpolate ODE output to input grid
modelPDF = interp1(ti, modelPDF, t); % same length as input
% renormalize (or: normalize input to model??)
dt = t(2)-t(1);
modelPDF = modelPDF/sum(dt*modelPDF);
v = modelPDF - lftDist;


% Model:
%       k1
%    S1 --> S2
function dy = pop1Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -k(1)*y(1);
dy(2) = k(1)*y(1);


% Model:
%       k2     k3 
%    S1 --> S3 --> S4
% k1 |       
%    v
%    S2
function dy = pop2Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -k(3)*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);


% Model:
%       k2     k4     k5
%    S1 --> S3 --> S5 --> S6
% k1 |   k3 | 
%    v      v
%    S2     S4
function dy = pop3Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -(k(3)+k(4))*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);
dy(5) = -k(5)*y(5) + k(4)*y(3);
dy(6) = k(5)*y(5);


% Model:
%       k2     k4     k6     k7
%    S1 --> S3 --> S5 --> S7 --> S8
% k1 |   k3 |   k5 |
%    v      v      V
%    S2     S4     S6
function dy = pop4Model(~, y, k)
N = numel(y);
dy = zeros(N,1);
dy(1) = -(k(1)+k(2))*y(1);
dy(2) = k(1)*y(1);
dy(3) = -(k(3)+k(4))*y(3) + k(2)*y(1);
dy(4) = k(3)*y(3);
dy(5) = -(k(5)+k(6))*y(5) + k(4)*y(3);
dy(6) = k(5)*y(5);
dy(7) = -k(7)*y(7) + k(6)*y(5);
dy(8) = k(7)*y(7);


function K = corrFromC(C)
n = size(C,1);
K = zeros(n,n);

idx = pcombs(1:n);
i = idx(:,1);
j = idx(:,2);
ij = i+n*(j-1);
ii = i+n*(i-1);
jj = j+n*(j-1);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));