%res = fitLifetimeModel(lftData, varargin)
%
% Inputs:
%      lftData : structure returned by 'runLifetimeAnalysis'
%
% Options: 
%       'Mode' : 'PDF | {'CDF'} fit to the histogram or empirical distribution
%   'AlphaBIC' : Threshold for BIC selection; default: 0.95
%       'MaxP' : Maximum number of populations to fit. Default: 3

% Francois Aguet (last modified 10/25/2011)

function res = fitLifetimeModel(lftData, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftData');
ip.addParamValue('Mode', 'CDF', @(x) any(strcmpi(x, {'PDF', 'CDF'})));
ip.addParamValue('MaxP', 3, @(x) any(ismember(x, 1:4)));
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('AlphaBIC', 0.95);
ip.parse(lftData, varargin{:});
dBIC = 2*log(ip.Results.AlphaBIC/(1-ip.Results.AlphaBIC));


% cutoff at last non-zero data point
endIdx = find(lftData.meanHist~=0, 1, 'last');

t = lftData.t(1:endIdx);
lftHist = lftData.meanHist(1:endIdx);
lftECDF = lftData.meanECDF(1:endIdx);

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', 1e5, ...
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
tfont = {'FontName', 'Helvetica', 'FontSize', 16};
sfont = {'FontName', 'Helvetica', 'FontSize', 20};
lfont = {'FontName', 'Helvetica', 'FontSize', 24};


dt = t(2)-t(1);
dti = dt/10;
t_fine = 0:dti:t(end);
n = numel(lftHist);


for i = 1:ip.Results.MaxP

    % # states
    ns = i*2;
    
    % Intializations & bounds
    S0 = [1 zeros(1,ns-1)];
    k0 = 0.05 * ones(1,ns-1);
    lb = zeros(1,ns-1);
    ub = Inf(1,ns-1);

    if strcmpi(ip.Results.Mode, 'CDF')
        [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costCDF, k0, lb, ub, opts, t, lftECDF, S0, i);
    else
        [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costPDF, k0, lb, ub, opts, t, lftHist, S0, i);
    end
    
    % BIC and correlation matrix
    J = full(J);
    C = resnorm/(n-numel(k)-1)*inv(J'*J);
    k_std = sqrt(diag(C))';
    K = corrFromC(C)';

    res.k{i} = k;
    res.k_std{i} = k_std;
    res.corr{i} = K;
    res.BIC(i) = n*log(resnorm/n) + numel(k)*log(n);
end


% Only differences in BIC with probability ? alpha (0.95 default)
% are considered significant
sortBIC = sort(res.BIC);
minIdx = find(res.BIC==min(sortBIC(diff(sortBIC) > dBIC)));

ns = minIdx*2;
% Intializations & bounds
S0 = [1 zeros(1,ns-1)];


% compute scaling factor/renormalize (or: normalize input to model??)
hf = str2func(['pop' num2str(minIdx) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, res.k{minIdx}), [0 t(end)], S0);

modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);
nf = sum(dt*interp1(t_ode, modelPDF, t));

% interpolate at fine scale for display
Y_fine = interp1(t_ode, Y, t_fine);

popMat = bsxfun(@times, Y_fine(end,2:2:end), Y_fine(:,1:2:end)) / nf;
modelPDF = sum(popMat, 2);
np = size(popMat,2);

pECDF = arrayfun(@(i) Y_fine(:,i)/Y_fine(end,i), 2:2:ns, 'UniformOutput', false);
res.pPercentiles = arrayfun(@(i) interp1(pECDF{i}, t_fine, [0.05 0.25 0.5 0.75 0.95]), 1:np, 'UniformOutput', false);
res.pA = Y_fine(end,2:2:end);
res.pMean = arrayfun(@(i) sum(t_fine.*popMat(:,i)'*dti) / sum(popMat(:,i)*dti), 1:np);

%------------------------------------
% Display result of best fit
%------------------------------------
figure('Position', [440 378 560+150 360], 'PaperPositionMode', 'auto');
% axes('Position', [0.15 0.18 0.8 0.75]);
axes('Units', 'Pixels', 'Position', [85 65 450 270]);

hold on;
hp(1) = plot(t, lftHist, '.', 'MarkerSize', 20, 'Color', [0 0 0]);
hi = plot(t_fine, popMat, 'Color', [0 0.7 1], 'LineWidth', 2);
hp(2) = hi(1);
hp(3) = plot(t_fine, modelPDF, '--', 'Color', [0.8 0 0], 'LineWidth', 4);

axis([0 100 0 0.12]);
set(gca, 'LineWidth', 2, 'Layer', 'top', sfont{:});
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});

% hl = legend(hp, 'Meas. lifetime', 'Pop. lifetimes', 'Model');
% set(hl, 'Box', 'off');

% Insets with population percentiles
xlabels = arrayfun(@(i) ['P' num2str(i)], 1:np, 'UniformOutput', false);

ha = axes('Units', 'Pixels', 'Position', [85+450+60 270 110 65]);
barplot2(res.pA, 'AdjustFigure', false, 'XLabels', cell(1,np),...
    'EdgeColor', [0 0.7 1], 'Color', [0.7 0.9 1],...
    'BarWidth', 0.5, 'GroupDistance', 0.5);
set(ha, tfont{:}, 'YTick', 0:0.2:0.8, 'YLim', [0 0.8]);
ylabel('Contrib.', sfont{:})
pos = get(get(ha, 'YLabel'), 'Position');

dx = pos(1);% = -0.4;
% set(get(ha, 'YLabel'), 'Position', pos);




ha = axes('Units', 'Pixels', 'Position', [85+450+60 65 110 180]);
pct = vertcat(res.pPercentiles{:})';
M = [pct(3,:); pct(2,:); pct(4,:); pct(1,:); pct(5,:)];
boxplot2(M, 'AdjustFigure', false, 'XLabels', xlabels,...
    'BarWidth', 0.5, 'GroupDistance', 0.5);
set(ha, tfont{:});
ylabel('Lifetime (s)', sfont{:})
pos = get(get(ha, 'YLabel'), 'Position');
pos(1) = dx;
set(get(ha, 'YLabel'), 'Position', pos);




%======================================

% test CDF fit
% figure; plot(t, lftData.meanECDF(1:endIdx), 'k.-');
% hold on;

% compute ECDF on 't'
% Y = interp1(t_ode, Y, t);
% modelCDF = sum(Y(:,2:2:end),2);
% T = modelCDF(1);
% modelCDF = (modelCDF-T)/(1-T) * (lftData.meanECDF(endIdx) - lftData.meanECDF(1)) + lftData.meanECDF(1);

% equiv. to CDF:
% modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);
% mx = cumsum(modelPDF)*dt/nf;
% plot(t, mx, 'g--');

% plot(t, modelCDF, 'r.--');
% plot(lftData.t_ecdf{1}, lftData.f_ecdf{1}, 'b--');

% attempt fit
% k0 = [0.05 0.05 0.05];
% [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costCDF, k0, lb, ub, opts, t, lftData.meanECDF(1:endIdx), S0, i);

% [t_ode, Y] = ode45(@(t,y) hf(t, y, k), [0 t(end)], S0);
% Y = interp1(t_ode, Y, t);

% modelCDF = sum(Y(:,2:2:end),2);
% T = modelCDF(1);
% modelCDF = (modelCDF-T)/(1-T) * (lftData.meanECDF(endIdx) - lftData.meanECDF(1)) + lftData.meanECDF(1);
% plot(t, modelCDF, 'm');

% modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);
% 
% p1 = Y(end,2)*Y(:,1)/sum(modelPDF*dt);
% p2 = Y(end,4)*Y(:,3)/sum(modelPDF*dt);
% 
% modelPDF = modelPDF/sum(modelPDF*dt);



% Plot control metrics: BIC, correlation btw. rates
if ip.Results.PlotAll
    figure;
    hold on;
    plot(res.BIC, 'r.', 'MarkerSize', 20);
    set(gca, 'LineWidth', 2, 'Layer', 'top', sfont{:}, 'XLim', [0.5 numel(res.BIC)+0.5], 'XTick', 1:numel(res.BIC));
    xlabel('# populations', lfont{:});
    ylabel('BIC', lfont{:});
    
    plotKineticModelRates(res.k{minIdx}, res.k_std{minIdx}, res.corr{minIdx});
end



% M: model #
function v = costPDF(kVect, t, lftHist, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[ti, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);
modelPDF = sum(bsxfun(@times, Y(end,2:2:end), Y(:,1:2:end)), 2);

% interpolate ODE output to input grid
modelPDF = interp1(ti, modelPDF, t); % same length as input
% renormalize (or: normalize input to model??)
dt = t(2)-t(1);
modelPDF = modelPDF/sum(dt*modelPDF);
v = modelPDF - lftHist;




function v = costCDF(kVect, t, lftECDF, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);
modelCDF = sum(Y(:,2:2:end),2);

% interpolate ODE output to input grid
modelCDF = interp1(t_ode, modelCDF, t);

% normalize to [0..1]
T = modelCDF(1);

% modelCDF = (modelCDF - modelCDF(1)) / (modelCDF(end)-modelCDF(1));
modelCDF = (modelCDF-T)/(1-T) * (lftECDF(end) - lftECDF(1)) + lftECDF(1);

v = modelCDF - lftECDF;








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