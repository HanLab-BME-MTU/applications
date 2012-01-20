%res = fitLifetimeModel(lftData, varargin)
%
% Inputs:
%         lftData : structure returned by 'runLifetimeAnalysis'
%
% Options:
%          'Mode' : 'PDF | {'CDF'} fit to the histogram or empirical distribution
%          'NumP' : Number of populations to test/fit. Default: 1:3
%  'ConstrainBIC' : Smallest BIC with probability > AlphaBIC than next candidate is chosen
%      'AlphaBIC' : Threshold for BIC selection; default: 0.95
%       'PlotAll' : Display correlation matrix and BIC values
%       'PlotCDF' : Display CDF and fitted model
%
% Output:
%             res : output structure with fields:
%                .k            : vector of rates from the model
%                .k_std        : standard deviation (error propagated) of the rates
%                .corr         : correlation matrix for the rates
%                .BIC          : BIC for the populations tested
%                .pPercentiles : percentiles of each subpopulations in the optimal model
%                .pA           : constributions of each subpopulation in the optimal model
%                .pMean        : means of each subpopulation in the optimal model

% Francois Aguet (last modified 10/25/2011)

function res = fitLifetimeModel(lftData, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftData');
ip.addParamValue('Mode', 'PDF', @(x) any(strcmpi(x, {'PDF', 'CDF'})));
ip.addParamValue('NumP', 3, @(x) all(ismember(x, 1:4)));
ip.addParamValue('ConstrainBIC', true, @islogical);
ip.addParamValue('AlphaBIC', 0.95);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('PlotCDF', false, @islogical);
% ip.addParamValue('ShowInset', false, @islogical);
ip.addParamValue('EndIdx', find(lftData.meanHist~=0, 1, 'last'));
ip.parse(lftData, varargin{:});
dBIC = 2*log(ip.Results.AlphaBIC/(1-ip.Results.AlphaBIC));


% cutoff at last non-zero data point
endIdx = ip.Results.EndIdx;

t = lftData.t(1:endIdx);
dt = t(2)-t(1);

lftHist = lftData.meanHist(1:endIdx);
lftHist = lftHist/sum(lftHist)/dt;
lftECDF = lftData.meanECDF(1:endIdx);

a = lftECDF(1);
lftECDF = (lftECDF-a)/(1-a);

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', 1e5, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

fset = loadFigureSettings();

dti = dt/10;
t_fine = 0:dti:t(end);
n = numel(lftHist);

for i = ip.Results.NumP
    
    % # states
    ns = i*2;
    
    % Intializations & bounds
    S0 = [1 zeros(1,ns-1)];
    switch i
        case 1
            k0 = 0.01;
        case 2
            k0 = [0.2 0.2 0.05];
        case 3
            k0 = [0.2 0.2 0.05 0.01 0.02];
            %k0 = 0.01 * ones(1,ns-1);
        case 4
            %k0 = [0.2 0.2 0.05 0.01 0.02 0.01 0.02];
            k0 = 0.02 * ones(1,ns-1);
    end
    lb = zeros(1,ns-1);
    ub = Inf(1,ns-1);
    
    switch ip.Results.Mode
        case 'PDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costPDF, k0, lb, ub, opts, t, lftHist, S0, i);
        case 'CDF'
            [k, resnorm, ~, ~, ~, ~, J] = lsqnonlin(@costCDF, k0, lb, ub, opts, t, lftECDF, S0, i);
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


% Only significant differences in BIC (alpha = 0.05 default) are considered
if ip.Results.ConstrainBIC && numel(res.BIC)>1
    sortBIC = sort(res.BIC);
    minIdx = find(res.BIC==min(sortBIC(diff(sortBIC) > dBIC)));
else
    minIdx = find(res.BIC==min(res.BIC));
end

np = minIdx;

% Intializations & bounds
ns = np*2;
S0 = [1 zeros(1,ns-1)];


hf = str2func(['pop' num2str(np) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, res.k{np}), [0 t(end)], S0);

% interpolate result over full time vector
popPDF = zeros(np,numel(t_fine));
Y = interp1(t_ode, Y, t_fine);

% normalize subpopulation PDFs, weigh by output
for k = 1:np
    p = Y(:,2*k-1);
    popPDF(k,:) = p/sum(p)/dti * Y(end,2*k);
end
pdf = sum(popPDF, 1);

popCDF = Y(:,2:2:end);
cdf = sum(Y(:,2:2:end),2);
a = interp1(t_fine, cdf, t(1));

% Compute population percentiles, mean, and contribution
pECDF = arrayfun(@(i) Y(:,i)/Y(end,i), 2:2:ns, 'UniformOutput', false);
for p = 1:np
    [u, uidx] = unique(pECDF{p});
    res.pPercentiles{p} = interp1(u, t_fine(uidx), [0.05 0.25 0.5 0.75 0.95]);
end
res.pA = Y(end,2:2:end);
res.pMean = arrayfun(@(i) sum(t_fine.*popPDF(i,:)*dti) / sum(popPDF(i,:)*dti), 1:np);

%------------------------------------
% Display result of best fit
%------------------------------------
switch np
    case 1
        colorOrder = fset.ceG;
    case 2
        colorOrder = [fset.ceR; fset.ceG];
    case 3
        colorOrder = [fset.ceR; fset.ceB2; fset.ceG];
    case 4
        colorOrder = [fset.ceR; fset.ceB2; fset.ceB; fset.ceG];
end
colorOrderFill = rgb2hsv(colorOrder);
colorOrderFill(:,2) = colorOrderFill(:,2)*0.3;
colorOrderFill = hsv2rgb(colorOrderFill);


dx = 85; % spacing between panels
% layout: width: 85 450 dx 200 dx 240

if ip.Results.PlotAll
    xt = 250;
else
    xt = 0;
end

figure('Position', [240 378 850+xt 360], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
%---------------------------------
% Lifetime histogram
%---------------------------------
axes('Units', 'Pixels', 'Position', [85 65 450 270]);
set(gca, 'ColorOrder', colorOrder);
hold on;
hp(1) = plot(t, lftHist*(1-a), '.', 'MarkerSize', 20, 'Color', [0 0 0]);
hi = plot(t_fine, popPDF, 'LineWidth', 2);
hp(2) = hi(1);
hp(3) = plot(t_fine, pdf, '--', 'Color', fset.ceB, 'LineWidth', 4);

axis([0 100 0 getYAxisBound(max(lftHist)*(1-a))]);
set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:});
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% hl = legend(hp, 'Meas. lifetime', 'Pop. lifetimes', 'Model');
% set(hl, 'Box', 'off');

%---------------------------------
% Inset with amplitudes
%---------------------------------


% main axes: [85 65 450 270]
% ha = axes('Units', 'Pixels', 'Position', [85+450-110 270 110 65]); % matches edge
ha = axes('Units', 'Pixels', 'Position', [85+450-32*np-10 260 32*np 65]);
% ha = axes('Units', 'Pixels', 'Position', [85+450+60 270 110 65]);
xlabels = arrayfun(@(i) ['P' num2str(i)], 1:np, 'UniformOutput', false);

barplot2(res.pA, 'AdjustFigure', false, 'XLabels', xlabels,...
    'FaceColor', colorOrderFill, 'EdgeColor', colorOrder,...
    'BarWidth', 0.6, 'GroupDistance', 0.5, 'Angle', 45);
set(ha, fset.tfont{:}, 'YTick', 0:0.2:0.8, 'YLim', [0 0.8]);
ylabel('Contrib.', fset.sfont{:})

% pos = get(get(ha, 'YLabel'), 'Position');
% % dx = pos(1);

%---------------------------------
% Lifetime histogram zoom
%---------------------------------
% if ip.Results.ShowInset
%     % Inset with zoom
%     axes('Units', 'Pixels', 'Position', [300 200 220 120]);
%     set(gca, 'ColorOrder', colorOrder);
%     
%     hold on;
%     hp(1) = plot(t, lftHist*(1-a), '.', 'MarkerSize', 20, 'Color', [0 0 0]);
%     hi = plot(t_fine, popMat, 'LineWidth', 2);
%     hp(2) = hi(1);
%     hp(3) = plot(t_fine, pdf, '--', 'Color', fset.ceB, 'LineWidth', 4);
%     axis([10 40 0.005 0.035]);
%     set(gca, 'LineWidth', 1.5, 'Layer', 'top', fset.tfont{:});
% end

ha = axes('Units', 'Pixels', 'Position', [85+450+dx 65 200 270]);
rateLabels = plotKineticModelRates(res.k{np}, res.k_std{np}, 'Handle', ha);

if ip.Results.PlotAll
 
    % 85 450 dx 200 dx 240
    nk = numel(res.k{np})-1;
    ha = axes('Units', 'Pixels', 'Position', [85+450+200+dx+60 65+270-30*nk 30*nk 30*nk]);
    plotCorrelationMatrix(res.corr{np}, 'Handle', ha, 'TickLabels', rateLabels, 'ColorBar', false);
    axis off;
    
    axes('Units', 'Pixels', 'Position', [85+450+200+dx+60+30*nk+10 65+270-30*4 1 30*4]);
    axis off;
    
    values = -1:1/100:1;
    N = length(values);
    map = zeros(N,3);
    
    ridx = values<0;
    map(ridx,1) = -values(ridx);
    gidx = values>0;
    map(gidx,2) = values(gidx);
    colormap(map);
    caxis([-1 1]);
    hc = colorbar('Units', 'pixels', 'YTick', -1:0.2:1);
    pos = get(hc, 'Position');
    pos(3) = 15;
    set(hc, 'Position', pos);
    
    % Plot control metrics: BIC, correlation btw. rates
    if numel(ip.Results.NumP)>1
        figure('Position', [440 378 400 300], 'PaperPositionMode', 'auto');
        np = numel(res.BIC);
        %axes('Units', 'pixels', 'Position', [85+450+200+dx+60 65 60*np 120]);
        axes('Units', 'pixels', 'Position', [100 60 60*np 220])
        
        hold on;
        plot(res.BIC, 'r.', 'MarkerSize', 40);
        set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:}, 'XLim', [0.5 np+0.5], 'XTick', 1:np);
        xlabel('# populations', fset.lfont{:});
        ylabel('BIC', fset.lfont{:});
    end
end

%---------------------------------
% Lifetime CDF
%---------------------------------
if ip.Results.PlotCDF
    figure('Position', [440 378 550 360], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
    axes('Units', 'Pixels', 'Position', [85 65 450 270]);
    set(gca, 'ColorOrder', colorOrder);
    hold on;
    plot(t, lftECDF*(1-a)+a, 'k.', 'MarkerSize', 20)
    plot(t_fine, popCDF, 'LineWidth', 2);
    plot(t_fine, cdf, '--', 'Color', fset.ceB, 'LineWidth', 4);
    
    axis([0 100 0 1]);
    set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
end

% ha = axes('Units', 'Pixels', 'Position', [85+450+60 65 110 180]);
% pct = vertcat(res.pPercentiles{:})';
% M = [pct(3,:); pct(2,:); pct(4,:); pct(1,:); pct(5,:)];
% boxplot2(M, 'AdjustFigure', false, 'XLabels', xlabels,...
%     'FaceColor', colorOrderFill, 'EdgeColor', colorOrder,...
%     'BarWidth', 0.5, 'GroupDistance', 0.5);
% set(ha, fset.tfont{:});
% ylabel('Lifetime (s)', fset.sfont{:})
% pos = get(get(ha, 'YLabel'), 'Position');
% pos(1) = dx;
% set(get(ha, 'YLabel'), 'Position', pos);







% M: model #
function v = costPDF(kVect, t, lftHist, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);

% interpolate result over full time vector
dt = t(2)-t(1);
t_full = 0:dt:t(end);

popPDF = zeros(M,numel(t_full));
Y = interp1(t_ode, Y, t_full);

% normalize subpopulation PDFs, weigh by output
for k = 1:M
    p = Y(:,2*k-1);
    popPDF(k,:) = p/sum(p)/dt * Y(end,2*k);
end
pdf = sum(popPDF, 1);

%normalize pdf to 1 over 't'
pdf = interp1(t_full, pdf, t);
n = sum(pdf)*dt;
pdf = pdf/n;

v = pdf - lftHist;



function v = costCDF(kVect, t, lftECDF, S0, M)
hf = str2func(['pop' num2str(M) 'Model']);
[t_ode, Y] = ode45(@(t,y) hf(t, y, kVect), [0 t(end)], S0);

% interpolate ODE output to input grid
CDF = interp1(t_ode, sum(Y(:,2:2:end),2), t);

% normalize to [0..1]
T = CDF(1);
CDF = (CDF-T)/(1-T);

v = CDF - lftECDF;




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


function y = getYAxisBound(vmax)
d = floor(log10(vmax));
% y-axis unit
yunit = round(vmax ./ 10.^d) .* 10.^(d-1);
y = ceil(vmax ./ yunit) .* yunit;