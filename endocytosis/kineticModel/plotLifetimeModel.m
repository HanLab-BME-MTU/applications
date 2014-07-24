% Francois Aguet (last modified 01/23/2012)

function plotLifetimeModel(lftFit, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftFit');
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('PlotCDF', false, @islogical);
ip.addParamValue('fYLim', []);
ip.addParamValue('rYLim', []);
% ip.addParamValue('ShowInset', false, @islogical);
ip.parse(lftFit, varargin{:});

np = lftFit.np;
a = lftFit.a;

fYLim = ip.Results.fYLim;
if isempty(fYLim)
    fYLim = getYAxisBound(max(lftFit.lftHist)*(1-a));
end

fset = loadFigureSettings();

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
axes('Units', 'Pixels', 'Position', [95 65 450 270]);
set(gca, 'ColorOrder', colorOrder);
hold on;
hp(1) = plot(lftFit.t, lftFit.lftHist*(1-a), '.', 'MarkerSize', 20, 'Color', [0 0 0]);
hi = plot(lftFit.t_fine, lftFit.popPDF, 'LineWidth', 2);
hp(2) = hi(1);
hp(3) = plot(lftFit.t_fine, lftFit.pdf, '--', 'Color', fset.ceB, 'LineWidth', 4);

axis([0 100 0 fYLim]);
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
ha = axes('Units', 'Pixels', 'Position', [95+450-32*np-10 260 32*np 65]);
% ha = axes('Units', 'Pixels', 'Position', [85+450+60 270 110 65]);
xlabels = arrayfun(@(i) ['P' num2str(i)], 1:np, 'UniformOutput', false);

barplot2(lftFit.pA, 'AdjustFigure', false, 'XLabels', xlabels,...
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

ha = axes('Units', 'Pixels', 'Position', [95+450+dx 65 200 270]);
rateLabels = plotKineticModelRates(lftFit.k{lftFit.np}, lftFit.k_std, 'Handle', ha);
if ~isempty(ip.Results.rYLim)
    set(ha, 'YLim', [0 ip.Results.rYLim]);
end

if ip.Results.PlotAll
 
    % 85 450 dx 200 dx 240
    nk = numel(lftFit.k{np})-1;
    ha = axes('Units', 'Pixels', 'Position', [85+450+200+dx+60 65+270-30*nk 30*nk 30*nk]);
    plotCorrelationMatrix(lftFit.corr{np}, 'Handle', ha, 'TickLabels', rateLabels, 'ColorBar', false);
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
    if sum(lftFit.BIC~=0)>1
        figure('Position', [440 378 400 300], 'PaperPositionMode', 'auto');
        %np = numel(lftFit.BIC);
        px = find(lftFit.BIC~=0);
        %axes('Units', 'pixels', 'Position', [85+450+200+dx+60 65 60*np 120]);
        axes('Units', 'pixels', 'Position', [100 60 60*numel(px) 220])
        
        hold on;
        plot(px, lftFit.BIC(px), 'r.', 'MarkerSize', 40);
        set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:}, 'XLim', [px(1)-0.5 px(end)+0.5], 'XTick', 1:np);
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
    plot(lftFit.t, lftFit.lftECDF*(1-a)+a, 'k.', 'MarkerSize', 20)
    plot(lftFit.t_fine, lftFit.popCDF, 'LineWidth', 2);
    plot(lftFit.t_fine, lftFit.cdf, '--', 'Color', fset.ceB, 'LineWidth', 4);
    
    axis([0 100 0 1]);
    set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
end

% ha = axes('Units', 'Pixels', 'Position', [85+450+60 65 110 180]);
% pct = vertcat(lftFit.pPercentiles{:})';
% M = [pct(3,:); pct(2,:); pct(4,:); pct(1,:); pct(5,:)];
% boxplot2(M, 'AdjustFigure', false, 'XLabels', xlabels,...
%     'FaceColor', colorOrderFill, 'EdgeColor', colorOrder,...
%     'BarWidth', 0.5, 'GroupDistance', 0.5);
% set(ha, fset.tfont{:});
% ylabel('Lifetime (s)', fset.sfont{:})
% pos = get(get(ha, 'YLabel'), 'Position');
% pos(1) = dx;
% set(get(ha, 'YLabel'), 'Position', pos);

function y = getYAxisBound(vmax)
d = floor(log10(vmax));
% y-axis unit
yunit = round(vmax ./ 10.^d) .* 10.^(d-1);
y = ceil(vmax ./ yunit) .* yunit;
