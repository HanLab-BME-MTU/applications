% Francois Aguet (last modified 06/01/2012)

function plotLifetimeDistModel(lftData, fitRes, varargin)

% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('lftFit');
% ip.addParamValue('PlotAll', false, @islogical);
% ip.addParamValue('PlotCDF', false, @islogical);
% ip.addParamValue('fYLim', []);
% ip.addParamValue('rYLim', []);
% % ip.addParamValue('ShowInset', false, @islogical);
% ip.parse(lftFit, varargin{:});

np = fitRes.N;
a = fitRes.a;

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




figure('Name', [fitRes.ModelType ' distribution fit (to ' fitRes.FitMode ')'])%('Position', [240 378 850 360], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
%---------------------------------
% Lifetime histogram
%---------------------------------
%axes('Units', 'Pixels', 'Position', [95 65 450 270]);
set(gca, 'ColorOrder', colorOrder);
hold on;
%hp(1) = plot(lftData.t_hist, lftData.meanLftHist_A*(1-a), '.', 'MarkerSize', 20, 'Color', [0 0 0]);
hp(1) = plot(lftData.t, lftData.meanLftHist_A*(1-a), '.-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 18);
YLim = get(gca, 'YLim');
hi = plot(fitRes.t, fitRes.popPDF, 'LineWidth', 2);
hp(2) = hi(1);
hp(3) = plot(fitRes.t, fitRes.PDF, '--', 'Color', fset.ceB, 'LineWidth', 4);

%axis([0 100 0 fYLim]);
set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:}, 'XLim', [0 120], 'YLim', YLim);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});

% hl = legend(hp, 'Meas. lifetime', 'Pop. lifetimes', 'Model');
% set(hl, 'Box', 'off');

%---------------------------------
% Inset with amplitudes
%---------------------------------


% main axes: [85 65 450 270]
% ha = axes('Units', 'Pixels', 'Position', [95+450-32*np-10 260 32*np 65]);
ha = axes('Position', [0.7 0.7 0.2 0.2]);
xlabels = arrayfun(@(i) ['P' num2str(i)], 1:np, 'UniformOutput', false);

barplot2(fitRes.A, 'AdjustFigure', false, 'XLabels', xlabels,...
    'FaceColor', colorOrderFill, 'EdgeColor', colorOrder,...
    'BarWidth', 0.6, 'GroupDistance', 0.5, 'Angle', 45);
set(ha, fset.tfont{:}, 'YTick', 0:0.2:0.8, 'YLim', [0 0.8]);
ylabel('Contrib.', fset.sfont{:})


% if ~isempty(ip.Results.rYLim)
%     set(ha, 'YLim', [0 ip.Results.rYLim]);
% end

