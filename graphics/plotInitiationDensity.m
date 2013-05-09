function ha = plotInitiationDensity(lftRes, xlabels, cv, cf, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('SigLink', cell(1,2), @iscell)
ip.addParamValue('YLim', {[0 0.65], [0 0.325]});
ip.addParamValue('YTick', {0:0.1:0.6, 0:0.05:0.3});
ip.parse(varargin{:});
YLim = ip.Results.YLim;
YTick = ip.Results.YTick;

fset = loadFigureSettings('print');

M_Ia = cellfun(@(i) i.initDensity_Ia(:,1), lftRes, 'UniformOutput', false);
M_above = cellfun(@(i) i.initDensity_above(:,1), lftRes, 'UniformOutput', false);


figure(fset.fOpts{:}, 'Position', [10 10 7 7]);

ha(1) = axes(fset.axOpts{:}, 'Position', [1.5 3 2.25 3.5], 'TickLength', fset.TickLength*6/3.5);
boxplot2({M_Ia}, ip.Results.SigLink{1}, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BarWidth', 0.6, 'LineWidth', 1, 'GroupDistance', 0, 'BorderWidth', 0.5,...
    'FaceColor', {cf, ones(size(cf))}, 'EdgeColor', {cv,cv}, 'ErrorbarColor', {cv,cv},...
    'YLim', YLim{1}, 'Angle', 0);
ylabel(['Initiations (' char(181) 'm^{-2} min^{-1})'], fset.lfont{:});
rotateXTickLabels(gca, 'Angle', 45, 'AdjustFigure', true, 'YOffset', 0.03);
set(gca, 'YTick', YTick{1});

ha(2) = axes(fset.axOpts{:}, 'Position', [4.6 3 2.25 3.5], 'TickLength', fset.TickLength*6/3.5);
boxplot2({M_above}, ip.Results.SigLink{2}, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BarWidth', 0.6, 'LineWidth', 1, 'GroupDistance', 0, 'BorderWidth', 0.5,...
    'FaceColor', {cf, ones(size(cf))}, 'EdgeColor', {cv,cv}, 'ErrorbarColor', {cv,cv},...
    'YLim', YLim{2}, 'Angle', 0);
rotateXTickLabels(gca, 'Angle', 45, 'AdjustFigure', true, 'YOffset', 0.03);
set(gca, 'YTick', YTick{2});

formatTickLabels(ha(1));
formatTickLabels(ha(2));
