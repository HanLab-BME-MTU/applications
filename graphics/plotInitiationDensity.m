%plotInitiationDensity(lftRes, xlabels, cv, cf, varargin) compares CS vs. CCP initiation density under different conditions
%
% Inputs:
%
%   lftRes : output structure from runLifetimeAnalysis()
%  xlabels : legend; cell array of strings
%       cv : Nx3 color matrix for edges
%       cv : Nx3 color matrix

% Francois Aguet

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

M_Ia = cellfun(@(i) i.initDensityIa(:,1), lftRes, 'unif', 0);
M_above = cellfun(@(i) i.initDensityCCP(:,1), lftRes, 'unif', 0);

figure(fset.fOpts{:}, 'Position', [10 10 7 7]);

ha(1) = axes(fset.axOpts{:}, 'Position', [1.5 3 2.25 3.5], 'TickLength', fset.TickLength*6/3.5);
boxplot2(M_Ia, ip.Results.SigLink{1}, 'AdjustFigure', false, 'XTickLabel', xlabels,...
   'BarWidth', 0.6, 'LineWidth', 1, 'GroupDistance', 0, 'BorderWidth', 0.5,...
    'YLim', YLim{1}, 'Angle', 0, 'DetectOutliers', false,...
    'FaceColor', cf, 'EdgeColor', cv, 'ErrorbarColor', cv);
ylabel(['Initiations (' char(181) 'm^{-2} min^{-1})'], fset.lfont{:});
rotateXTickLabels(gca, 'Angle', 45, 'AdjustFigure', true);
set(gca, 'YTick', YTick{1});


ha(2) = axes(fset.axOpts{:}, 'Position', [4.6 3 2.25 3.5], 'TickLength', fset.TickLength*6/3.5);
boxplot2(M_above, ip.Results.SigLink{2}, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BarWidth', 0.6, 'LineWidth', 1, 'GroupDistance', 0, 'BorderWidth', 0.5,...
    'FaceColor', cf, 'EdgeColor', cv, 'ErrorbarColor', cv,...
    'YLim', YLim{2}, 'Angle', 0, 'DetectOutliers', false);
rotateXTickLabels(gca, 'Angle', 45, 'AdjustFigure', true);
set(gca, 'YTick', YTick{2});

formatTickLabels(ha(1));
formatTickLabels(ha(2));
