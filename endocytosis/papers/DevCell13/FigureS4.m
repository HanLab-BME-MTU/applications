%=========================================================================================
% Figure S4
%=========================================================================================
% This script generates the panels for Figure S4 from Aguet et al., Dev. Cell, 2013.

%=========================================================================================
% Panel A: Intensity cohorts
%=========================================================================================
% load data using Figure5.m first

[cohorts, res] = plotIntensityCohorts(dataDC3, [1 2], 'ShowBackground', false,...
    'DisplayMode', 'print', 'ScaleSlaveChannel', false, 'ScalingFactor', [5 1],...
    'ShowLegend', false, 'ShowPct', false, 'YTick', 0:20:80,...
    'ExcludeVisitors', true, 'SlaveName', {'enDyn2-EGFP'}, 'AvgFun', @nanmean);

fset = loadFigureSettings('print');
figure(fset.fOpts{:}, 'Position', [2 2 17 5.5]);
axes(fset.axOpts{:})
hold on;
nc = 6;
hues = 0.33;
v = mod(hues+linspace(-0.1, 0.1, nc)', 1);
cmap = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
cv = hsv2rgb([v 0.4*ones(nc,1) ones(nc,1)]);
xc = 2;
b = 5;
framerate = 2;
for c = nc:-1:xc
    fill([cohorts(1).t{c} cohorts(1).t{c}(end:-1:1)], [cohorts(1).Amin{2,c} cohorts(1).Aplus{2,c}(end:-1:1)],...
        cv(c,:), 'EdgeColor', cmap(c,:), 'HandleVisibility', 'off');
    plot(cohorts(1).t{c}, cohorts(1).A{2,c}, '-', 'Color', cmap(c,:), 'LineWidth', 1);
end
set(gca, 'XLim', [-b*framerate-5 120], 'XTick', 0:20:200);
set(gca, 'YTick', 0:10:50, 'YLim', [0 50]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
text(XLim(1)+0.025*diff(XLim), YLim(2), 'Appearance-aligned', fset.sfont{:}, 'VerticalAlignment', 'bottom');


axes(fset.axOpts{:}, 'Position', [8.25 1.5 6 3.5]);
hold on;
for c = nc:-1:xc
    fill([cohorts(1).t{c} cohorts(1).t{c}(end:-1:1)]-cohorts(1).t{c}(end)+10, [cohorts(1).Amin{2,c} cohorts(1).Aplus{2,c}(end:-1:1)], cv(c,:), 'EdgeColor', cmap(c,:), 'HandleVisibility', 'off');
    hp(c-xc+1) = plot(cohorts(1).t{c}-cohorts(1).t{c}(end)+10, cohorts(1).A{2,c}, '-', 'Color', cmap(c,:), 'LineWidth', 1);
end
set(gca, 'XLim', [-b*framerate-5-110 10], 'XTick', -200:20:0);
set(gca, 'YTick', 0:10:50, 'YLim', [0 50]);
xlabel('Time (s)', fset.lfont{:});
set(gca, 'YTickLabel', []);
xlabel('Time (s)', fset.lfont{:});
XLim = get(gca, 'XLim');
text(XLim(1)+0.025*diff(XLim), YLim(2), 'Disappearance-aligned', fset.sfont{:}, 'VerticalAlignment', 'bottom');

cohortLabels = arrayfun(@(i) [' ' num2str(cohorts(1).bounds(i)) '-' num2str(cohorts(1).bounds(i+1)-framerate) ' s'], xc:nc, 'unif', 0);
hl = legend(hp, cohortLabels, 'Location', 'SouthEast');
set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [6.75+7.65 1.5 1.25 1.5]);
% print('-depsc2', '-loose', 'dynCohortsOX2aligned.eps');

%%
%=========================================================================================
% Panel B: Lifetimes for all CS/CCP classes
%=========================================================================================
lftResDC3 = runLifetimeAnalysis(dataDC3, 'MaxIntensityThreshold', 15.59,...
    'RemoveOutliers', false, 'LifetimeData', 'lifetimeData.mat', 'Cutoff_f',5,...
    'DisplayMode', 'print', 'Display', 'off');

plotLifetimes(lftResDC3, 'ShowCargoDependent', true, 'SlaveNames', 'enDyn2-EGFP',...
    'DisplayMode', 'print', 'YTick', 0:0.005:0.02, 'PlotAll', true);


