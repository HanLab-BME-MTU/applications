function plotCellArea(lftRes, labels)

fset = loadFigureSettings('print');

mu = cellfun(@(i) mean(i.cellArea), lftRes);
%sigma = cellfun(@(i) std(i.cellArea), lftRes);
sigma = cellfun(@(i) std(i.cellArea)/sqrt(numel(i.cellArea)), lftRes);

figure(fset.fOpts{:}, 'Position', [5 5 5.5 7]);
axes(fset.axOpts{:}, 'Position', [2 3 3 3.5], 'TickLength', fset.TickLength/3.5*6);
barplot2(mu', sigma', 'Angle', 0, 'BarWidth', 1, 'GroupDistance', 1,...
    'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.4*[1 1 1], 'AxisFontSize', 8,...
    'LineWidth', 1);
ylabel(['Cell area (' char(181) 'm^2)'], fset.lfont{:})

set(gca, 'XTickLabel', labels);
rotateXTickLabels(gca, 'AdjustFigure', false);
