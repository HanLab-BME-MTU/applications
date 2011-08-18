% Overlay of multiple mean lifetime distributions

function plotLifetimeDistCmp(res)

cut = 2;

tfont = {'FontName', 'Helvetica', 'FontSize', 14};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

N = length(res);
colors = jet(N);

figure;
hold on;

for k = 1:length(res)
    plot(res(k).t(cut:end), res(k).meanHist(cut:end), '-', 'Color', colors(k,:), 'LineWidth', 2);
end

% legend
h = legend(arrayfun(@(k) ['Movie ' num2str(k)], 1:N, 'UniformOutput', false), tfont{:});
legend(h, 'boxoff');


t_max = max(arrayfun(@(x) x.t(find(x.meanHist~=0, 1, 'last')+1), res));
f_max = max(arrayfun(@(x) max(x.meanHist(2:end)), res));
axis([0 t_max 0 1.05*f_max]);

set(gca, 'Layer', 'top', 'LineWidth', 1.5, sfont{:});
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});