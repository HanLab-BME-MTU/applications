% Francois Aguet, August 2010

function R = scatterplot2(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('y');
ip.addOptional('xname', 'x', @ischar);
ip.addOptional('yname', 'y', @ischar);
ip.addParamValue('EqualAxes', false, @islogical);
ip.parse(x, y, varargin{:});

mux = mean(x);
muy = mean(y);
N = length(x);

sxx = sum(x.^2) - N*mux^2; % N*var(x,1)
syy = sum(y.^2) - N*muy^2;
sxy = sum(x.*y) - N*mux*muy; % N*cov(x,y)

R = sxy^2/(sxx*syy);

figure;
hold on;
plot(x, y, 'k.');
axis square;

XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
plot([0 max(XLim(2),YLim(2))], [0 max(XLim(2),YLim(2))], 'r--');

set(gca, 'FontName', 'Helvetica', 'FontSize', 12, 'LineWidth', 1);
xlabel(ip.Results.xname, 'FontName', 'Helvetica', 'FontSize', 14);
ylabel(ip.Results.yname, 'FontName', 'Helvetica', 'FontSize', 14);
title(['R^2 = ' num2str(R, '%.4f')]);

if ip.Results.EqualAxes
    lim = [min(XLim(1),YLim(1)) min(XLim(2),YLim(2))];
    set(gca, 'XLim', lim, 'YLim', lim);
end
set(gca, 'Layer', 'top', 'TickDir', 'out');