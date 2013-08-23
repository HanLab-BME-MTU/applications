% Francois Aguet, August 2010

function R = scatterplot2(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('y');
ip.addOptional('xname', 'x', @ischar);
ip.addOptional('yname', 'y', @ischar);
ip.addParamValue('Handle', []);
ip.addParamValue('EqualAxes', false, @islogical);
ip.addParamValue('Color', []);
ip.parse(x, y, varargin{:});
ha = ip.Results.Handle;

mux = mean(x);
muy = mean(y);
N = length(x);

sxx = sum(x.^2) - N*mux^2; % N*var(x,1)
syy = sum(y.^2) - N*muy^2;
sxy = sum(x.*y) - N*mux*muy; % N*cov(x,y)

R = sxy^2/(sxx*syy);

if isempty(ha)
    figure('Color', 'w');
    hold on;
    ha = gca;
end

if ~isempty(ip.Results.Color) && size(ip.Results.Color,1)>1
    colormap(ip.Results.Color); % needed for vectorial EPS printing
    scatter(x, y, 30, 1:size(ip.Results.Color,1), 'filled', 'Parent', gca);
elseif ~isempty(ip.Results.Color)
    plot(ha, x, y, '.', 'Color', ip.Results.Color, 'MarkerSize', 13);
else
    plot(ha, x, y, 'k.', 'MarkerSize', 13);
end
axis square;

XLim = get(ha, 'XLim');
YLim = get(ha, 'YLim');
XLim(1) = 0;
YLim(1) = 0;
% b1 = sxy/sxx;
% b0 = mean(y)-mean(x)*b1;
% x0 = [0 max(x)];
% plot(x0, b1*x0+b0, 'r--');

set(ha, 'FontName', 'Helvetica', 'FontSize', 12, 'LineWidth', 1);
xlabel(ha, ip.Results.xname, 'FontName', 'Helvetica', 'FontSize', 14);
ylabel(ha, ip.Results.yname, 'FontName', 'Helvetica', 'FontSize', 14);
% title(['R^2 = ' num2str(R, '%.4f')]);

if ip.Results.EqualAxes
    lim = [min(XLim(1),YLim(1)) min(XLim(2),YLim(2))];
    set(ha, 'XLim', lim, 'YLim', lim);
else
    set(ha, 'XLim', XLim, 'YLim', YLim);
end
set(ha, 'Layer', 'top', 'TickDir', 'out');
