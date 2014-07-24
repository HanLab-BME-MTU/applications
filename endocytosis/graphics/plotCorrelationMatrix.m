%plotCorrelationMatrix(M, varargin) displays a correlation matrix as a green-red colormapped image

% Francois Aguet, October 2010 (Last modified 01/19/2012)

function plotCorrelationMatrix(M, varargin)

np = size(M,1);
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('M');
ip.addParamValue('AxisLabels', 1:np);
ip.addParamValue('Handle', []);
ip.addParamValue('TickLabels', arrayfun(@(i) num2str(i), 1:np, 'unif', 0));
ip.addParamValue('ScalePrint', 10);
ip.addParamValue('Colorbar', true, @islogical);
ip.addParamValue('ShowTitle', false, @islogical);
ip.addParamValue('ShowDiagonal', false, @islogical);
ip.parse(M, varargin{:});
ha = ip.Results.Handle;
scale = ip.Results.ScalePrint;

fset = loadFigureSettings();

% find strong correlations
[x,y] = ind2sub(size(M), find(abs(M)>0.8 & diag(diag(M))~=1));

% convert values [-1..1] to RGB
M = cat(3, -M.*(M<0), M.*(M>0), zeros(size(M)));
M = uint8(255*M);


if isempty(ha)
    figure('Position', [440 378 300 200], 'Color', 'w', 'PaperPositionMode', 'auto');
    ha = gca;
end

imagesc(imresize(M, scale, 'nearest'), 'Parent', ha); colormap(getMap());
axis image; caxis([-1 1]);


np = size(M,1);
ta = (1:np)*scale - scale/2 + 0.5;
set(gca, 'XTick', ta, 'YTick', ta, 'XAxisLocation', 'bottom',...
    'TickLength', [0 0], 'XTickLabel', [], 'YTickLabel', [],...
    fset.sfont{:}, 'LineWidth', 1.5);%, 'XColor', [1 1 1], 'YColor', [1 1 1]);

if ip.Results.Colorbar
    hc = colorbar('YTick', -1:0.2:1);
    set(hc, fset.tfont{:});
end

if ip.Results.ShowTitle
    title('Correlation matrix', fset.sfont{:})
end

XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');

% plot x-axis tick labels
arrayfun(@(k) text(ta(k), YLim(2)+0.05*diff(YLim), ip.Results.TickLabels{k},...
    fset.sfont{:},...
    'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
    'Interpreter', 'TeX'),...
    1:np, 'unif', 0);


% plot y-axis tick labels
arrayfun(@(k) text(XLim(1)-0.05*diff(XLim), ta(k), ip.Results.TickLabels{k},...
    fset.sfont{:},...
    'Units', 'data', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'Interpreter', 'TeX'),...
    1:np, 'unif', 0);

hold on;
plot(ta(y), ta(x), 'ws', 'LineWidth', 1, 'MarkerSize', 8);


function map = getMap()
values = -1:1/100:1;
N = length(values);
map = zeros(N,3);

ridx = values<0;
map(ridx,1) = -values(ridx);
gidx = values>0;
map(gidx,2) = values(gidx);
