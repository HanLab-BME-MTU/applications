% Francois Aguet, October 2010

function plotCorrelationMatrix(M, varargin)

np = size(M,1);
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('M');
ip.addParamValue('AxisLabels', 1:np);
ip.addParamValue('Handle', []);
ip.addParamValue('TickLabels', arrayfun(@(i) num2str(i), 1:np, 'UniformOutput', false));
ip.parse(M, varargin{:});
ha = ip.Results.Handle;

sfont = {'FontName', 'Helvetica', 'FontSize', 14};
lfont = {'FontName', 'Helvetica', 'FontSize', 18};


triuMask = triu(ones(size(M)));
[x,y] = ind2sub(size(M), find(abs(M)>0.8 & triuMask==0));

% remove upper triangular part
M = M-triu(M);

% convert values [-1..1] to RGB
M = cat(3, -M.*(M<0), M.*(M>0), zeros(size(M)));

M(repmat(triuMask, [1 1 3])==1)=1;


if isempty(ha)
    figure('Position', [440 378 300 250]);
    ha = gca;
end

imagesc(M, 'Parent', ha); colormap(getMap()); axis image; caxis([-1 1]);
%hc = colorbar('YTick', -1:0.2:1, 'YTickLabel', arrayfun(@(i) num2str(i, '%.2f'), -1:0.2:1, 'UniformOutput', false));
colorbar('YTick', -1:0.2:1);
% hc = colorbar('SouthOutside', 'XTick', -1:0.2:1);
% cpos = get(hc, 'Position');
% cpos(2) = cpos(2)*0.8;
% cpos(1) = cpos(1)*1.1;
% set(hc, 'Position', cpos);
% get(hc)
% set(hc, 'YAxisLocation', 'left');


np = size(M,1);
ta = 1:np;
set(gca, 'XTick', ta, 'YTick', ta, 'XAxisLocation', 'bottom',...
    'TickLength', [0 0], 'XTickLabel', [], 'YTickLabel', [],...
    sfont{:}, 'LineWidth', 1.5);%, 'XColor', [1 1 1], 'YColor', [1 1 1]);

title('Correlation matrix', sfont{:})

XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');

% plot x-axis tick labels
arrayfun(@(k) text(ta(k), YLim(2)+0.05*diff(YLim), ip.Results.TickLabels{k},...
    sfont{:},...
    'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
    'Interpreter', 'TeX'),...
    ta, 'UniformOutput', false);


% plot y-axis tick labels
arrayfun(@(k) text(XLim(1)-0.05*diff(XLim), ta(k), ip.Results.TickLabels{k},...
    sfont{:},...
    'Units', 'data', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'Interpreter', 'TeX'),...
    ta, 'UniformOutput', false);


hold on;
plot(y, x, 'w.', 'LineWidth', 3, 'MarkerSize', 20);



function map = getMap()
values = -1:1/100:1;
N = length(values);
map = zeros(N,3);

ridx = values<0;
map(ridx,1) = -values(ridx);
gidx = values>0;
map(gidx,2) = values(gidx);
