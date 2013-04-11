% Francois Aguet, 02/19/2013

function [ha2, dRange] = densityplot(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('xv', 1:numel(x));
ip.addOptional('yv', 1:numel(y));
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('AdjustBorder', false, @islogical);
ip.addParamValue('Div', 1, @isscalar);
ip.addParamValue('DisplayFunction', @(x) x, @(x) isa(x, 'function_handle'));
ip.parse(varargin{:});
xv = ip.Results.xv(:)';
yv = ip.Results.yv(:)';
ha = ip.Results.Handle;

set(gca, 'Color', 'r');
XLim = xv([1 end]);
YLim = yv([1 end]);
dx = xv(2)-xv(1);
dy = yv(2)-yv(1);
if ip.Results.AdjustBorder % bug in eps files? png works fine
    XLim = [xv(1)-dx/2 xv(end)+dx/2];
    YLim = [yv(1)-dx/2 yv(end)+dy/2];
end
xv = [xv(1)-dx xv xv(end)+dx];
yv = [yv(1)-dy yv yv(end)+dy];
M = hist3([y(:) x(:)], {yv, xv});
% M = M(2:end-1,2:end-1);
M = ip.Results.DisplayFunction(M/ip.Results.Div);
imagesc(xv, yv, M, 'Parent', ha);
dRange = [min(M(:)) max(M(:))];

% re-plot axes on top to create box w/o ticks
ha2 = copyobj(ha, get(ha, 'Parent'));
set(ha2, 'XTick', [], 'YTick', [], 'Color', 'b', 'Box', 'on');
axis([ha ha2], [XLim YLim]);
% set(gcf, 'CurrentAxes', ha);% doesn't work
