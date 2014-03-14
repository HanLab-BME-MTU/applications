

% Francois Aguet, 03/08/2014

function [] = scattercontour(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('y');
ip.addOptional('xv', []);
ip.addOptional('yv', []);
ip.addParamValue('Parent', [], @ishandle);
ip.addParamValue('npoints', 100);
% ip.addParamValue('Color', []);
ip.parse(x, y, varargin{:});
xv = ip.Results.xv;
yv = ip.Results.yv;
ha = ip.Results.Parent;

p.N = ip.Results.npoints;
if ~isempty(xv) && ~isempty(yv)
    p.xy = [xv(:) yv(:)];
end
p = gkde3([x(:) y(:)]);%, p);

xv = p.x(1,:);
yv = p.y(:,1);
dx = xv(2)-xv(1);
dy = yv(2)-yv(1);

% re-normalize
p.pdf = p.pdf/(sum(p.pdf(:))*dx*dy);

% value of pdf at sample locations -> percentiles
val = interp2(p.x, p.y, p.pdf, x, y);
c = prctile(val, [5 25 50 75 95]);

% figure; surf(xv, yv, p.pdf); hold on; plot3(x, y, val, 'r.')
if isempty(ha)
    figure;
    hold on;
    plot(x,y,'k.');
    ha = gca;
    axis equal;
else
%     plot(ha, x, y, 'k.', 'MarkerSize', 4);
end
nc = numel(c);
cmap = colormap(jet(nc));
for i = 1:nc
    [~,hc] = contour(ha, p.x, p.y, p.pdf, c(i), 'LineColor', cmap(i,:),...
        'LineWidth', 1, 'HitTest', 'off');
end
