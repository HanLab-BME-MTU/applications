function [ha, hf] = setupFigure(nh, nw, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('nh', @isposint);
ip.addRequired('nw', @isposint);
ip.addOptional('na', nh*nw, @isposint);
ip.addParamValue('SameAxes', false, @islogical);
ip.addParamValue('AspectRatio', []);
ip.addParamValue('AxesWidth', 6);
ip.addParamValue('AxesHeight', 3.5);
ip.addParamValue('DisplayMode', 'print', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.parse(nh, nw, varargin{:});
na = ip.Results.na;

ar = ip.Results.AspectRatio;
if isempty(ar)
    ar = ip.Results.AxesWidth/ip.Results.AxesHeight;
end
w0 = 8;
h0 = 1.5+0.5+ar*6;


% default proportions: left/bottom: 1.5, width: 6, height: 3.5, top/right: 0.5
aw = 6/8;
xl = 1.5/8; % left spacing (relative to single axes)
if ip.Results.SameAxes
    xc = xl/2; % spacing btw axes
else
    xc = xl;
end
xr = 0.5/8;

ah = ar*6/h0;
yb = 1.5/h0;
if ip.Results.SameAxes
    yc = yb/2;
else
    yc = yb;
end
yt = 0.5/h0;

% width
w = xl + nw*aw + (nw-1)*xc + xr;
% height
h = yb + nh*ah + (nh-1)*yc + yt;


% normalize
aw = aw/w;
xl = xl/w;
xc = xc/w;

ah = ah/h;
yb = yb/h;
yc = yc/h;

% resize figure window
fset = loadFigureSettings(ip.Results.DisplayMode);
fpos = fset.fPos;
% if w>h
%     fpos(4) = fpos(3)*h/w;
% else
%     fpos(3) = fpos(4)*w/h;
% end

% f = 5.5/8;
% if w/fpos(3) > h/fpos(4) % figure width limiting
%     fpos(4) = fpos(3)*h/w*f;
% else % height limiting
%     fpos(3) = fpos(4)*w/h/f;
% end

fpos(3) = w*w0;
fpos(4) = h*h0;
if strcmpi(ip.Results.DisplayMode, 'print')
    units = 'centimeters';
else
    units = 'pixels';
end

hf = figure('PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off',...
    'Units', units, 'Position', fpos);

ha = zeros(na,1);
x0 = zeros(na,1);
y0 = zeros(na,1);
for i = 1:na
    y0(i) = ceil(i/nw)-1;
    x0(i) = mod(i-1,nw);
    ha(i) = axes('Position', [xl+x0(i)*(aw+xc) yb+y0(i)*(ah+yc) aw ah]); %#ok<LAXES>
    hold(ha(i), 'on');
end

if ip.Results.SameAxes
    set(ha(y0>0), 'XTickLabel', []);
    set(ha(x0>0), 'YTickLabel', []);
end
set(ha, 'TickDir', 'out', 'TickLength', fset.TickLength, 'LineWidth', 1);
