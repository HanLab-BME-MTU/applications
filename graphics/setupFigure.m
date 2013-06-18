function [ha, hf] = setupFigure(nh, nw, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('nh', @isposint);
ip.addRequired('nw', @isposint);
ip.addOptional('na', nh*nw, @isposint);
ip.addParamValue('SameAxes', false, @islogical);
% ip.addParamValue('Units', 'pixels', @ischar);
ip.addParamValue('DisplayMode', 'print', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.parse(nh, nw, varargin{:});
na = ip.Results.na;

% default proportions:
aw = 0.75;
xl = 3/16; % left spacing (relative to single axes)
if ip.Results.SameAxes
    xc = 3/32;
else
    xc = xl;
end
xr = 1/16;

ah = 7/11;
yb = 3/11;
if ip.Results.SameAxes
    yc = 3/23;
else
    yc = yb;
end
yt = 1/11;

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
% fpos = get(0, 'DefaultFigurePosition');
% f = 5.5/8;
% if w/fpos(3) > h/fpos(4) % figure width limiting
%     fpos(4) = fpos(3)*h/w*f;
% else % height limiting
%     fpos(3) = fpos(4)*w/h/f;
% end

fset = loadFigureSettings(ip.Results.DisplayMode);
fpos = fset.fPos;
fpos(3) = w*fpos(3);
fpos(4) = h*fpos(4);
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
