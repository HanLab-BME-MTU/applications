function [ha, hi, hf] = setupFigure(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('nh', 1, @isposint);
ip.addOptional('nw', 1, @isposint);
ip.addOptional('na', [], @isposint);
ip.addParamValue('SameAxes', false, @islogical);
ip.addParamValue('AspectRatio', []);
ip.addParamValue('AxesWidth', 6);
ip.addParamValue('AxesHeight', 3.5);
ip.addParamValue('XSpace', [1.5 0.75 0.5]);
ip.addParamValue('YSpace', [1.5 0.75 0.5]);
ip.addParamValue('DisplayMode', 'print', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.addParamValue('InsetPosition', []);
ip.parse(varargin{:});
nh = ip.Results.nh;
nw = ip.Results.nw;
na = ip.Results.na;
if isempty(na)
    na = nh*nw;
end

w0 = ip.Results.AxesWidth + sum(ip.Results.XSpace);

ah0 = ip.Results.AxesHeight;
h0 = ah0 + sum(ip.Results.YSpace);
if ~isempty(ip.Results.AspectRatio)
    ah0 = ip.Results.AspectRatio*ip.Results.AxesWidth;
end

% default proportions: left/bottom: 1.5, width: 6, height: 3.5, top/right: 0.5
aw = ip.Results.AxesWidth/w0;
xl = ip.Results.XSpace(1)/w0; % left spacing (relative to single axes)
if ip.Results.SameAxes
    xc = ip.Results.XSpace(2)/w0;  % spacing btw axes
else
    xc = xl;
end
xr = ip.Results.XSpace(3)/w0;

ah = ah0/h0;
yb = ip.Results.YSpace(1)/h0;
if ip.Results.SameAxes
    yc = ip.Results.YSpace(2)/h0;
else
    yc = yb;
end
yt = ip.Results.YSpace(3)/h0;

% width (relative to normalized single axes)
w = xl + nw*aw + (nw-1)*xc + xr;
% height
h = yb + nh*ah + (nh-1)*yc + yt;


% convert to normalized units
aw = aw/w;
xl = xl/w;
xc = xc/w;

ah = ah/h;
yb = yb/h;
yc = yc/h;

% resize figure window
fset = loadFigureSettings(ip.Results.DisplayMode);
fpos = fset.fPos;

fpos(3) = w*w0;
fpos(4) = h*h0;
if strcmpi(ip.Results.DisplayMode, 'print')
    units = 'centimeters';
else
    units = 'pixels';
end

fposPx = fpos/2.54*get(0,'ScreenPixelsPerInch');
fpos0 = get(0, 'ScreenSize');
c = min(0.8*fpos0(3:4)./fposPx(3:4));
if c<1
    fpos(3:4) = c*fpos(3:4);
end

hf = figure('PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off',...
    'Units', units, 'Position', fpos, 'Units', 'pixels');

ha = zeros(na,1);
x0 = zeros(na,1);
y0 = zeros(na,1);
hi = zeros(na,1);
ipos = ip.Results.InsetPosition;
if numel(ipos)==2
    ipos = [ipos 0.95-ipos(1) 0.95-ipos(2)];
end
for i = 1:na
    y0(i) = nh-ceil(i/nw);
    x0(i) = mod(i-1,nw);
    pos = [xl+x0(i)*(aw+xc) yb+y0(i)*(ah+yc) aw ah];
    ha(i) = axes('Position', pos); %#ok<LAXES>
    hold(ha(i), 'on');
    
    if ~isempty(ipos);
        hi(i) = axes('Position', [pos(1)+pos(3)*ipos(3) pos(2)+pos(4)*ipos(4)...
            ipos(1)*pos(3) ipos(2)*pos(4)]); %#ok<LAXES>
        hold(hi(i), 'on');
    end
end

if ip.Results.SameAxes
    set(ha(y0>0), 'XTickLabel', []);
    set(ha(x0>0), 'YTickLabel', []);
end

set(ha, 'TickDir', 'out', 'TickLength', fset.TickLength*6/max(ip.Results.AxesWidth, ah0),...
    'LineWidth', 1, 'Layer', 'top');

if ~isempty(ipos)
    set(hi, 'TickDir', 'out', 'TickLength', fset.TickLength*max(aw,ah)/max(ipos(1)*pos(3), ipos(2)*pos(4)),...
    'LineWidth', 1, 'Layer', 'top');
end
