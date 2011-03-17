function cellViewer(data, varargin)


if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

%======================================
% Parse inputs, set defaults
%======================================
idx = find(strcmpi(varargin, 'Frame'));
if ~isempty(idx)
    frameIdx = varargin{idx+1};
else
    frameIdx = 1;
end

idx = find(strcmpi(varargin, 'Channel'));
if ~isempty(idx)
    ch = varargin{idx+1};
else
    ch = find(strcmp(data(1).channels, data(1).source));
end

idx = find(strcmpi(varargin, 'Scale'));
if ~isempty(idx)
    scale = varargin{idx+1};
else
    scale = 5;
end

idx = find(strcmpi(varargin, 'Units'));
if ~isempty(idx)
    units = varargin{idx+1};
else
    units = 'µm';
end

idx = find(strcmpi(varargin, 'FontName'));
if ~isempty(idx)
    fontName = varargin{idx+1};
else
    fontName = 'Helvetica';
end


idx = find(strcmpi(varargin, 'FontSize'));
if ~isempty(idx)
    fontSize = varargin{idx+1};
else
    fontSize = 10;
end

idx = find(strcmpi(varargin, 'Mode'));
if ~isempty(idx)
    mode = varargin{idx+1};
else
    mode = 'raw';
end


if nargin<2 || isempty(frameIdx)
    frameIdx = 1;
end

N = numel(data);


ny = ceil(sqrt(2*N/(1+sqrt(5))));
nx = ceil(N/ny);

if nx > 1
    dx = 0.05/(nx-1);
    wx = 0.95/nx;    
else
    dx = 0;
    wx = 1;    
end
if ny > 1
    dy = 0.05/(ny-1);
    wy = 0.95/ny;
else
    dy = 0;
    wy = 1;
end

figure;
for k = 1:N
    y = ceil(k/nx);
    x = k-(y-1)*nx;
    h = axes('Position', [(dx+wx)*(x-1) 1-wy*y-(y-1)*dy wx wy]);
    
    plotFrame(data(k), [], frameIdx, ch, 'Handle', h, 'Units', 1e-6, 'Mode', mode);

    axis image off;
    plotScaleBar(scale, 'Handle', h, 'Label', [num2str(scale) ' ' units], 'FontName', fontName, 'FontSize', fontSize);
    r = data(k).imagesize(2)/data(k).imagesize(1);
    d = 0.02;
    text(1-d, 1-d*r, getCellDirectory(data(k)), 'Color', 'w',...
                'Units', 'normalized',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Right',...
                'FontName', fontName, 'FontSize', fontSize, 'interpreter', 'none');

end