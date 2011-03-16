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
    ch = find(strcmp(data.channels, data.source));
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
    frame = imread(data(k).framePaths{ch}{frameIdx});
    
    [sy sx] = size(frame);
    psize = data(k).pixelSize/data(k).M*1e6; % in µm
    xa = (0:sx-1)*psize;
    ya = (0:sy-1)*psize;
    if sy>sx
        frame = imrotate(frame, 90);
        sx = size(frame,2);
    end
    
    dx = sx*psize/40;
    
    imagesc(xa,ya,frame);
    colormap(gray(256));
    axis image off;
    plotScaleBar(scale, h, 'Label', [num2str(scale) ' ' units], 'FontName', fontName, 'FontSize', fontSize);
    
    text(xa(end)-dx, dx, getCellDirectory(data(k)), 'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Right',...
                'FontName', fontName, 'FontSize', fontSize, 'interpreter', 'none');

end