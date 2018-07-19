%

% Francois Aguet, 2011

function cellViewer(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('Frame', [], @isposint);
ip.addParameter('Channel', 1, @isscalar);
ip.addParameter('Scale', 5, @isscalar);
ip.addParameter('Units', '�m', @ischar);
ip.addParameter('FontName', 'Helvetica', @ischar);
ip.addParameter('FontSize', 10, @isscalar);
ip.addParameter('Mode', 'raw', @(x) any(strcmpi(x, {'raw', 'mask'})));
ip.parse(data, varargin{:});

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

% load/generate mask
if strcmpi(ip.Results.Mode, 'mask')
    mask = getCellMask(data);
    if N==1
        mask = {mask};
    end
end

figure;
colormap(gray(256));
for k = 1:N
    y = ceil(k/nx);
    x = k-(y-1)*nx;
    ha = axes('Position', [(dx+wx)*(x-1) 1-wy*y-(y-1)*dy wx wy]);
    
    if ~isempty(ip.Results.Frame)
        if iscell(data(k).framePaths{ip.Results.Channel})
            frame = imread(data(k).framePaths{ip.Results.Channel}{ip.Results.Frame});
        else
            frame = imread(data(k).framePaths{ip.Results.Channel}, ip.Results.Frame);
        end
    else
        %tmp = load([data(k).source 'Detection' filesep 'avgProj.mat']);
        %frame = scaleContrast(tmp.aip);
        %frame = scaleContrast(sqrt(frame));
        frame = scaleContrast(double(imread([data(k).source 'Detection' filesep 'maxproj.tif'])));
    end
    
    imagesc(frame, 'Parent', ha);
    if strcmpi(ip.Results.Mode, 'mask')
        hold on;
        B = bwboundaries(mask{k});
        cellfun(@(i) plot(i(:,2),i(:,1), 'r'), B);
    end
    
    psize = data(k).pixelSize/data(k).M;
    plotScaleBar(5e-6/psize, 'Label', '5 �m');
    
    axis image off;
    r = data(k).imagesize(2)/data(k).imagesize(1);
    d = 0.02;
    text(1-d, 1-d*r, getCellDir(data(k)), 'Color', 'w',...
        'Units', 'normalized',...
        'VerticalAlignment', 'Top',...
        'HorizontalAlignment', 'Right',...
        'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize, 'interpreter', 'none');
    
end