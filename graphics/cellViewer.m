function cellViewer(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Frame', 1, @isscalar);
ip.addParamValue('Channel', strcmp(data(1).channels, data(1).source), @isscalar);
ip.addParamValue('Scale', 5, @isscalar);
ip.addParamValue('Units', 'µm', @ischar);
ip.addParamValue('FontName', 'Helvetica', @ischar);
ip.addParamValue('FontSize', 10, @isscalar);
ip.addParamValue('Mode', 'raw', @(x) any(strcmpi(x, {'raw', 'mask'})));
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

figure;
colormap(gray(256));
for k = 1:N
    y = ceil(k/nx);
    x = k-(y-1)*nx;
    ha = axes('Position', [(dx+wx)*(x-1) 1-wy*y-(y-1)*dy wx wy]);

    %frame = scaleContrast(double(imread(data(k).framePaths{ip.Results.Channel}{ip.Results.Frame})));
    frame = double(imread([data(k).source 'Detection' filesep 'avgProj.tif']));
    if strcmpi(ip.Results.Mode, 'mask')
        mask = double(imread([data(k).source 'Detection' filesep 'cellmask.tif']));
        
        B = bwboundaries(mask);
        B = sub2ind(size(mask), B{1}(:,1), B{1}(:,2));
        mask = zeros(size(mask));
        mask(B) = 1;
        mask = bwmorph(mask, 'dilate');

        frame(mask==1) = 0;
        overlay = frame;
        overlay(mask==1) = 255;
        frame = uint8(cat(3, overlay, frame, frame));
    end
    
    imagesc(frame, 'Parent', ha);
    psize = data(k).pixelSize/data(k).M;
    plotScaleBar(5e-6/psize, 'Label', '5 µm');

    axis image off;
    r = data(k).imagesize(2)/data(k).imagesize(1);
    d = 0.02;
    text(1-d, 1-d*r, getCellDir(data(k)), 'Color', 'w',...
                'Units', 'normalized',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Right',...
                'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize, 'interpreter', 'none');

end