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
ip.addParamValue('Mode', 'raw', @ischar);
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
for k = 1:N
    y = ceil(k/nx);
    x = k-(y-1)*nx;
    h = axes('Position', [(dx+wx)*(x-1) 1-wy*y-(y-1)*dy wx wy]);
    
    if strcmpi(ip.Results.Mode, 'cellmask')
        imagesc(imread([data(k).source 'Detection' filesep 'cellmask.tif']))
        axis image tight;
    elseif strcmpi(ip.Results.Mode, 'cellproj')
        imagesc(imread([data(k).source 'Detection' filesep 'cellAIP.tif']))
        axis image tight;        
    else
        plotFrame(data(k), [], ip.Results.Frame, ip.Results.Channel, 'Handle', h, 'Mode', ip.Results.Mode, 'ScaleBar', 5e-6);
    end
    
    axis image off;
    r = data(k).imagesize(2)/data(k).imagesize(1);
    d = 0.02;
    text(1-d, 1-d*r, getCellDirectory(data(k)), 'Color', 'w',...
                'Units', 'normalized',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Right',...
                'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize, 'interpreter', 'none');

end