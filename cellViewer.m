function cellViewer(data, frameIdx)

if nargin<2 || isempty(frameIdx)
    frameIdx = 1;
end

N = numel(data);


ny = ceil(sqrt(2*N/(1+sqrt(5))));
nx = ceil(N/ny);

dx = 0.1/(nx-1);
dy = 0.1/(ny-1);

wx = 0.9/nx;
wy = 0.9/ny;

figure;
for k = 1:N
    y = ceil(k/nx);
    x = k-(y-1)*nx;
    h = axes('Position', [(dx+wx)*(x-1) 1-wy*y-(y-1)*dy wx wy]);
    frameList = dir([data(k).source '*.tif*']);
    frame = imread([data(k).source frameList(frameIdx).name]);
    [sy sx] = size(frame);
    psize = data(k).pixelSize/data(k).M*1e6; % in µm
    xa = (0:sx-1)*psize;
    ya = (0:sy-1)*psize;
    if sx>sy
        imagesc(xa,ya,frame);
    else
        imagesc(ya,xa,imrotate(frame, 90));
    end
    colormap(gray(256));
    axis image off;
    title(getCellDirectory(data(k)), 'interpreter', 'none');
    plotScaleBar(5, h, 'Label', '5 µm', 'FontSize', 8);

end