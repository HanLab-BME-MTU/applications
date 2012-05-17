function plotTrackXY(data, track, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('splitIdx', []);
ip.addParamValue('Legend', 'hide', @(x) any(strcmpi(x, {'show','hide'})));
ip.addParamValue('Units', 'pixels', @(x) any(strcmpi(x, {'um','pixels'})));
ip.parse(varargin{:});


mCh = find(strcmp(data.channels, data.source));

x0 = nanmean(track.x(mCh,:));
y0 = nanmean(track.y(mCh,:));

x = track.x(mCh,:)-x0;
y = track.y(mCh,:)-y0;
np = numel(x);

switch ip.Results.Units
    case 'pixels'
        %w = ceil(max([-min(x) max(x) -min(y) max(y)])/2)*2;
        w = 10;
        xa = -w:2:w;
    case 'um'
        x = x*data.pixelSize/data.M*1e6;
        y = y*data.pixelSize/data.M*1e6;
        max([-min(x) max(x) -min(y) max(y)])
        w = 0.6;
        xa = -w:0.2:w;
end


figure;
hold on;

cmap = jet(np);
colormap(cmap);
mesh([x(:) x(:)], [y(:) y(:)], zeros(np,2), repmat((1:np)', [1 2]),...
    'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.5);

gapIdx = find(track.gapVect~=0);
for k = 1:numel(gapIdx)
    plot(x(gapIdx(k)), y(gapIdx(k)), 'o', 'Color', cmap(gapIdx(k),:), 'MarkerFaceColor', 0.999*[1 1 1], 'LineWidth', 1.5);
end

axis([xa(1) xa(end) xa(1) xa(end)]);
axis square;

fset = loadFigureSettings();
set(gca, fset.axOpts{:}, 'XTick', xa, 'YTick', xa);
xlabel(['x (' ip.Results.Units ')'], fset.lfont{:});
ylabel(['y (' ip.Results.Units ')'], fset.lfont{:});

hc = colorbar;
set(get(hc, 'YLabel'), 'String', 'Time (s)', 'FontSize', 16);
set(hc, 'TickLength', [0 0]);

if strcmpi(ip.Results.Legend, 'show')
    h = plot(0,0, 'ko-', 'MarkerFaceColor', 0.999*[1 1 1], 'LineWidth', 1.5);
    hl = legend(h, 'gap');
    set(hl, 'Box', 'off', fset.tfont{:});
    delete(h);
end

if ~isempty(ip.Results.splitIdx)
    splitIdx = ip.Results.splitIdx;

    plot(x(splitIdx), y(splitIdx), 'o', 'Color', cmap(splitIdx,:), 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5);

    x1 = x(mCh, 1:splitIdx-1);
    y1 = y(mCh, 1:splitIdx-1);
    x2 = x(mCh, splitIdx+1:end);
    y2 = y(mCh, splitIdx+1:end);
    mux1 = median(x1);
    muy1 = median(y1);
    mux2 = median(x2);
    muy2 = median(y2);
    
    % projections
    v = [mux2-mux1 muy2-muy1];
    %v(2) = v(1); % test final projection
    v = v/norm(v);
    
    % scalar products
    sp1 = sum(repmat(v', [1 numel(x1)]) .* [x1; y1], 1);
    sp2 = sum(repmat(v', [1 numel(x2)]) .* [x2; y2], 1);
    
    x1proj = v(1)*sp1;
    y1proj = v(2)*sp1;
    
    x2proj = v(1)*sp2;
    y2proj = v(2)*sp2;
   
    
    plot(x1proj, y1proj, 'r.');
    plot(x2proj, y2proj, 'b.');
    
    plot(mux1, muy1, 'rx', 'MarkerSize', 20, 'LineWidth', 2);
    plot(mux2, muy2, 'bx', 'MarkerSize', 20, 'LineWidth', 2);
    
    % vectors on connecting line
    V1 = [x1proj-mux1; y1proj-muy1];
    % sign: scalar product with 'v'
    N1 = sqrt(sum(V1.^2,1));
    d1 = mux1 - N1 .* sum(V1./repmat(N1,[2 1]) .* repmat(v', [1 numel(x1)]), 1);

    V2 = [x2proj-mux2; y2proj-muy2];
    N2 = sqrt(sum(V2.^2,1));
    d2 = mux2 - N2 .* sum(V2./repmat(N2,[2 1]) .* repmat(v', [1 numel(x2)]), 1);

    plot(d1, zeros(size(d1)), 'm.');
    plot(d2, zeros(size(d2)), 'c.');
end

