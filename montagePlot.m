% Montage of multichannel time-lapses with adjustable window size.
%
% INPUT:   inputCell : cell array containing image frames
%
% François Aguet, last modified March 22, 2011

function montagePlot(inputCell, varargin)

[nc, nf] = size(inputCell);

if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

%======================================
% Parse inputs, set defaults
%======================================
% idx = find(strcmpi(varargin, 'Color'));
% if ~isempty(idx)
%     colorV = varargin{idx+1};
% else
%     colorV = ones(nc,3);
% end

idx = find(strcmpi(varargin, 'Labels'));
if ~isempty(idx)
    labels = varargin{idx+1};
    if length(labels)~=nc
        error('The number of labels must match the number of channels.');
    end
else
    labels = [];
end

if ~isempty(labels) && sum(isnan(name2wavelength(labels)))==0
    rgbColors = arrayfun(@(x) hsv2rgb([x 1 0.9]), getHuesFromMarkers(labels), 'UniformOutput', false);
else
    rgbColors = [];
end

% dynamic range adjustments
maxI = zeros(1,nc);
minI = zeros(1,nc);
for k = 1:nc
    cCat = [inputCell{k,:}];
    maxI(k) = max(cCat(:));
    minI(k) = min(cCat(:));
end

% loop through cells, adjust contrast, convert to 8-bit
for c = 1:nc
    inputCell(c,:) = cellfun(@(x) uint8(scaleContrast(x, [minI(c) maxI(c)])), inputCell(c,:), 'UniformOutput', false);
end

idx = find(strcmpi(varargin, 'Mode'));
if ~isempty(idx)
    rgbIdx = getRGBindex(labels);
    [ny,nx] = size(inputCell{1,1});
    rgbCells = cell(3,nf);
    rgbCells(rgbIdx,:) = inputCell;
    zcells = squeeze(num2cell(zeros(ny,nx,nf,'uint8'), [1 2]));
    for k = setdiff(1:3, rgbIdx)
        rgbCells(k,:) = zcells;
    end
    rgbCells = arrayfun(@(x) cat(3,rgbCells{:,x}), 1:nf, 'UniformOutput', false);
    
    switch varargin{idx+1}
        case 'color'
            for c = 1:nc
                inputCell(c,:) = cellfun(@(x) circshift(cat(3, x, zeros(ny,nx,2)), [0 0 rgbIdx(c)-1]), inputCell(c,:), 'UniformOutput', false);
            end
            rgbColors = repmat({zeros(1,3)}, [1 nc]);
        case 'RGB+color'
            for c = 1:nc
                inputCell(c,:) = cellfun(@(x) circshift(cat(3, x, zeros(ny,nx,2)), [0 0 rgbIdx(c)-1]), inputCell(c,:), 'UniformOutput', false);
            end
            inputCell = cat(1, rgbCells, inputCell);
            nc = nc + 1;
            labels = [' ', labels];
            rgbColors = [zeros(1,3) rgbColors];
        case 'RGB+gray'
            inputCell = cat(1, rgbCells, inputCell);
            nc = nc + 1;
            labels = [' ', labels];
            rgbColors = [zeros(1,3) rgbColors];
        case 'RGB'
            inputCell = rgbCells;
            nc = 1;
            labels = [];
    end
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
    fontSize = 14;
end

idx = find(strcmpi(varargin, 'Visible'));
if ~isempty(idx)
    visible = varargin{idx+1};
else
    visible = 'on';
end

idx = find(strcmpi(varargin, 'Print'));
if ~isempty(idx)
    printPath = varargin{idx+1};
else
    printPath = [];
end


height = 250;
width = height * (1+sqrt(5))/2;
[wxi, dxi, dci, nx, nr, width] = getProportions(width, height, nf, nc);


hf = figure('Visible', 'off', 'PaperPositionMode', 'auto', 'Position', [50, 200, width, height]);
ha = axes('Position', [0 0 1 1], 'XLim', [0 width], 'YLim', [0 height]);
set(hf, 'Units', 'pixels');

% compute max. text length
textWidth = zeros(1,nc);
if ~isempty(labels)
    for c = 1:nc
        ht = text(0, 0, labels{c}, 'FontName', fontName, 'FontSize', fontSize);
        extent = get(ht, 'extent');
        textWidth(c) = extent(3);
        delete(ht);
    end
    maxTextWidth = max(textWidth);
    offset = maxTextWidth + 2*dci;
else
    offset = 0;
end
delete(ha);

ha = zeros(nc,nf);
set(hf, 'Position', [50, 100, width+offset, height], 'Visible', visible, 'ResizeFcn', {@resizeCallback});
for rowi = 1:nr
    for c = 1:nc
        for x = 1:nx
            fi = x+(rowi-1)*nx;
            if fi<=nf
                ha(c,fi) = axes('Units', 'pixels',...
                    'Position', [offset+(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
                    'XLim', [0 wxi], 'YLim', [0 wxi]);
                imagesc([0 wxi], [0 wxi], inputCell{c, fi}); axis off; caxis([0 255]);%caxis([minI(c) maxI(c)]);
                if x==1 && rowi==1 && ~isempty(labels)
                    ht(c) = text(-dci, wxi/2, labels{c}, 'Units', 'pixels',...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline',...
                        'FontUnits', 'pixels', 'FontName', fontName, 'FontSize', wxi/2.5, 'Color', rgbColors{c});
                end
            end
        end
    end
end
colormap(gray(256));


if ~isempty(printPath)
    print(hf, '-depsc2', printPath);
end

if strcmp(visible, 'off')
    close(hf);
end

    function resizeCallback(src, ~)
        [nc,nf] = size(inputCell);
        pos = get(src, 'Position');
        width = pos(3);
        height = pos(4);
        
        [wxi, dxi, dci, nx, nr] = getProportions(width-offset, height, nf, nc);
        
        for rowi = 1:nr
            for c = 1:nc
                for x = 1:nx
                    fi = x+(rowi-1)*nx;
                    if fi<=nf
                        set(ha(c, fi),...
                            'Position', [offset+(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi]);
                        
                    end
                end
                set(ht(c), 'Position', [-dci, wxi/2], 'FontSize', wxi/2.5);
            end
        end
    end

end



function [wxi, dxi, dci, nx, nr, width] = getProportions(width, height, nf, nc)

% fixed proportions
wx = 1;
dx = 1/15; % gap between frames, relative to frame width
dc = 1/5; % gap between channels

hwRatio = height/width;

nr = 1:20;

% find arrangement that comes closest to < ratio
iheight = nr*(nc*wx+(nc-1)*dx) + (nr-1)*dc;
nx = ceil(nf./nr);
iwidth = nx*wx + (nx-1)*dx;

% distance
idx = (iheight./iwidth - hwRatio).^2;
idx = idx == min(idx);

% corrected width/height:
ratio = iheight(idx)/iwidth(idx);
nr = nr(idx);
nx = nx(idx);

% initial width/height
width = height/ratio;
wxi = width / ((nx/dx + nx-1)*dx);
dxi = wxi*dx;
dci = wxi*dc;

end


