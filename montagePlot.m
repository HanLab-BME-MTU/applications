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



% figure size
wx = 1;
dx = 1/15; % gap between frames, relative to frame width
dc = 1/5; % gap between channels

nr = 1:5;

height = nr*(nc*wx+(nc-1)*dx) + (nr-1)*dc;
nx = ceil(nf./nr);
width = nx*wx + (nx-1)*dx;

idx = (width./height - (1+sqrt(5))/2).^2;
idx = idx == min(idx);

% final values
width = width(idx);
height = height(idx);
nx = nx(idx);
nr = nr(idx);

ratio = height/width;

% initial width/height
width = 500;
height = ratio*width;
wxi = width / ((nx/dx + nx-1)*dx);
dxi = wxi*dx;
dci = wxi*dc;




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
    offset = maxTextWidth + dci;
else
    offset = 0;
end
delete(ha);


set(hf, 'Position', [50, 200, width+offset, height], 'Visible', visible);
for rowi = 1:nr
    for c = 1:nc
        for x = 1:nx
            fi = x+(rowi-1)*nx;
            if fi<=nf
                axes('Units', 'pixels',...
                    'Position', [offset+(x-1)*(wxi+dxi) height-wxi - ((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
                    'XLim', [0 wxi], 'YLim', [0 wxi]);
                imagesc([0 wxi], [0 wxi], inputCell{c, fi}); axis off; caxis([0 255]);%caxis([minI(c) maxI(c)]);
                hold on;
                if x==1 && rowi==1 && ~isempty(labels)
                    ht = text(-dci, wxi/2, labels{c}, 'Units', 'pixels',...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline',...
                        'FontName', fontName, 'FontSize', fontSize, 'Color', rgbColors{c});
                    get(ht, 'extent');
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


% % [nx ny wx wy dx dy offset dc] = getLayout(width, height, nf, nc);
% 
% 
% 
% 
% figure('ResizeFcn', @resizeCallback, 'Position', [50, 200, width, height]);
% colormap(gray(256));
% if nc==1
%     ha = zeros(1,nf);
%     for yi = 1:ny
%         for xi = 1:nx
%             % frame index
%             fi = xi+(yi-1)*nx;
%             if fi<=nf
%                 ha(fi) = axes('Position', [(xi-1)*(wx+dx) offset+(ny-yi)*(wy+dy) wx wy]);
%                 imagesc(inputCell{1, fi}); axis off; caxis([minI(1) maxI(1)]);
%             end
%         end
%     end
% else
%     ha = zeros(nc,nf);
%     for yi = 1:ny
%         for ci = 1:nc
%             for xi = 1:nx
%                 % axes index
%                 fi = xi + (yi-1)*nx;
%                 if xi+(yi-1)*nx <= nf
%                     ha(ci,fi) = axes('Position', [(xi-1)*(wx+dx) offset+(ny-yi)*(nc*wy+(nc-1)*dy+dc) + (nc-ci)*(wy+dy) wx wy]);
%                     frame = scaleContrast(double(inputCell{ci, xi+(yi-1)*nx}), [minI(ci) maxI(ci)]);
%                     frame = cat(3, colorV(ci,1)*frame, colorV(ci,2)*frame, colorV(ci,3)*frame);
%                     imagesc(uint8(frame)); axis off;
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% %     function resizeCallback(src, ~)
% %         pos = get(src, 'Position');
% %         [nx ny wx wy dx dy offset dc] = getLayout(pos(3), pos(4), nf, nc);
% %         
% %         if nc==1
% %             for yi2 = 1:ny
% %                 for xi2 = 1:nx
% %                     % frame index
% %                     fi = xi2+(yi2-1)*nx;
% %                     if fi<=nf
% %                         set(ha(fi), 'Position', [(xi2-1)*(wx+dx) (ny-yi2)*(wy+dy) wx wy]);
% %                     end
% %                 end
% %             end
% %         else
% %             for yi2 = 1:ny
% %                 for ci2 = 1:nc
% %                     for xi2 = 1:nx
% %                         fi = xi2 + (yi2-1)*nx;
% %                         if fi<=nf
% %                             set(ha(ci2,fi), 'Position', [(xi2-1)*(wx+dx) offset+(ny-yi2)*(nc*wy + (nc-1)*dy + dc) + (nc-ci2)*(wy+dy) wx wy]);
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end
% end
% 
% 
% function [nx ny wx wy dx dy offset dc] = getLayout(width, height, nf, nc)
% df = 1/15; % gap between frames, relative to frame width
% dc = 3/15; % gap between channels
% 
% nxVect = 1:nf*nc;
% wVect = width ./ (nxVect + df*(nxVect-1));
% 
% if nc==1
%     nyVect = floor((height+df*wVect)./(wVect*(1+df)));
%     
% else
%     nyVect = floor((height+dc*wVect)./(wVect*(nc + df*(nc-1) + dc)));
% end
% ntot = nxVect .* nyVect*nc;
% idx = find(ntot > nf*nc, 1, 'first');
% 
% nx = nxVect(idx);
% ny = nyVect(idx);
% 
% % width of frame
% w = width / (nx + df*(nx-1));
% 
% % fraction of height used
% if nc==1
%     p = w*(ny + df*(ny-1))/height;
% else
%     p = w*(ny*nc + dc*(ny-1) + df*(nc-1)*ny)/height;
% end
% 
% % proportions relative to 1
% wx = w/width;
% 
% if nc==1
%     wy = p / (ny + df*(ny-1));
%     dc = 0;
% else
%     wy = p / (ny*nc + dc*(ny-1) + df*(nc-1)*ny);
%     dc = wy*dc;
% end
% 
% dx = wx*df;
% dy = wy*df;
% offset = 1-p;
% 
% end
