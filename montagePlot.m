% Montage of multichannel time-lapses with adjustable window size.
%
% INPUT:   inputCell : cell array containing image frames
%
% François Aguet, last modified October 19, 2010

function montagePlot(inputCell, varargin)
[nc, nf] = size(inputCell);

phi = (1+sqrt(5))/2;
width = 500;
height = floor(width/phi);

[nx ny wx wy dx dy offset dc] = getLayout(width, height, nf, nc);

maxI = zeros(1,nc);
minI = zeros(1,nc);
for k = 1:nc
    cCat = [inputCell{k,:}];
    maxI(k) = max(cCat(:));
    minI(k) = min(cCat(:));
end

if nargin>1
    colorV = varargin{1};
else
    colorV = ones(nc,3);
end

figure('ResizeFcn', @resizeCallback, 'Position', [50, 200, width, height], 'PaperPositionMode', 'auto');
colormap(gray(256));
if nc==1
    ha = zeros(1,nf);
    for yi = 1:ny
        for xi = 1:nx
            % frame index
            fi = xi+(yi-1)*nx;
            if fi<=nf
                ha(fi) = axes('Position', [(xi-1)*(wx+dx) offset+(ny-yi)*(wy+dy) wx wy]);
                imagesc(inputCell{1, fi}); axis off; caxis([minI(1) maxI(1)]);
            end
        end
    end
else
    ha = zeros(nc,nf);
    for yi = 1:ny
        for ci = 1:nc
            for xi = 1:nx
                % axes index
                fi = xi + (yi-1)*nx;
                if xi+(yi-1)*nx <= nf
                    ha(ci,fi) = axes('Position', [(xi-1)*(wx+dx) offset+(ny-yi)*(nc*wy+(nc-1)*dy+dc) + (nc-ci)*(wy+dy) wx wy]);
                    frame = scaleContrast(double(inputCell{ci, xi+(yi-1)*nx}), [minI(ci) maxI(ci)]);
                    frame = cat(3, colorV(ci,1)*frame, colorV(ci,2)*frame, colorV(ci,3)*frame);
                    imagesc(uint8(frame)); axis off;
                end
            end
        end
    end
end

% if (nargin == 3)
%     print('-depsc2', '-loose', '-tiff', fileName);
% end;

    function resizeCallback(src, ~)
        pos = get(src, 'Position');
        [nx ny wx wy dx dy offset dc] = getLayout(pos(3), pos(4), nf, nc);
        
        if nc==1
            for yi2 = 1:ny
                for xi2 = 1:nx
                    % frame index
                    fi = xi2+(yi2-1)*nx;
                    if fi<=nf
                        set(ha(fi), 'Position', [(xi2-1)*(wx+dx) (ny-yi2)*(wy+dy) wx wy]);
                    end
                end
            end
        else
            for yi2 = 1:ny
                for ci2 = 1:nc
                    for xi2 = 1:nx
                        fi = xi2 + (yi2-1)*nx;
                        if fi<=nf
                            set(ha(ci2,fi), 'Position', [(xi2-1)*(wx+dx) offset+(ny-yi2)*(nc*wy + (nc-1)*dy + dc) + (nc-ci2)*(wy+dy) wx wy]);
                        end
                    end
                end
            end
        end
    end
end


function [nx ny wx wy dx dy offset dc] = getLayout(width, height, nf, nc)
df = 1/15; % gap between frames, relative to frame width
dc = 3/15; % gap between channels

nxVect = 1:nf*nc;
wVect = width ./ (nxVect + df*(nxVect-1));

if nc==1
    nyVect = floor((height+df*wVect)./(wVect*(1+df)));
    
else
    nyVect = floor((height+dc*wVect)./(wVect*(nc + df*(nc-1) + dc)));
end
ntot = nxVect .* nyVect*nc;
idx = find(ntot > nf*nc, 1, 'first');

nx = nxVect(idx);
ny = nyVect(idx);

% width of frame
w = width / (nx + df*(nx-1));

% fraction of height used
if nc==1
    p = w*(ny + df*(ny-1))/height;
else
    p = w*(ny*nc + dc*(ny-1) + df*(nc-1)*ny)/height;
end

% proportions relative to 1
wx = w/width;

if nc==1
    wy = p / (ny + df*(ny-1));
    dc = 0;
else
    wy = p / (ny*nc + dc*(ny-1) + df*(nc-1)*ny);
    dc = wy*dc;
end

dx = wx*df;
dy = wy*df;
offset = 1-p;

end
