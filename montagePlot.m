% Montage of multichannel time-lapses with adjustable window size.
%
% INPUT:   inputCell : cell array containing image frames
%
% François Aguet, last modified October 16, 2010

function montagePlot(inputCell, varargin)
[nc, nf] = size(inputCell);

phi = (1+sqrt(5))/2;
width = 500;

nx = ceil(sqrt(nc*nf/phi)*phi);
ny = ceil(nf/nx);


% fix dimensions
df = 1/15; % gap between frames, relative to frame width
dc = 3/15; % gap between channels


% width, gaps relative to '1'
wx = 1 / (nx + df*(nx-1));

if nc==1
    wy = 1 / (ny + df*(ny-1));
    % recompute height based on arrangement
    height = width * (ny+(ny-1)*df) / (nx+(nx-1)*df);
    dc = 0;
    dx = wx*df;
    dy = wy*df;
else
    wy = 1 / (ny*nc + df*ny + (ny-1)*dc);
    height = width * (ny*nc + (nc-1)*ny*df + (ny-1)*dc) / (nx+(nx-1)*df);
    dc = wy*dc;
    dx = wx*df;
    dy = wy*df;
end


maxI = zeros(1,nc);
minI = zeros(1,nc);


%ratio = height/width;
%nx = (((1-ratio) + sqrt((ratio-1)^2*df^2+4*nf*nc*(1+df)^2*ratio)) / (2*ratio*(1+df)))
%ny = (nf*nc/nx)
%return

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
                ha(fi) = axes('Position', [(xi-1)*(wx+dx) (ny-yi)*(wy+dy) wx wy]);
                imagesc(inputCell{1, fi}); axis  off; caxis([minI(1) maxI(1)]);
            end
        end
    end
else
    ha = zeros(1,nf*nc);
    for yi = 1:ny
        for ci = 1:nc
            for xi = 1:nx
                % axes index
                fi = xi + (yi-1)*nc*nx + (ci-1)*nx;
                if xi+(yi-1)*nx <= nf
                    ha(fi) = axes('Position', [(xi-1)*(wx+dx) (ny-yi)*(nc*wy + (nc-1)*dy + dc) + (nc-ci)*(wy+dy) wx wy]);
                    frame = scaleContrast(double(inputCell{ci, xi+(yi-1)*nx}), [minI(ci) maxI(ci)]);
                    frame = cat(3, colorV(ci,1)*frame, colorV(ci,2)*frame, colorV(ci,3)*frame);
                    imagesc(uint8(frame)); axis  off;
                end
            end
        end
    end
end

% if (nargin == 3)
%     print('-depsc2', '-loose', '-tiff', fileName);
% end;

    function resizeCallback(src, ~)
%         pos = get(src, 'Position');
%         ratio = pos(4)/pos(3); %H/W
% 
%         nx = ceil(((1-ratio) + sqrt((ratio-1)^2*df^2+4*nf*nc*(1+df)^2*ratio)) / (2*ratio*(1+df)));
%         ny = ceil(nf/nx);       
%                 
%         wx = 1 / (nx + df*(nx-1));
%         if nc==1
%             wy = 1/(nx*(1+df)-df)*pos(3)/pos(4);
%             dc = 0;
%         else
%             %wy = 1 / (ny*nc + df*ny + (ny-1)*dc);
%             wy = (1 - (ny-1)*dc - ny*(nc-1)*dy) / (nc*ny);
%             dc = wy*dc;
%         end
%         
%         dx = wx*df;
%         dy = wy*df;
%         
%         % y offset
%         yOffset = 1-(ny*(1+df)-df)/(nx*(1+df)-df)*pos(3)/pos(4);
%                 
%         if nc==1
%             for yi2 = 1:ny
%                 for xi2 = 1:nx
%                     % frame index
%                     fi = xi2+(yi2-1)*nx;
%                     if fi<=nf
%                         set(ha(fi), 'Position', [(xi2-1)*(wx+dx) yOffset+(ny-yi2)*(wy+dy) wx wy]);
%                     end
%                 end
%             end
%         else
%             for yi2 = 1:ny
%                 for ci2 = 1:nc
%                     for xi2 = 1:nx
%                         % axes index
%                         fi = xi2 + (yi2-1)*nc*nx + (ci2-1)*nx;
%                         if xi2+(yi2-1)*nx <= nf
%                             set(ha(fi), 'Position', [(xi2-1)*(wx+dx) (ny-yi2)*(nc*wy + (nc-1)*dy + dc) + (nc-ci)*(wy+dy) wx wy]);
%                         end
%                     end
%                 end
%             end
%         end 
    end

end