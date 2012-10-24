% Francois Aguet, 02/17/2012

function mask = getCellMask(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Connect', true, @islogical);
ip.addParamValue('Display', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('ShowHistogram', false, @islogical);
ip.addParamValue('modeRatio', 0.6, @isscalar);
ip.parse(data, varargin{:});

nd = numel(data);
mask = cell(1,nd);
for i = 1:nd
    aipPath = [data(i).source 'Detection' filesep 'avgProj.mat'];
    if ~(exist(aipPath, 'file')==2)
        df = floor(data(i).movieLength/min(data(i).movieLength,50));
        fv = 1:df:data(i).movieLength;
        nf = numel(fv);
        aip = zeros(data(i).imagesize);
        mCh = strcmpi(data(i).source, data(i).channels);
        for k = 1:nf
            aip = aip + double(imread(data(i).framePaths{mCh}{fv(k)}));
        end
        aip = aip/nf;
        save(aipPath, 'aip');
    else
        load(aipPath);
    end
    
    maskPath = [data(i).source 'Detection' filesep 'cellmask.tif'];
    if ~(exist(maskPath, 'file') == 2) || ip.Results.Overwrite
        mask{i} = computeMask(data(i), aip, ip.Results.Connect, ip.Results.ShowHistogram, ip.Results.modeRatio);
        % save
        imwrite(uint8(mask{i}), maskPath, 'tif', 'compression' , 'lzw');
    else
        %fprintf('Cell mask has already been computed for %s\n', getShortPath(data(i)));
        mask{i} = double(imread(maskPath));
    end
end

if strcmpi(ip.Results.Display, 'on')
    for i = 1:nd
        if ~isempty(mask{i})
            [ny,nx] = size(mask{i});
            B = bwboundaries(mask{i});
            B = cellfun(@(var)sub2ind([ny nx], var(:,1), var(:,2)),B,'UniformOutput',false);
            B = cell2mat(B);
            bmask = zeros([ny nx]);
            bmask(B) = 1;
            bmask = bwmorph(bmask, 'dilate');
            
            aipPath = [data(i).source 'Detection' filesep 'avgProj.mat'];
            load(aipPath);
            aip = scaleContrast(aip);
            aip(bmask==1) = 0;
            overlay = aip;
            overlay(bmask==1) = 255;
            overlay = uint8(cat(3, overlay, aip, aip));
            figure; imagesc(overlay); axis image; colormap(gray(256)); colorbar;
        end
    end
end

if nd==1
    mask = mask{1};
end



function mask = computeMask(data, aip, connect, showHist,modeRatio)

aip = scaleContrast(aip, [], [0 1]);
g = filterGauss2D(aip, 5);

v = aip(:);

% [f_ecdf, x_ecdf] = ecdf(v);
% x_ecdf = x_ecdf(2:end)';
% f_ecdf = f_ecdf(2:end)';

% x1 = interp1(f_ecdf, x_ecdf, 0.99);
pct = prctile(v, [1 99]);
v(v<pct(1) | pct(2)<v) = [];

[f,xi] = ksdensity(v, 'npoints', 100);

% local max/min
lmax = locmax1d(f, 3);
lmin = locmin1d(f, 3);

% max value
% hmax = find(f==max(f), 1, 'first');
dxi = xi(2)-xi(1);

% identify min after first mode
if ~isempty(lmin)
    idx = find(lmin>lmax(1), 1, 'first');
    if ~isempty(idx) && sum(f(1:lmin(idx(1))))*dxi < modeRatio
        min0 = lmin(idx);
        T = xi(min0);
        mask = g>T;
        
        nx = data.imagesize(2);
        ny = data.imagesize(1);
        borderIdx = [1:ny (nx-1)*ny+(1:ny) ny+1:ny:(nx-2)*ny+1 2*ny:ny:(nx-1)*ny];
        
        % retain largest connected component
        if connect
            CC = bwconncomp(mask, 8);
            compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
            mask = zeros(data.imagesize);
            mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
        end
        
%         % fill holes (retain largest boundary)
%         B = bwboundaries(mask);
%         nb = cellfun(@(i) numel(i), B);
%         B = B{nb==max(nb)};
%         boundary = zeros(data.imagesize);
%         boundary(sub2ind(data.imagesize, B(:,1), B(:,2))) = 1;
%         % boundary(borderIdx) = 1;
%         CC = bwconncomp(1-boundary, 4);
%         mask(CC.PixelIdxList{1}) = mode(double(mask(CC.PixelIdxList{1})));
%         mask(CC.PixelIdxList{2}) = mode(double(mask(CC.PixelIdxList{2})));
%         
%         % mask indexes
%         labels = double(labelmatrix(CC));
%         idx = unique(mask.*labels);
%         mask = boundary | ismember(labels, idx(2:end));%labels==idx(2);
           
    else
        mask = ones(data.imagesize);
    end
else
    mask = ones(data.imagesize);
end

if showHist
    dx = xi(2)-xi(1);
    ni = hist(v, xi);
    ni = ni/sum(ni)/dx;
    
    figure;
    hold on;
    plot(xi, ni, 'k.-', 'LineWidth', 3, 'MarkerSize', 20);
    plot(xi, f, 'r-', 'LineWidth', 1);
    set(gca, 'LineWidth', 2, 'FontSize', 18);
end


% figure; imagesc(mask); axis image; colormap(gray(256)); colorbar;

% % vector bounds for projection line
% x0 = [interp1(f_ecdf, x_ecdf, 0.01) 0.01]';
% xN = [interp1(f_ecdf, x_ecdf, 0.99) 0.99]';
%
% p = xN; % proj. vector
% p = p/norm(p);
%
% N = 1e4;
%
% xv = linspace(x0(1), xN(1), N);
% yv = interp1(x_ecdf, f_ecdf, xv);
%
% X = [xv; yv];
%
% xnorm = sqrt(sum(X.^2,1));
% theta = atan2(X(2,:),X(1,:)) - atan2(p(2),p(1));
% ynorm = sin(theta).*xnorm;
%
% T = X(1,find(ynorm==min(ynorm), 1, 'first'));
%
% mask = bf>T;
%
% % close with PSF
% mCh = strcmpi(data.source, data.channels);
% if isempty(sigma)
%     sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});
% end
% w = ceil(4*sigma);
% se = strel('disk', w, 0);
% mask = imclose(mask, se);
%
% % retain largest connected component
% CC = bwconncomp(mask, 8);
% compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
% mask = zeros(data.imagesize);
% mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;

% figure; imagesc(mask); axis image; colormap(gray(256)); colorbar;

% load([data.source 'Detection' filesep 'detection_v2.mat']);
% ny = data.imagesize(1);
% nx = data.imagesize(2);
%
%
%
%
% % concatenate all positions
% X = [frameInfo.x];
% X = X(mCh,:);
% Y = [frameInfo.y];
% Y = Y(mCh,:);
%
% mask = zeros(data.imagesize);
% mask(sub2ind(data.imagesize, round(Y), round(X))) = 1;
%
% se = strel('disk', w, 0);
% mask = imclose(mask, se);
% mask = bwmorph(mask, 'clean');
% mask = imdilate(mask, se);
%
%
% % retain largest connected component
% CC = bwconncomp(mask, 8);
% compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
% mask = zeros(ny,nx);
% mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
%
% % add border within 'w'
% [yi,xi] = ind2sub([ny nx], find(mask==1));
% [yb,xb] = ind2sub([ny nx], borderIdx);
% idx = KDTreeBallQuery([xi yi], [xb' yb'], w);
% mask(borderIdx(cellfun(@(x) ~isempty(x), idx))) = 1;
%
% mask = imclose(mask, se);
%
% % fill holes
% CC = bwconncomp(~mask, 8);
% M = labelmatrix(CC);
%
% % hole labels in image border
% borderLabels = unique(M(borderIdx));
% labels = setdiff(1:CC.NumObjects, borderLabels);
% idx = vertcat(CC.PixelIdxList{labels});
% mask(idx) = 1;
