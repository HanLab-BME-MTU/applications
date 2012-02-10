% Francois Aguet, 02/08/2012

function [mask] = getCellMask(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Sigma', []);
ip.addParamValue('Display', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
ip.parse(data, varargin{:});

for i = 1:length(data)
    if ~(exist([data(i).source 'Detection' filesep 'detection_v2.mat'], 'file')==2)
        error('Detection must be run before using this function.');
    elseif ~(exist([data(i).source 'Detection' filesep 'cellmask.tif'], 'file') == 2) || ip.Results.Overwrite
        mask = computeMask(data(i), ip.Results.Sigma);
        
        if strcmpi(ip.Results.Display, 'on')
            [ny,nx] = size(mask);
            B = bwboundaries(mask);
            B = sub2ind([ny nx], B{1}(:,1), B{1}(:,2));
            bmask = zeros([ny nx]);
            bmask(B) = 1;
            bmask = bwmorph(bmask, 'dilate');
            
            frame1 = scaleContrast(double(imread(data.framePaths{1}{1})));
            frame1(bmask==1) = 0;
            overlay = frame1;
            overlay(bmask==1) = 255;
            overlay = uint8(cat(3, overlay, frame1, frame1));
            figure; imagesc(overlay); axis image; colormap(gray(256)); colorbar;
        end
        % save
        imwrite(uint8(255*mask), [data.source 'Detection' filesep 'cellmask.tif'], 'tif', 'compression' , 'lzw');
    else
        fprintf('Cell mask has already been computed for %s\n', getShortPath(data(i)));
    end
end


function mask = computeMask(data, sigma)
load([data.source 'Detection' filesep 'detection_v2.mat']);
mCh = find(strcmpi(data.source, data.channels));

ny = data.imagesize(1);
nx = data.imagesize(2);
borderIdx = [1:ny (nx-1)*ny+(1:ny) ny+1:ny:(nx-2)*ny+1 2*ny:ny:(nx-1)*ny];

% concatenate all positions
X = [frameInfo.x];
X = X(mCh,:);
Y = [frameInfo.y];
Y = Y(mCh,:);

mask = zeros(data.imagesize);
mask(sub2ind(data.imagesize, round(Y), round(X))) = 1;

if isempty(sigma)
    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});
end
w = ceil(4*sigma);

se = strel('disk', w, 0);
mask = imclose(mask, se);

mask = bwmorph(mask, 'clean');
CC = bwconncomp(mask, 8);
compsize = cellfun(@(i) numel(i), CC.PixelIdxList);

% retain largest component
mask = zeros(ny,nx);
mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
mask = imdilate(mask, se);

% add border within 'w'
% [yi,xi] = ind2sub([ny nx], find(mask==1));
% [yb,xb] = ind2sub([ny nx], borderIdx);
% idx = KDTreeBallQuery([xi yi], [xb' yb'], w);
% mask(borderIdx(cellfun(@(x) ~isempty(x), idx))) = 1;

mask = imclose(mask, se);

% fill holes
CC = bwconncomp(~mask, 8);
M = labelmatrix(CC);

% hole labels in image border
borderLabels = unique(M(borderIdx));
labels = setdiff(1:CC.NumObjects, borderLabels);
idx = vertcat(CC.PixelIdxList{labels});
mask(idx) = 1;

% average intensity projection
aip = zeros(ny,nx);
nf = 0;
for k = 1:10:data.movieLength
    aip = aip + double(imread(data.framePaths{mCh}{k}));
    nf = nf+1;
end
aip = aip / nf;

% smooth projection
% sigma2 = 2*sigma;
% w = ceil(4*sigma2);
% g = exp(-(-w:w).^2/(2*sigma2^2));
% g = g/sum(g);
% aip = conv2(g, g, padarrayXT(aip, [w w], 'Symmetric'), 'valid');

mf = 1/norminv(0.75);

bgMean = median(aip(mask==0));
fgMean = median(aip(mask==1));
bgVar = (mf * mad(aip(mask==0),1))^2;
fgVar = (mf * mad(aip(mask==1),1))^2;

obj1 = gmdistribution.fit(aip(:), 1, 'CovType', 'diagonal');
S.mu = [bgMean fgMean]';
S.Sigma = reshape([bgVar fgVar], [1 1 2]);
obj2 = gmdistribution.fit(aip(:), 2, 'start', S);

if obj2.BIC < obj1.BIC % threshold found
%     mu1 = obj2.mu(1);
%     mu2 = obj2.mu(2);
    stdV = sqrt(obj2.Sigma);
%     s1 = stdV(1);
%     s2 = stdV(2);
%     a1 = obj2.PComponents(1);
%     a2 = obj2.PComponents(2);
%     
%     a = s2^2-s1^2;
%     b = 2.0*(mu2*s1^2 - mu1*s2^2);
%     c = mu1^2*s2^2 - mu2^2*s1^2 - 2.0*log(a1/a2);
%     
%     T = (-b + sqrt(b^2-4.0*a*c))/(2*a);
    T = norminv(0.99, obj2.mu(1), stdV(1));
%     T = norminv(0.05, obj2.mu(2), stdV(2));
    mask = aip>T;
    mask = bwmorph(mask, 'clean');
    mask = medfilt2(mask,[9 9]);
    se = strel('disk', w, 0);
    mask = imclose(mask, se);
    
    % close any holes not connected to border
    CC = bwconncomp(~mask, 8);
    M = labelmatrix(CC);
    borderLabels = unique(M(borderIdx));
    labels = setdiff(1:CC.NumObjects, borderLabels);
    idx = vertcat(CC.PixelIdxList{labels});
    mask(idx) = 1;
    
    % keep largest connected component
    CC = bwconncomp(mask, 8);
    mask = zeros(ny,nx);
    mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
    
else
    mask = ones(ny,nx);
end
