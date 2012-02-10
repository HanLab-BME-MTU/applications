function getBackgroundIntensity(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Sigma', []);
% ip.addParamValue('Visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
% ip.addParamValue('Mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.parse(data, varargin{:});
sigma = ip.Results.Sigma;

mCh = find(strcmpi(data.channels, data.source));

load([data.source 'Detection' filesep 'detection_v2.mat']);
mask = double(imread([data.source 'Detection' filesep 'cellmask.tif']));

ny = data.imagesize(1);
nx = data.imagesize(2);

if isempty(sigma)
    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});
end

w = ceil(4*sigma);
se = strel('disk', 2*w, 0);
mask = imerode(mask, se);

xi = round(frameInfo(1).x(mCh,:));
yi = round(frameInfo(1).y(mCh,:));

posMask = zeros(ny,nx);
posMask(sub2ind([ny nx], yi,xi)) = 1;
posMask = posMask & mask;
mask = mask - (mask & imdilate(posMask, se));

mask([1:1+w end-w:end],:) = 0;
mask(:,[1:1+w end-w:end]) = 0;

% figure; imagesc(mask); axis image; colormap(gray(256)); colorbar;

nf = data.movieLength;
fi = ceil([1 nf/4 nf/2 3*nf/4 nf]);
for f = 1:numel(fi)
    
    % draw 1000 points in mask
    N = 1e3;
    xacc = [];
    yacc = [];
    while numel(xacc)<N
        xcand = (nx-1)*rand(1,N)+1;
        ycand = (ny-1)*rand(1,N)+1;
        idx = mask(sub2ind([ny,nx], round(ycand), round(xcand)));
        xacc = [xacc xcand(idx==1)]; %#ok<AGROW>
        yacc = [yacc ycand(idx==1)]; %#ok<AGROW>
    end
    xacc = xacc(1:N);
    yacc = yacc(1:N);
    
    
    frame = double(imread(data.framePaths{mCh}{fi(f)}));
    
    pStruct(f) = fitGaussians2D(frame, xacc, yacc, [], sigma, [], 'Ac');
end



A = arrayfun(@(s) mean(s.A), pStruct);
A_pstd = arrayfun(@(s) mean(s.A_pstd), pStruct);

stdA = arrayfun(@(s) std(s.A), pStruct);
stdA_pstd = arrayfun(@(s) std(s.A_pstd), pStruct);


figure;
hold on;
plot(fi, A, 'k.', 'MarkerSize', 24);
he = errorbar(fi, A, A_pstd, 'k', 'LineWidth', 2, 'LineStyle', 'none');

he = errorbar(fi, A, stdA, 'r', 'LineWidth', 2, 'LineStyle', 'none');

he = errorbar(fi, A+A_pstd, stdA_pstd, 'b', 'LineWidth', 2, 'LineStyle', 'none');
he = errorbar(fi, A-A_pstd, stdA_pstd, 'b', 'LineWidth', 2, 'LineStyle', 'none');

% setErrorbarStyle(he);
set(gca, 'XTick', fi);
xlabel('Frame #');

