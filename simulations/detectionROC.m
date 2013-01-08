



nx = 500;
ny = 500;
sigma = 1.4;
w = ceil(4*sigma);

np = 1000;

x = [];
y = [];
while numel(x)<np
    xcand = (nx-2*w-1)*rand(np,1) + w+1;
    ycand = (ny-2*w-1)*rand(np,1) + w+1;
    idx = KDTreeBallQuery([xcand ycand; x y], [xcand ycand], 2*w);
    idx(cellfun(@numel, idx)>1) = [];
    idx = vertcat(idx{:});
    x = [x; xcand(idx)]; %#ok<AGROW>
    y = [y; ycand(idx)]; %#ok<AGROW>
end
x = x(1:np);
y = y(1:np);

% generate point source image
A = 1;
ni = 2*w+1; % width of PSF support
img = simGaussianSpots(nx, ny, sigma, 'x', x, 'y', y);

figure;
imagesc(img); colormap(gray(256)); axis image;
hold on;
plot(x, y, 'ro');
%%

% loop through PSNR values
psnrV = 0:1:25;
N = numel(psnrV);

psnrV = 15;
i = 1;

sigma_n = sqrt(A^2 * ni / (ni-1) / 10^(psnrV(i)/10));
img_n = img + sigma_n*randn(ny,nx);

[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img_n, sigma, 'Alpha', 0.2);
[ly,lx] = find(imgLM~=0);


figure;
imagesc(img_n); colormap(gray(256)); axis image;
hold on;
plot(x, y, 'go');

%%

sigma_n = sqrt(1^2 * ni / (ni-1) / 10^(10/10));

av = [0.5:0.5:5 6:20];
N = numel(av);
meanStd = zeros(N,1);
meanA = zeros(N,1);
parfor i = 1:N
    img_n = av(i)*img + sigma_n*randn(ny,nx);
    pstruct = pointSourceDetection(img_n, sigma);
    meanStd(i) = median(pstruct.A_pstd);
    meanA(i) = median(pstruct.A);
end

%%
figure;
hold on;
plot(meanA, meanStd, 'ko-');
% plot(meanA,%);
xlabel('A');
ylabel('\sigma_A')


%%
img_n = randn(ny,nx);

[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img_n, sigma, 'Alpha', 0.2);
[ly,lx] = find(imgLM~=0);


if ~isempty(pstruct)
    figure;
    imagesc(img_n); colormap(gray(256)); axis image;
    hold on;
    plot(lx, ly, 'go');
    plot(pstruct.x, pstruct.y, 'rx');
end

%%
% Fisher information matrix (parameters: x,y,A,c)
A = 1;
c = 0;

[x,y] = meshgrid(-w:w);
x0 = 0;
y0 = 0;
c = 0;
g = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));
dA = g;
dc = ones(size(g));
dx0 = (x-x0)/sigma^2*A.*g; 
dy0 = (y-y0)/sigma^2*A.*g;

F = [sum(dx0(:).^2) 0 0 0;
     0 sum(dy0(:).^2) 0 0;
     0 0 sum(dA(:).^2) sum(dA(:).*dc(:));
     0 0 sum(dA(:).*dc(:)) sum(dc(:).^2)];

 
crbV = sqrt(diag(inv(F)));


%%
g = simGaussianSpots(21, 21, sigma, 'x', 11, 'y', 11);

A = 1;
N = 1e3;

psnrV = 1:30;
ns = numel(psnrV);
A_crb = zeros(ns,1);
A_std = zeros(ns,1);
A_stdP = zeros(ns,1);
for k = 1:ns

    sigma_n = sqrt(A^2 * ni / (ni-1) / 10^(psnrV(k)/10));

    A_est = NaN(N,1);
    A_pstd = NaN(N,1);
    x_pstd = NaN(N,1);
    for i = 1:N
        gn = g + sigma_n*randn(21,21);
        %pstruct = fitGaussians2D(gn, 11, 11, 1, sigma, 0);
        pstruct = pointSourceDetection(gn, sigma);
        if ~isempty(pstruct)
            A_est(i) = pstruct.A;
            A_pstd(i) = pstruct.A_pstd;
            x_pstd(i) = pstruct.x_pstd;
        end
    end
    A_crb(k) = sigma_n*crbV(3);
    A_std(k) = nanstd(A_est);
    A_stdP(k) = nanmean(A_pstd);
end

figure;
hold on;
plot(psnrV, A_crb, 'k-');
plot(psnrV, A_std, 'r--');
plot(psnrV, A_stdP, 'g--');
xlabel('PSNR [dB]');
legend('CRB', 'Measured s.d.', 'Propagated s.d.');
