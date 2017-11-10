%========================================================
% CRB comparison: A vs. sigma_x
%========================================================
% Validation of localization precision for fitGaussian2D
% Gaussian and Poisson noise cases are compared

x0 = 0;
y0 = 0;
w = 6;
sigma = 1.4;
c = 10;

ni = (2*ceil(4*sigma)+1)^2; % support of the fit

[x,y] = meshgrid(-w:w);
g = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));


% Gaussian noise std
sigma_n = 3;

% lowest amplitude: PSNR = 1
psnr0 = 1;
psnrN = 30;
a0 = sigma_n*sqrt(10^(psnr0/10)*(ni-1)/ni);
aN = sigma_n*sqrt(10^(psnrN/10)*(ni-1)/ni);
na = 50;
% amplitude vector
av = linspace(a0, aN, na);
psnrG = 10*log10(av.^2/sigma_n^2*ni/(ni-1));

% expected PSNR for Poisson noise
psnrP = 10*log10(av.^2 ./ (av*mean(g(:))+c));

% ai = 40;
% gnG = av(ai)*g+c + sigma_n*randn(size(g));
% gnP = poissrnd((av(ai)*g+c));
% figure;
% colormap(gray(256));
% subplot(1,3,1); imagesc(av(ai)*g+c); axis image; colorbar;
% subplot(1,3,2); imagesc(gnG); axis image; colorbar;
% subplot(1,3,3); imagesc(gnP); axis image; colorbar;


crbG = sigma_n*sqrt(2/pi)./av;
% crbP = 1./sqrt(av*2*pi); % when c=0
crbP = sqrt(arrayfun(@(i) 1/sum(sum( ((x-x0)./sigma^2*i.*g).^2 ./ (i*g+c))), av) );

meanxpstdG = zeros(na,1);
meanxpstdP = zeros(na,1);
psnrEstG = zeros(na,1);
psnrEstP = zeros(na,1);
xestG = zeros(na,1);
xestP = zeros(na,1);
N = 1e3;

for ai = 1:na
    % Gaussian
    ixestG = NaN(N,1);
    ixpstdG = NaN(N,1);
    ipsnr = NaN(N,1);
    parfor i = 1:N
        gn = av(ai)*g+c + sigma_n*randn(size(g));
        pstruct = fitGaussians2D(gn, w+1, w+1, av(ai), sigma, c);
        ixestG(i) = pstruct.x;
        ixpstdG(i) = pstruct.x_pstd;
        ipsnr(i) = 10*log10(pstruct.A^2*ni/pstruct.RSS);
    end
    xestG(ai) = nanmean(ixestG);
    xstdG(ai) = nanstd(ixestG);
    meanxpstdG(ai) = nanmean(ixpstdG);
    psnrEstG(ai) = nanmean(ipsnr);
    
    % Poisson
    ixestP = NaN(N,1);
    ixpstdP = NaN(N,1);
    ipsnr = NaN(N,1);
    parfor i = 1:N
        gn = poissrnd(av(ai)*g+c);
        pstruct = fitGaussians2D(gn, w+1, w+1, av(ai), sigma, c);
        ixestP(i) = pstruct.x;
        ixpstdP(i) = pstruct.x_pstd;
        ipsnr(i) = 10*log10(pstruct.A^2/mean(gn(:)));
    end
    xestP(ai) = nanmean(ixestP);
    xstdP(ai) = nanstd(ixestP);
    meanxpstdP(ai) = nanmean(ixpstdP);
    psnrEstP(ai) = nanmean(ipsnr);
end

fset = loadFigureSettings();
figure(fset.fOpts{:});
axes(fset.axOpts{:});
hold on;

plot(av, crbG, 'k', 'LineWidth', 1, 'MarkerSize', 16);
plot(av, xstdG, 'k.', 'LineWidth', 1, 'MarkerSize', 16);
plot(av, meanxpstdG, 'ko', 'LineWidth', 1, 'MarkerSize', 8);
plot(av, crbP, 'b', 'LineWidth', 1, 'MarkerSize', 10);
plot(av, xstdP, 'b.', 'LineWidth', 1, 'MarkerSize', 16);
plot(av, meanxpstdP, 'bo', 'LineWidth', 1, 'MarkerSize', 8);
axis([0 70 0 0.8])
legend('CRB(x) Gaussian noise', '\sigma_x', '\sigma_x (propagated)',...
    'CRB(x) Poisson noise', '\sigma_x', '\sigma_x (propagated)');
xlabel('A');
ylabel('\sigma_x (pixels)');
print('-depsc2', ['GaussVsPoissonCRB_c=' num2str(c) '.eps']);

figure(fset.fOpts{:});
axes(fset.axOpts{:});
hold on;
plot([0 30], [0 30], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(psnrG, psnrEstG, 'r.-', 'LineWidth', 1, 'MarkerSize', 16);
plot(psnrP, psnrEstP, 'g.-', 'LineWidth', 1, 'MarkerSize', 16);
axis([0 30 0 30]); axis square;

xlabel('Reference PSNR (dB)');
ylabel('Calculated PSNR (dB)');
hl = legend('Gaussian noise', 'Poisson noise', 'Location', 'NorthWest');
set(hl, 'Box', 'off');
print('-depsc2', ['PSNRest_c=' num2str(c) '.eps']);
