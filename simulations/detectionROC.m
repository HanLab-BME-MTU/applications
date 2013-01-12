
% Image parameters

nx = 400;
ny = 400;
sigma = 1.4;
w = ceil(4*sigma);
ni = 2*w+1; % width of PSF support

np = 500; % # points

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
img = simGaussianSpots(nx, ny, sigma, 'x', x, 'y', y);

figure;
imagesc(img); colormap(gray(256)); axis image;
hold on;
% plot(x, y, 'ro');
%%

% loop through PSNR values
psnrV = 0:1:25;
ns = numel(psnrV);

alphaV = [0.05 0.10 0.2 0.5 1];
na = numel(alphaV);

% nTP = NaN(ns,na);
% nFP = NaN(ns,na);
% nTN = NaN(ns,na);
% nFN = NaN(ns,na);
tpr = NaN(ns,na);
fpr = NaN(ns,na);

N = 200;

for ai = 1:na
    for k = 1:ns
        sigma_n = sqrt(A^2 * ni / (ni-1) / 10^(psnrV(k)/10));

        itpr = NaN(N,1);
        ifpr = NaN(N,1);
        parfor i = 1:N
        
            img_n = img + sigma_n*randn(ny,nx);
            
            [pstruct, ~, imgLM] = pointSourceDetection(img_n, sigma, 'Alpha', alphaV(ai));
            if ~isempty(pstruct)
                [ly,lx] = find(imgLM~=0);
                % remove loc. max. within image border
                rmIdx = lx<=w | ly<=w | lx>nx-w | ly>ny-w;
                lx(rmIdx) = [];
                ly(rmIdx) = [];
                
%                 figure; imagesc(img_n); colormap(gray(256)); colorbar; axis image;
%                 hold on; 
%                 plot(x, y, 'go');
%                 %plot(lx, ly, 'ws'); 
%                 plot(pstruct.x, pstruct.y, 'rx');
                
                
                % match local maxima to closest true position
                D = createSparseDistanceMatrix([lx ly], [x y], 10);
                [linkMax2In, ~] = lap(D, [], [], 1);
                linkMax2In = double(linkMax2In(1:numel(lx)));
                linkMax2In(linkMax2In>numel(x)) = NaN;
                
                % classify loc. max. wrt input
                tfidx = ~isnan(linkMax2In);
                
                % match detections to closest local maximum
                D = createSparseDistanceMatrix([pstruct.x' pstruct.y'], [lx ly], 10);
                [linkDet2Max, ~] = lap(D, [], [], 1);
                linkDet2Max = double(linkDet2Max(1:numel(pstruct.x)));
                didx = tfidx(linkDet2Max);
                didx2 = false(size(tfidx));
                didx2(linkDet2Max) = true;
                
                
                % display match
                %         lxi = [lx(tfidx) x(linkMax2In(tfidx))];
                %         lyi = [ly(tfidx) y(linkMax2In(tfidx))];
                %         dxi = [lx(linkDet2Max(didx)) pstruct.x(didx)'];
                %         dyi = [ly(linkDet2Max(didx)) pstruct.y(didx)'];
                %
                %         figure;
                %         imagesc(img_n); colormap(gray(256)); axis image;
                %         hold on;
                %         plot(lxi', lyi', 'r');
                %         plot(dxi', dyi', 'm');
                %         plot(lx, ly, 'wo');
                %         plot(x, y, 'go');
                %         plot(pstruct.x, pstruct.y, 'ys');
                
                
                tp = didx2 & tfidx;
                fp = didx2 & ~tfidx;
                tn = ~didx2 & ~tfidx;
                fn = ~didx2 & tfidx;
                ntp = sum(tp);
                nfp = sum(fp);
                ntn = sum(tn);
                nfn = sum(fn);
                
                itpr(i) = ntp/(ntp+nfn);
                ifpr(i) = nfp/(nfp+ntn);
            end
        end
        tpr(k,ai) = nanmedian(itpr);
        fpr(k,ai) = nanmedian(ifpr);
            %nTP(k,ai) = sum(tp);
            %nFP(k,ai) = sum(fp);
            %nTN(k,ai) = sum(tn);
            %FN(k,ai) = sum(fn);
    end
end

save roc
%%
% figure;
% hold on;
% plot(psnrV, nTP, 'g');
% plot(psnrV, nFP, 'r');
% plot(psnrV, nTN, 'k');
% plot(psnrV, nFN, 'b');
% legend('TP', 'FP', 'TN', 'FN');

%%
cmap = jet(ns);
fset = loadFigureSettings();

cv = jet(na);

figure(fset.fOpts{:});
axes(fset.axOpts{:});
colormap(cmap);
hold on;
for ai = 1:na

    %tpr = nTP(:,ai) ./ (nTP(:,ai)+nFN(:,ai));
    %fpr = nFP(:,ai) ./ (nFP(:,ai)+nTN(:,ai));

    %mesh([fpr fpr], [tpr tpr], zeros(ns,2), repmat((1:ns)', [1 2]),...
    %    'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.5);
    
    plot(fpr(:,ai), tpr(:,ai), 'Color', cv(ai,:), 'LineWidth', 2);
    text(fpr(1,ai)+0.02, tpr(1,ai), num2str(alphaV(ai)))
    
%     set(gca, 'ColorOrder', cmap);
%     plot([fpr fpr]', [tpr tpr]', 'o');
end
axis([0 1 0 1]);
axis square;
xlabel('False positive rate', fset.lfont{:});
ylabel('True positive rate', fset.lfont{:});

% print('-depsc2', 'detectionROC.eps');

% figure; hold on; plot(psnrV, tpr, 'g'); plot(psnrV, fpr, 'r');

%%
% figure;
% imagesc(img_n); colormap(gray(256)); axis image;
% hold on;
% plot(x, y, 'go');
% plot(lx(tp), ly(tp), 'yx');
% plot(lx(fp), ly(fp), 'rx');
% plot(lx(tn), ly(tn), 'wx');
% plot(lx(fn), ly(fn), 'cx');
% 
% plot(pstruct.x(unmatchedIdx), pstruct.y(unmatchedIdx), 'ms');

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
%========================================================
% CRB comparison: A vs. sigma_x
%========================================================
x0 = 0;
y0 = 0;
w = 6;
sigma = 1.4;
c = 0;

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

figure;
hold on;
plot(av, crbG, 'k');
plot(av, xstdG, 'k.');
plot(av, meanxpstdG, 'ko');
plot(av, crbP, 'b');
plot(av, xstdP, 'b.');
plot(av, meanxpstdP, 'bo');
axis([0 70 0 0.8])
legend('CRB(x) Gaussian noise', '\sigma_x', '\sigma_x (propagated)',...
    'CRB(x) Poisson noise', '\sigma_x', '\sigma_x (propagated)');
xlabel('A');
ylabel('\sigma_x');

figure;
hold on;
plot([0 30], [0 30], 'k--');
plot(psnrG, psnrEstG, 'r.-');
plot(psnrP, psnrEstP, 'g.-');
axis([0 30 0 30]); axis square;


%%
%========================================================
% PSNR estimation with Gaussian and Poisson noise
%========================================================
x0 = 0;
y0 = 0;
sigma = 1.4;
c = 10;

w = ceil(5*sigma);
[x,y] = meshgrid(-w:w);

g = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));
ni = (2*ceil(4*sigma)+1)^2;

sigma_n = 3;

psnr0 = 1;
psnrN = 30;
a0 = sqrt(sigma_n^2 * 10^(psnr0/10) * (ni-1)/ni);
aN = sqrt(sigma_n^2 * 10^(psnrN/10) * (ni-1)/ni);
np = 30;
av = linspace(a0, aN, np);

% reference PSNR for Gaussian noise
psnrG = 10*log10(av.^2/sigma_n^2*ni/(ni-1));

% expected PSNR for Poisson noise
psnrP = 10*log10(av.^2 ./ (av*mean(g(:))+c));

N = 1e2;

psnrEstG = NaN(np,1);
psnrEstP = NaN(np,1);
for k = 1:np
    % Gaussian
    iPSNR = NaN(N,1);
    parfor i = 1:N
        gn = av(k)*g + c + sigma_n*randn(2*w+1,2*w+1);
        pstruct = fitGaussians2D(gn, w+1, w+1, av(k), sigma, c);
        iPSNR(i) = 10*log10(pstruct.A^2*ni/pstruct.RSS);
    end
    psnrEstG(k) = nanmean(iPSNR);
    
    % Poisson
    iPSNR = NaN(N,1);
    parfor i = 1:N
        gn = poissrnd(av(k)*g + c);
        pstruct = fitGaussians2D(gn, w+1, w+1, av(k), sigma, c);
        iPSNR(i) = 10*log10(pstruct.A^2/(mean(gn(:))));
    end
    psnrEstP(k) = nanmean(iPSNR);
end

figure;
hold on;
plot([0 30], [0 30], 'k--', 'HandleVisibility', 'off');
plot(psnrG, psnrEstG, 'r.-');
plot(psnrP, psnrEstP, 'g.-');
axis([0 30 0 30]);
axis square;
xlabel('PSNR')
ylabel('Estimated PSNR');
legend('Gaussian noise', 'Poisson noise', 'Location', 'SouthEast');

