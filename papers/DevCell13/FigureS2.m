%=========================================================================================
% Figure S2
%=========================================================================================
% This script generates the panels for Figure S2 from Aguet et al., Dev. Cell, 2013.

spath = '/Users/aguet/Documents/Papers/DevCell13/FigureS2/';
fset = loadFigureSettings('print');
printFigures = false;
loadData = false;

%%
%====================================================
% Panel A: 2-D illustration of the detection approach
%====================================================

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [180 180];

np = 30;
if loadData
     load([spath 'FigS2a.mat']); %#ok<*UNRCH>
else
    sigma = 1.5;
    sigma_n = 1.25;
    [frame, xv, yv, sv, Av] = simGaussianSpots(160, 160, sigma, 'npoints', np, 'A', 4+0.5*randn(np,1), 'Border', 'padded');
    frame = frame + 10 + sigma_n*randn(size(frame));
end

dRange = [min(frame(:)) max(frame(:))];

figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frame); colormap(gray(256)); axis image; colormap;
axis off;
if printFigures
    print('-depsc', '-loose', [spath 'detectionSteps_1_raw.eps']);
end


% pointSourceDetection
img = frame;
w = ceil(4*sigma);
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
u = ones(1,length(x));

% convolutions
imgXT = padarrayXT(img, [w w], 'symmetric');
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');
fu2 = conv2(u', u, imgXT.^2, 'valid');

% Laplacian of Gaussian
gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

% 2-D kernel
g = g'*g;
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% solution to linear system
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

J = [g(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
sigma_e2 = RSS/(n-3);

sigma_A = sqrt(sigma_e2*C(1,1));

% standard deviation of residuals
sigma_res = sqrt((RSS - (A_est*gsum+n*c_est - fu)/n)/(n-1));

kLevel = norminv(1-0.05/2, 0, 1);

SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
T = (A_est - sigma_res*kLevel) ./ scomb;
pval = tcdf(-T, df2);

% mask of admissible positions for local maxima
mask = pval < 0.05;

mask0RGB = rgbOverlay(frame, 255*mask, [1 0 0], dRange);

% all local max
allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);

% local maxima above threshold in image domain

imgLM = allMax .* mask;
[lmAy, lmAx] = find(allMax~=0);
[lm0y, lm0x] = find(imgLM~=0);

% -> set threshold in LoG domain
logThreshold = min(imgLoG(imgLM~=0));
logMask = imgLoG >= logThreshold;

% combine masks
mask = mask | logMask;

% re-select local maxima
imgLM = allMax .* mask;

[lmy, lmx] = find(imgLM~=0);


% figure('Position', pos, 'PaperPositionMode', 'auto');
% axes('Position', [0 0 1 1]);
cla;
imagesc(imgLoG); colormap(gray(256)); axis image; colormap;
axis off;
if printFigures
    print('-depsc', '-loose', [spath 'detectionSteps_2a_LoG.eps']);
end
hold on;
hp = plot(lm0x, lm0y, 'rx');
if printFigures
print('-depsc', '-loose', [spath 'detectionSteps_3_LoG+locmax.eps']);
end
delete(hp);
hp = plot(lmAx, lmAy, 'rx');
% print('-depsc', '-loose', [rpath 'detectionSteps_2_LoG+locmax_all.eps']);

cla
imagesc(mask0RGB); axis off;
if printFigures
    print('-depsc', '-loose', [spath 'detectionSteps_2b_mask0.eps']);
end

cla
imagesc(img); 
hold on;
hp = plot(lmx, lmy, 'rx');
axis off;
% print('-depsc', '-loose', [rpath 'detectionSteps_4_locmaxFinal.eps']);
% delete(hp);
[pstruct, maskPSD, imgLM, imgLoG] = pointSourceDetection(frame, sigma);
plot(pstruct.x, pstruct.y, 'go', 'MarkerSize', 8, 'LineWidth', 1);
% plot(xv, yv, 'go', 'MarkerSize', 8, 'LineWidth', 1);
if printFigures
    print('-depsc', '-loose', [spath 'detectionSteps_4_Localization.eps']);
end

cla;
imagesc(img); 
plot(xv, yv, 'yo', 'MarkerSize', 8, 'LineWidth', 1);
if printFigures
    print('-depsc', '-loose', [spath 'detectionSteps_5_GroundTruth.eps']);
end
% plot(lmx, lmy, 'go');

%%
%----------------------------------------------------------------- 
% Compare vs: t-test of A±A_pstd vs. 0±sigma_r (background)
%-----------------------------------------------------------------

df2 = (n-1) * (sigma_A.^2 + sigma_res.^2).^2 ./ (sigma_A.^4 + sigma_res.^4);
scomb = sqrt((sigma_A.^2 + sigma_res.^2)/n);
T = A_est ./ scomb;
pval = tcdf(-T, df2);
mask = pval < 0.05;
maskRGB = rgbOverlay(frame, 255*mask, [1 0 0], dRange);
cla
imagesc(maskRGB);
% print('-depsc', '-loose', [rpath 'detectionSteps_Mask_Avs0.eps']);



%%
%====================================================
% Panel B: 1-D illustration of the detection approach
%====================================================
alpha = 0.05;
kLevel = norminv(1-alpha, 0, 1); % ~2 std above background


c = 3;
sBgr = 1;
A = 2.5;
sigma = 3;
ymax = 7;

w = ceil(4*sigma);
x = -50:50;
nx = numel(x);
% x0 = x(end)/2;
x0 = 0;
x_fine = x(1):0.1:x(end);

g = A*exp(-(x-x0).^2/(2*sigma^2));
g_fine = A*exp(-(x_fine-x0).^2/(2*sigma^2));

noise = sBgr*randn(1,nx);
signal = c + g + noise;
% save FigS2bNEW x signal noise
% tmp = load('FigS2b.mat');
%============================
load('/Users/aguet/Documents/MATLAB/endocytosis/CMEpaper/Figure S2 - Detection/FigS2bNEW.mat');
%============================



[prmVect, prmStd, C, res, J] = fitGaussian1D(x, signal, [0 max(signal) sigma min(signal)], 'xAc');
x_est = prmVect(1);
x_pstd = prmStd(1);
A_est = prmVect(2);
A_pstd = prmStd(2);
c_est = prmVect(4);

sBgr_est = std(res);
SE = sBgr_est/sqrt(2*(nx-1));
SE_sigma_r = SE*kLevel;
npx = nx;
sigma_A = A_pstd;

% perform test:
% df2 = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
% scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
% T = (A_est - res.std*kLevel) ./ scomb;
% % 1-sided t-test: A_est must be greater than k*sigma_r
% pStruct.pval_Ar(p) = tcdf(-T, df2);
df2 = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
T = (A_est - sBgr_est*kLevel) ./ scomb;
pval = tcdf(-T, df2);
hval = pval<0.05

% gaussian convolved w/ std. dev on x
a = 2*x_pstd;
gx = (erf((a-2*(x_fine-x_est))/(2*sqrt(2)*sigma)) + erf((a+2*(x_fine-x_est))/(2*sqrt(2)*sigma))) /...
    (2*erf(a/(2*sqrt(2)*sigma)));

g_est = exp(-(x_fine-x_est).^2/(2*sigma^2));
%g_plus = (A_est+A_pstd)*g_est + c_est;
%g_min = (A_est-A_pstd)*g_est + c_est;
g_plus = (A_est+A_pstd)*gx + c_est;
g_min = (A_est-A_pstd)*gx + c_est;
g_est = A_est*g_est + c_est;


% % plot exact signal first, then noise
% pos = get(0, 'DefaultFigurePosition');
% pos(3:4) = [600 300];
% figure(fset.fOpts{:});
% axes(fset.axOpts{:});
% plot(x_fine, g_fine+c, 'k', 'LineWidth', 1);
% axis([x(1) x(end)+5 0 ymax]);
% set(gca, fset.axOpts{:}, 'Layer', 'bottom', 'XTickLabel', [], 'YTickLabel', [], 'box', 'off');
% ylabel('Intensity (A.U.)', fset.lfont{:});
% xlabel('x (A.U.)', fset.lfont{:});
% % print('-depsc2', [spath 'detection1D_model.eps']);
% 
% figure(fset.fOpts{:});
% axes(fset.axOpts{:});
% plot(x, signal, 'k', 'LineWidth', 1);
% axis([x(1) x(end)+5 0 ymax]);
% set(gca, fset.axOpts{:}, 'Layer', 'bottom', 'XTickLabel', [], 'YTickLabel', [], 'box', 'off');
% ylabel('Intensity (A.U.)', fset.lfont{:});
% xlabel('x (A.U.)', fset.lfont{:});
% % print('-depsc2', [spath 'detection1D_signal_H' num2str(hval) '.eps']);

%%
%------------------------------------
% Plot residual/noise
%------------------------------------

figure(fset.fOpts{:}, 'Position', [2 2 11 5.5000]);
axes(fset.axOpts{:})
hold on;
% setupFigure();
plot(x, noise+c, 'k', 'LineWidth', 1);
% g00 = A_est*exp(-(x-x_est).^2/(2*sigma^2));
% plot(x, signal-g00, 'r--', 'LineWidth', 1);
% plot(x, signal, 'b', 'LineWidth', 1);


plot([-w -w w w], 6.5+[0 0.5 0.5 0], 'k', 'LineWidth', 1);
axis([x(1) x(end)+5 0 ymax]);
text(0, 7.25, 'Support of fit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

set(gca, fset.axOpts{:}, 'Layer', 'bottom', 'XTick', -50:20:50, 'XTickLabel', [], 'YTickLabel', []);
ylabel('Intensity (A.U.)', fset.lfont{:});
hx = xlabel('x (A.U.)', fset.lfont{:});
posA = get(hx, 'Position');

print('-depsc2', '-loose', [spath 'residual' num2str(hval) '_alpha=' num2str(alpha) '.eps']);


%%

%------------------------------------
figure(fset.fOpts{:}, 'Position', [2 2 11 5.5000]);
axes(fset.axOpts{:})
hold on;
plot(x, signal, 'k', 'LineWidth', 1);
fill([x_fine x_fine(end:-1:1)], [g_plus g_min(end:-1:1)], fset.cfB, 'EdgeColor', 'none');
plot(x_fine, g_plus, '-', 'Color', fset.ceB, 'LineWidth', 0.5);
plot(x_fine, g_min, '-', 'Color', fset.ceB, 'LineWidth', 0.5);
plot(x_fine, g_est, '-', 'Color', fset.ceB, 'LineWidth', 1);
plot([x_est x_fine(end)+30], (A_est+c_est)*[1 1], '--', 'Color', fset.ceB, 'LineWidth', 1);
plot([x(1) 130], c_est+kLevel*sBgr_est*[1 1], 'r--', 'LineWidth', 1);
plot([-w -w w w], 6.5+[0 0.5 0.5 0], 'k', 'LineWidth', 1);
axis([x(1) x(end)+5 0 ymax]);
text(0, 7.25, 'Support of fit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

set(gca, fset.axOpts{:}, 'Layer', 'bottom', 'XTick', -50:20:50, 'XTickLabel', [], 'YTickLabel', []);
ylabel('Intensity (A.U.)', fset.lfont{:});
hx = xlabel('x (A.U.)', fset.lfont{:});
posA = get(hx, 'Position');



%------------------------------------
% PDF projetions
%------------------------------------
xx = 0:0.01:ymax;
axes(fset.axOpts{:}, 'Position', [7.5 1.5 3 3.5]);
hold on;

xf = linspace(c_est+kLevel*sBgr_est, ymax, 100);
fill([xf xf(end:-1:1)], [exp(-(xf-c_est).^2/(2*sBgr_est^2)) zeros(1,100)], 0.6*[1 1 1], 'EdgeColor', 'none')
% exp(-(lin-c).^2/(2*sBgr^2))

plot(xx, exp(-(xx-c_est).^2/(2*sBgr_est^2)), 'k', 'LineWidth', 1); %/sqrt(2*pi)/sBgr


% h = exp(-(kLevel*sBgr).^2/(2*sBgr^2)) / sqrt(2*pi)/sBgr;
plot((c_est+kLevel*sBgr_est)*[1 1], [0 1], 'r--', 'LineWidth', 1)


plot((A_est+c_est)*[1 1], [0 1], '--', 'Color', fset.ceB, 'LineWidth', 1);
plot(xx, exp(-(xx-(A_est+c_est)).^2/(2*A_pstd^2)), 'Color', fset.ceB, 'LineWidth', 1);%/sqrt(2*pi)/A_pstd

plot(xx, exp(-(xx-(c_est+kLevel*sBgr_est)).^2/(2*SE^2)), 'r', 'LineWidth', 1);
axis([0 ymax -0.02 1.02]);

view(90,-90)
set(gca, 'XTick', [], 'XColor', [1 1 1], 'Layer', 'bottom', 'TickLength', [0 0], 'YTick', []);

ylabel('Probability', fset.lfont{:})
print('-depsc2', '-loose', [spath 'detection1D_H' num2str(hval) '_alpha=' num2str(alpha) '.eps']);


%%
%====================================================
% Panel C: ROC curve
%====================================================
% detectionROC();
% load roc800_2


%%
%====================================================
% Panel D: SNR series
%====================================================
sigma = 1.4;

%----------------------------------------------------
% u-track parameters:
clear movieParam detectionParam;
movieParam.imageDir = [pwd filesep]; %directory where images are
movieParam.filenameBase = 'snrTest'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 1; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

detectionParam.psfSigma = sigma;
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.alphaLocMax = 0.1;%0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration
%----------------------------------------------------

psnrV = [-5:15 20:5:30];
% psnrV = [0:15 20:5:30];

% # noise iterations
N = 1e4;

% nx = 31;
% ny = 31;
nx = 101; % u-track has many FN with smaller templates (background estimation problem?)
ny = 101;

% PSF parameters
x0 = (nx+1)/2;
y0 = (ny+1)/2;
A = 1;
g0 = simGaussianSpots(nx, ny, sigma, 'x', x0, 'y', y0, 'A', A);
c = 10;

ni = (2*ceil(4*sigma)+1)^2; % # points used for fit

n = numel(g0);
ns = numel(psnrV);

A = 10;
sigmaV = sqrt(A^2 * n / (n-1) ./ 10.^(psnrV/10));
av = A*ones(size(sigmaV));
% sigma_n = 1;
% av = sqrt((n-1)/n*10.^(psnrV/10) * sigma_n^2);
% sigmaV = sigma_n*ones(size(av));
%%
tpr = zeros(ns,1);
tprUT = zeros(ns,1);
parfor i = 1:ns
    %sigma_n = sqrt(A^2 * n / (n-1) / 10^(psnrV(i)/10));
    psnrEst = NaN(1,N);
    imovieParam = movieParam;
    imovieParam.filenameBase = ['snr' num2str(psnrV(i)) 'Test'];

    iDetected = zeros(N,1);
    iDetectedUT = zeros(N,1);
    
    % loop through noise iterations
    for k = 1:N
        
        noise = sigmaV(i)*randn(ny,nx);        
        gn = av(i)*g0 + c + noise;
        
        % model-based detection
        [pstruct] = pointSourceDetection(gn, sigma, 'Alpha', 0.1);%, 'RefineMaskLoG', false);
        
        % true positive counted if within 2*sigma radius
        if ~isempty(pstruct)
            d = sqrt((pstruct.x-x0).^2 + (pstruct.y-y0).^2);
            idx = find(d==min(d), 1, 'first');
            if d(idx)<=2*sigma
                psnrEst(k) = 10*log10(pstruct.A(idx)^2*ni/pstruct.RSS(idx));
                iDetected(k) = 1;
            end
        end

        % u-track detection        
        imwrite(uint16(scaleContrast(gn,[],[0 65535])), [imovieParam.filenameBase '001.tif']);
        movieInfo = detectSubResFeatures2D_StandAlone(imovieParam,detectionParam,0,0);
        if ~isempty(movieInfo.xCoord)
            d = sqrt((movieInfo.xCoord(:,1)-x0).^2 + (movieInfo.yCoord(:,1)-y0).^2);
            idx = find(d==min(d), 1, 'first');
            if d(idx)<=2*sigma
                iDetectedUT(k) = 1;
            end
        end
        
        % wavelet-based detection (calculate ROC -> high sensitiviy but low selectivity)
        %[frameInfoW, imgDenoised] = main283AUTO_standalone(gn,1);
        %[frameInfoWavelets, maskWavelets] = detectSpotsWT(gn, 4, 5, 0);
        
        % many false positives in background: search center only
        %candIdx = x0-w<=frameInfoW.xav & frameInfoW.xav<=x0+w & y0-w<=frameInfoW.yav & frameInfoW.yav<=y0+w;
        %iNumDetectWT(k) = numel(frameInfoW.xav); % # detected objects
        %iNumTruePosWT(k) = sum(candIdx)==1; % true positives
        %iNumFalsePosWT(k) = sum(~candIdx);        
    end
    tpr(i) = sum(iDetected)/N;
    tprUT(i) = sum(iDetectedUT)/N;
    psnrEst(isnan(psnrEst)) = [];
    fprintf('Inputs PSNR: %.2f, est: %.2f\n', psnrV(i), median(psnrEst));
end
% save Fig2d.mat psnrV tpr tprUT;
%%
% load([spath 'Fig2d.mat']);
load([spath 'Fig2dx101.mat']);

snrV = 10.^(psnrV/10);

w = ceil(6*sigma);
[x,y] = meshgrid(-w:w);
[ny,nx] = size(x);
% PSF model
g0 = exp(-(x.^2+y.^2)/(2*sigma^2));

c = 0;
A = 1;

g = A*g0 + c;
n = numel(g);

% g_ex = cell(1,ns);
% for i = 1:ns
%     sigma_n = sqrt(A^2 * n / (n-1) / 10^(psnrV(i)/10));
%     g_ex{i} = g + sigma_n*randn(ny,nx);
% end
% gcat = cat(3, g_ex{:});
% dRange = [min(gcat(:)) max(gcat(:))];

% SNR inset
snrInset = [1 3 6 10 30 100];
g_ex = cell(1,numel(snrInset));
for i = 1:numel(snrInset)
    sigma_n = sqrt(A^2 * n / (n-1) / snrInset(i));
    g_ex{i} = g + sigma_n*randn(ny,nx);
end
gcat = cat(3, g_ex{:});



fset = loadFigureSettings('print');
figure(fset.fOpts{:}, 'Position', [2 2 8 6.5]);

% grid below data
axes(fset.axOpts{:});
set(gca, 'xscale', 'log');
grid on;
axis([1 100 0 1.02]);


set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', '-', 'XTickLabel', [], 'YTickLabel', [],...
    'XColor', 0.75*[1 1 1], 'YColor', 0.75*[1 1 1]);

% plot SNR curves
axes(fset.axOpts{:});
hold on;
set(gca, 'xscale', 'log');
axis([1 100 0 1.02]);


plot(snrV, tpr, 'Color', hsv2rgb([0.55 1 1]), 'LineWidth', 1);
plot(snrV, tprUT, 'Color', hsv2rgb([0.05 1 1]), 'LineWidth', 1);

set(gca, fset.axOpts{:}, 'Color', 'none', 'XTickLabel', [1 10 100]);
xlabel('SNR', fset.lfont{:});
ylabel('Detection sensitivity', fset.lfont{:});
hl = legend(' This work', ' u-track', 'Location', 'NorthEast');
set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [5.25 2.15 1.5 0.8]);


% panels with different SNRs
dx = 0.06;
w = 0.95;
% idx = find(ismember(psnrV, [1 5 8 10 15 20]));
for i = 1:numel(snrInset)
    axes(fset.axOpts{:}, 'Position', [1.5+(i-1)*(dx+w) 5.1 w w]);
    imagesc(g_ex{i}); colormap(gray(256)); axis image off;
    text(nx,0, num2str(snrInset(i)), 'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    if i==1
       text(-1,0, 'SNR', 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    end
end
% print('-depsc', '-loose', [spath 'detection_SNR.eps']);
print('-depsc2', '-loose', 'test.eps');



%%
% verify t-test
disp('----------');
n1 = 50;
sample1 = 3+0.5*randn(1,n1);
n2 = n1;
sample2 = 2+0.5*randn(1,n2);
s1 = std(sample1);
s2 = std(sample2);

df2 = (n1-1) * (s1.^2 + s2.^2).^2 ./ (s1.^4 + s2.^4);
scomb = sqrt((s1.^2 + s2.^2)/n1);

T = (mean(sample1) - mean(sample2)) ./ scomb
pval = tcdf(-T, df2)

[hval, pval, ci, tstats] = ttest2(sample1, sample2, alpha, 'right', 'unequal')



%%
%====================================================
% Comparison with tracks from wavelet detection
%====================================================
runWaveletDetection(dataOX, 'Overwrite', true);
runTracking(dataOX, loadTrackSettings('Radius', [4 7]), 'Overwrite', true, 'DetectionFile', 'detection_v1.mat', 'FileName', 'trackedFeatures_v1.mat');
runWaveletTrackProcessing(dataOX, 'TrackerOutput', 'trackedFeatures_v1.mat', 'FileName', 'ProcessedTracks_v1.mat', 'Overwrite', true);
% Last run on 082712, 041113
%%
%----------------------------------------------------
% Mean lifetime histogram
%----------------------------------------------------
data = dataOX;
nd = numel(data);
cutoff_f = 5;
buffer = 5;
Nmax = max([data.movieLength])-2*buffer;
framerate = 1;

for i = 1:nd
    load([data(i).source 'Tracking/ProcessedTracks_v1.mat']);
    trackLengths = [tracks.end] - [tracks.start] + 1;
    lifetimes_s = [tracks.lifetime_s];
    lifetimes_s = lifetimes_s(trackLengths>=cutoff_f & [tracks.catIdx]==1);
    
    N = data(i).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    lftHist_Ia = hist(lifetimes_s, t);
    
    w = N./(N-cutoff_f+1:-1:1);
    pad0 = zeros(1,Nmax-N);
    lftHist_Ia =  [lftHist_Ia.*w  pad0];
    t = (cutoff_f:Nmax)*framerate;
    res(i).lftHist_Ia = lftHist_Ia / sum(lftHist_Ia) / framerate;
end
lftRes = runLifetimeAnalysis(dataOX);
%%
plotLifetimes(lftRes);
hold on;
meanHist = mean(vertcat(res.lftHist_Ia),1);

plot(t, meanHist, 'r');

figure;
plot(t, vertcat(res.lftHist_Ia));
axis([0 120 0 0.05]);

% mean lifetime, truncated & re-normalized
% lft = mean(lftRes.lftHist_Ia, 1);

%%
%----------------------------------------------------
% Comparison of tracks
%----------------------------------------------------
k = 3;
tracksW = load([dataOX(k).source 'Tracking/ProcessedTracks_v1.mat']);
tracksW = tracksW.tracks;
tracksW = tracksW([tracksW.catIdx]==1);

tracks = loadTracks(dataOX(k), 'Category', 'Ia', 'Cutoff_f', 5, 'Mask', true);

[map, masterIdx, slaveIdx, unassignedMasterIdx, unassignedSlaveIdx] = assignSlaveTracksToMaster(tracksW, tracks);

idx60 = find([tracks.lifetime_s]==60);
%%
figure;
hold on;
i = idx60(10);
plot(tracks(i).t, tracks(i).A + tracks(i).c);
plot(tracksW(map{1}{i}).t, tracksW(map{1}{i}).A/8, 'r');

%%
%==================================================
% New Panel D
%==================================================
% The following files were generated with FigureDetectionROC.m
fpath = '/Users/aguet/Documents/MATLAB/endocytosis/DetectionComparison/';
load([fpath 'LC_pfa=15.mat']);
load([fpath 'ImarisDetectionData.mat']);
load([fpath 'detectionROCData.mat']);

fset = loadFigureSettings('print');
figure(fset.fOpts{:}, 'Position', [2 2 8 9.5]);

% grid below data
axes(fset.axOpts{:}, 'Position', [1.5 5 6 3]);
set(gca, 'xscale', 'log');
grid on;
axis([1 100 0 1.02]);


set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', '-', 'XTickLabel', [], 'YTickLabel', [],...
    'XColor', 0.75*[1 1 1], 'YColor', 0.75*[1 1 1]);

% plot SNR curves
axes(fset.axOpts{:}, 'Position', [1.5 5 6 3]);
hold on;
set(gca, 'xscale', 'log', 'Color', 'none');
axis([1 100 0 1.02]);


cv = hsv2rgb([0.33 1 0.9
              0.55 1 1;
              0.55+.15 1 1;
              0.11 1 1;
              0.0 1 1]);         

% total points/SNR level: nrep
hp = zeros(5,1);
% plot(psnrV, nanmean(ntpWT2,2)/nrep, 'Color', hsv2rgb([0.85 1 1]), 'LineWidth', 1);
% plot(psnrV, nanmean(ntpTH,2)/nrep, 'Color', hsv2rgb([0 1 0]), 'LineWidth', 1);
hp(5) = plot(psnrV, nanmean(ntpIM,2)/nrep, 'Color', cv(5,:), 'LineWidth', 1);
hp(4) = plot(psnrV, nanmean(ntpWT,2)/nrep, 'Color', cv(4,:), 'LineWidth', 1);
hp(3) = plot(psnrV, nanmean(ntpUT,2)/nrep, 'Color', cv(3,:), 'LineWidth', 1);
hp(2) = plot(psnrV, nanmean(ntpLC,2)/nrep, 'Color', cv(2,:) , 'LineWidth', 1);
hp(1) = plot(psnrV, nanmean(ntpPS,2)/nrep, 'Color', cv(1,:), 'LineWidth', 1);
set(gca, 'XTickLabel', [1 10 100], 'XTickLabel', []);
hy(1) = ylabel('Detection sensitivity', fset.lfont{:});
% hl = legend(hp, ' PSF model', ' u-track', ' Wavelet', ' Imaris', 'Location', 'NorthEast');
% set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [8.15 1.45 1.25 1.5]);
formatTickLabels();


axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 3]);
set(gca, 'xscale', 'log');
grid on;
axis([1 100 0 0.25]);

set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', '-', 'XTickLabel', [], 'YTickLabel', [],...
    'XColor', 0.75*[1 1 1], 'YColor', 0.75*[1 1 1]);

% plot SNR curves
axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 3]);
hold on;
set(gca, 'xscale', 'log', 'Color', 'none');
axis([1 100 0 1]);

tp = nanmedian(ntpIMf,2);
fp = nanmedian(nfpIMf,2);
plot(psnrV, fp./(tp+fp), 'Color', cv(5,:), 'LineWidth', 1);
tp = nanmedian(ntpWTf,2);
fp = nanmedian(nfpWTf,2);
plot(psnrV, fp./(tp+fp), 'Color', cv(4,:), 'LineWidth', 1);
tp = nanmedian(ntpLCf,2);
fp = nanmedian(nfpLCf,2);
plot(psnrV, fp./(tp+fp), 'Color', cv(2,:), 'LineWidth', 1);
tp = nanmedian(ntpUTf,2);
fp = nanmedian(nfpUTf,2);
plot(psnrV, fp./(tp+fp), 'Color', cv(3,:), 'LineWidth', 1);


tp = nanmedian(ntpPSf,2);
fp = nanmedian(nfpPSf,2);
plot(psnrV, fp./(tp+fp), '--', 'Color', cv(1,:), 'LineWidth', 1);

% 1) PSNR scale, as in paper:
% set(gca, 'XTickLabel', [1 10 100], 'Layer', 'bottom');
% xlabel('PSNR', fset.lfont{:});

% 2) A/sigma
ni = (ceil(4*sigma)*2+1)^2;
psnrv = [1 2 3  5 7 10 20 30 100];
asvec = sqrt((ni)/ni*psnrv);
set(gca, 'XTick', psnrv, 'XTickLabel', arrayfun(@(i) num2str(i, '%.1f'), asvec, 'unif', 0), 'Layer', 'bottom');
psnrv = [4 6 9];
asvec = sqrt((ni)/ni*psnrv);
arrayfun(@(i) text(psnrv(i), -0.16, num2str(asvec(i), '%.1f'), 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'top', fset.sfont{:}), 1:numel(psnrv))
hx = xlabel('A/\sigma_n', fset.lfont{:});
posx = get(hx, 'Position');
posx(2) = 1.5*posx(2);
set(hx, 'Position', posx);


hy(2) = ylabel('False discovery rate', fset.lfont{:});
% hl = legend(hp, ' PSF model', ' u-track', ' Wavelet', ' Imaris', 'Location', 'NorthEast');
% set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [9 1.25 1.25 1.5]);
formatTickLabels('FormatXY', [false true]);

% align ylabels
ypos2 = get(hy(2), 'Position');
ypos1 = get(hy(1), 'Position');
set(hy(1), 'Position', [ypos2(1) ypos1(2:end)]);



wg = ceil(6*sigma);
[xg,yg] = meshgrid(-wg:wg);
% PSF model
g0 = exp(-(xg.^2+yg.^2)/(2*sigma^2));
n = numel(g0);

% panels with different SNRs
dx = 0.06;
pw = 0.95;

% snrInset = [1 3 6 10 30 100];
snrInset = [1 1.5 2 2.5 3 5].^2;
g_ex = cell(1,numel(snrInset));
for i = 1:numel(snrInset)
    g_ex{i} = g0 + sqrt(1^2 * n/(n-1) / snrInset(i)) * randn(size(g0));
end
gcat = cat(3, g_ex{:});


for i = 1:numel(snrInset)
    axes(fset.axOpts{:}, 'Position', [1.5+(i-1)*(dx+pw) 8.2 pw pw]);
    imagesc(g_ex{i}); colormap(gray(256)); axis image off;
    text(size(g0,1),0, num2str(sqrt(snrInset(i)), '%.1f'), 'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    %text(size(g0,1),0, num2str(snrInset(i)), 'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    if i==1
       %text(-1,0, 'PSNR', 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
       text(-1,0, 'A/\sigma_n', 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    end
end

% print('-depsc2', '-loose', 'detectionComparison.eps');
%%
%==================================================
% New Panel C
%==================================================
ifile = [fpath 'LC_pfa=15.mat'];
if exist(ifile, 'file')==2
    load(ifile);
else
    alphaV = [0.05 0.1 0.2 0.5];
    na = numel(alphaV);
    % noise repeats            [pstruct, ~, imgLM] = pointSourceDetection(img_n, sigma, 'Alpha', alphaV(ai)*2);
    
    N = 50; % 100, 250?
    
    ntpPS = cell(1,na);
    nfpPS = cell(1,na);
    ntpPSf = cell(1,na);
    nfpPSf = cell(1,na);
    parfor ai = 1:na
        
        ntpPS{ai} = zeros(ns,N);
        nfpPS{ai} = zeros(ns,N);
        
        psopts = {'Alpha', alphaV(ai), 'RedundancyRadius', 1, 'RefineMaskLoG', true};
        
        % loop through images with variable SNR
        for i = 1:N
            fprintf('%d\n', i);
            img_n = img + sigma_n*randn(ny,nx);
            
            [pstruct, ~, imgLM, imgLoG] = pointSourceDetection(img_n, sigma, psopts{:});
            D = createSparseDistanceMatrix([pstruct.x' pstruct.y'], [x y; xn yn], R); % search radius < min distance
            [linkDet2Ref, ~] = lap(D, [], [], 1);
            linkDet2Ref = double(linkDet2Ref(1:numel(pstruct.x)));
            vidx = linkDet2Ref(linkDet2Ref<=np);
            ntpPS{ai}(:,i) = hist(ceil(vidx/nrep), 1:ns);
            nfpPS{ai}(:,i) = hist(ceil(linkDet2Ref(linkDet2Ref>np)/nrep), 1:ns); % detections w/o corresponding ref: FP
        end
        
        % loop through images with same SNR
        N2 = 100;
        ntpPSf{ai} = zeros(ns,N2);
        nfpPSf{ai} = zeros(ns,N2);
        for k = 1:ns
            fprintf('%d\n', k);
            for i = 1:N2 % noise repeats
                img_n = av0(k)*img0 + sigma_n*randn(ny,nx);
                [pstruct, ~, imgLM, imgLoG] = pointSourceDetection(img_n, sigma, psopts{:});
                D = createSparseDistanceMatrix([pstruct.x' pstruct.y'], [x y; xn yn], R); % search radius < min distance
                [linkDet2Ref, ~] = lap(D, [], [], 1);
                linkDet2Ref = double(linkDet2Ref(1:numel(pstruct.x)));
                ntpPSf{ai}(k,i) = sum(linkDet2Ref<=np);
                nfpPSf{ai}(k,i) = sum(linkDet2Ref>np);
            end
        end
    end
    save PanelC;
end


fset = loadFigureSettings('print');
hues = [0.33 0.55 0.15 0];
cv = hsv2rgb([hues' ones(na,2)]);

figure(fset.fOpts{:}, 'Position', [2 2 8 9.5]);

% grid below data
axes(fset.axOpts{:}, 'Position', [1.5 5 6 3]);
set(gca, 'xscale', 'log');
grid on;
axis([1 100 0 1.02]);


set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', '-', 'XTickLabel', [], 'YTickLabel', [],...
    'XColor', 0.75*[1 1 1], 'YColor', 0.75*[1 1 1]);

% plot SNR curves
axes(fset.axOpts{:}, 'Position', [1.5 5 6 3]);
hold on;
set(gca, 'xscale', 'log', 'Color', 'none');
axis([1 100 0 1.02]);
  

% total points/SNR level: nrep
hp = zeros(na,1);
for ai = 1:na
    hp(ai) = plot(psnrV, nanmean(ntpPS{ai},2)/nrep, 'Color', cv(ai,:), 'LineWidth', 1);
end
set(gca, 'XTickLabel', [1 10 100], 'XTickLabel', []);
hy(1) = ylabel('Detection sensitivity', fset.lfont{:});
% hl = legend(hp, ' PSF model', ' u-track', ' Wavelet', ' Imaris', 'Location', 'NorthEast');
% set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [8.15 1.45 1.25 1.5]);
formatTickLabels();


axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 3]);
set(gca, 'xscale', 'log');
grid on;
axis([1 100 0 0.25]);

set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', '-', 'XTickLabel', [], 'YTickLabel', [],...
    'XColor', 0.75*[1 1 1], 'YColor', 0.75*[1 1 1]);


% plot SNR curves
axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 3]);
hold on;
set(gca, 'xscale', 'log', 'Color', 'none');
axis([1 100 0 1]);

for ai = 1:na

    tp = nanmedian(ntpPSf{ai},2);
    fp = nanmedian(nfpPSf{ai},2);
    plot(psnrV, fp./(tp+fp), '-', 'Color', cv(ai,:), 'LineWidth', 1);

end
set(gca, 'XTickLabel', [1 10 100], 'Layer', 'bottom');
xlabel('PSNR', fset.lfont{:});
hy(2) = ylabel('False discovery rate', fset.lfont{:});
% hl = legend(hp, arrayfun(@(i) [' \alpha = ' num2str(i)], alphaV, 'unif', 0));
% set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [5.5 5.25 1.25 1.5]);
formatTickLabels();

% align ylabels
ypos2 = get(hy(2), 'Position');
ypos1 = get(hy(1), 'Position');
set(hy(1), 'Position', [ypos2(1) ypos1(2:end)]);




wg = ceil(6*sigma);
[xg,yg] = meshgrid(-wg:wg);
% PSF model
g0 = exp(-(xg.^2+yg.^2)/(2*sigma^2));
n = numel(g0);

% panels with different SNRs
dx = 0.06;
pw = 0.95;

snrInset = [1 3 6 10 30 100];
g_ex = cell(1,numel(snrInset));
for i = 1:numel(snrInset)
    g_ex{i} = g0 + sqrt(1^2 * n/(n-1) / snrInset(i)) * randn(size(g0));
end
gcat = cat(3, g_ex{:});


for i = 1:numel(snrInset)
    axes(fset.axOpts{:}, 'Position', [1.5+(i-1)*(dx+pw) 8.2 pw pw]);
    imagesc(g_ex{i}); colormap(gray(256)); axis image off;
    text(size(g0,1),0, num2str(snrInset(i)), 'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    if i==1
       text(-1,0, 'PSNR', 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.sfont{:});
    end
end
% print('-depsc2', '-loose', 'FigureS2_PanelC.eps');


