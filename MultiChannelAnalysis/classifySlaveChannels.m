% Francois Aguet, October 2010 (last modified: 06/30/2011)

function tracks = classifySlaveChannels(data, tracks)

% Significance thresholds
alpha = 0.05;
% sigmaT = icdf('normal', 1-alpha/2, 0, 1); % ~2 sigma
sigmaT = icdf('normal', 1-alpha, 0, 1); % weaker, single-tailed (~1.6 sigma)


%=================================
% Determine master/slave channels
%=================================
nCh = length(data.channels); % number of channels
masterChannel = strcmp(data.source, data.channels);
slaveChannels = setdiff(1:nCh, masterChannel);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{masterChannel});
w = ceil(4*sigma);

%===============================================
% Classification
%===============================================


%-----------------------------------------------
% Utilization of a generic mask is not sufficient:
%-----------------------------------------------
% define mask for localization
bgMask = ones(2*w+1);
[xg,yg] = meshgrid(-w:w);
r = sqrt(xg.^2+yg.^2);
bgMask(r<=2*sigma) = 0;
% figure; imagesc(bgMask); colormap(gray(256)); axis image; colorbar;




%-----------------------------------------------
% 1. Generate mask for background points
%-----------------------------------------------
maskList = dir([data.source 'Detection' filesep 'Masks' filesep '*.tif']);

i = 1;
for k = 1:10:data.movieLength
    
    % load mask, dilate, -mask
    mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(k).name]));
    mask(mask~=0) = 1;
    %dmask = imdilate(mask, strel('disk', 2*w)) - imdilate(mask, strel('disk', w));
    dmask = imdilate(mask, strel('disk', 1*w));
    dmask(dmask~=0) = 1;
    %dmask = dmask - mask;
    %dmask(dmask<0) = 0;
    
    % % in first frame, determine proportion of positive detections in background
    %
    % % Mask excluding CCPs in current frame
    % mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(10).name]));
    % mask(mask~=0) = 1;
    % mask = maskInt-imdilate(mask, strel('disk', w));
    % mask(mask<0) = 0;
    % % figure; imagesc(mask); colormap(gray(256)); axis image;
    %
    %
    % npRef = 3000;
    % % generate candidate points
    % x = (nx-2*w-1)*rand(1,npRef)+w+1;
    % y = (ny-2*w-1)*rand(1,npRef)+w+1;
    % xi = round(x);
    % yi = round(y);
    %
    % % remove points outside of mask or within border
    % idx = dmask(sub2ind([ny nx], yi,xi))==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
    % x(idx) = [];
    % y(idx) = [];
    % xi(idx) = [];
    % yi(idx) = [];
    % npRef = length(x);
    
    
    %
    % frameList = dir([data.channels{2} '*.tif*']);
    % frame = double(imread([data.channels{2} frameList(1).name]));
    % ctrlStatus = zeros(1, npRef);
    % for p = 1:npRef
    %     window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
    %     c = mean(window(lmask==1));
    %     cStd = std(window(lmask==1));
    %     prm = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) max(window(:))-c sigma c], 'xyA');
    %     A = prm(3);
    %     ctrlStatus(p) = A > sigmaT*cStd;
    % end
    % fprintf('Ch. %d background detections: %d/%d (%.2f %%)\n', 2, sum(ctrlStatus), npRef, 100*sum(ctrlStatus)/npRef);
    
    
    %-----------------------------------------------
    % Pick random mask from detection
    %-----------------------------------------------
    load([data.source 'Detection' filesep 'detectionResults.mat']);
    
    frameList = dir([data.channels{2} '*.tif*']);
    frame = double(imread([data.channels{2} frameList(10).name]));
    
    % ctrlStatus = zeros(1, npRef);
    N = sum(bgMask(:));
    
    %=========================================================================
    % Note: the following code is adapted from pointSourceDetection.m
    %=========================================================================
    
    % Gaussian kernel
    x = -w:w;
    g = exp(-x.^2/(2*sigma^2));
    u = ones(1,length(x));
    
    % convolutions
    imgXT = padarrayXT(frame, [w w], 'symmetric');
    fg = conv2(g', g, imgXT, 'valid');
    fu = conv2(u', u, imgXT, 'valid');
    fu2 = conv2(u', u, imgXT.^2, 'valid');
    
    % Laplacian of Gaussian
    %gx2 = g.*x.^2;
    %imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
    %imgLoG = imgLoG / (2*pi*sigma^2);
    
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
    
    kLevel = norminv(1-alpha/2.0, 0, 1);
    
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;
    pval = tcdf(real(T), df2);
    
    % mask of admissible positions for local maxima
    mask = pval > 0.95;
    
    % old threshold/test: background at each pixel
    %     E = conv2(padarray(frame, [w w], 'replicate'), bgMask/N, 'valid');
    %     E2 = conv2(padarray((frame).^2, [w w], 'replicate'), bgMask/N, 'valid');
    %     cStdMap = sqrt(N/(N-1) * (E2-E.^2));
    %     cMap = E;
    %
    %     % amplitude map
    %     gL2 = pi*sigma^2;
    %     g = exp(-(-w:w).^2/(2*sigma^2));
    %     aMap = conv2(g', g, padarray(frame, [w w], 'replicate'), 'valid');
    %     aMap = (aMap-2*pi*sigma^2*cMap) / gL2;
    %
    %     % figure; imagesc(aMap); colormap(gray(256)); axis image; colorbar;
    %     testMap = aMap > sigmaT * cStdMap;
    
    %     p(i) = sum(testMap(:)) / sum(dmask(:));
    p(i) = sum(mask(:)) / sum(dmask(:));
    i = i+1;
end
p = mean(p);


% for p = 1:npRef
%     window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
%     % pick random mask
%     mi = unidrnd(length(xm));
%     maskWindow = ~mask(ymi(mi)-w:ymi(mi)+w, xmi(mi)-w:xmi(mi)+w);
%
%     % background estimate
%     c = mean(window(maskWindow));
%     cStd = std(window(maskWindow));
%
%     prm = fitGaussian2D(window, [0 0 max(window(:))-c sigma c], 'xyA');
%     %prm = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) max(window(:))-c sigma c], 'A');
%     %prm = fitGaussian2D(window, [xm(mi)-xmi(mi) ym(mi)-ymi(mi) max(window(:))-c sigma c], 'A');
%     A = prm(3);
%     ctrlStatus(p) = A > sigmaT*cStd;
%
%     if A > sigmaT*cStd
%         prm
%         figure; imagesc(window); colormap(gray(256)); axis image;
%         hold on;
%         plot(xi(p) + prm(1), yi(p) + prm(2), 'rx');
%         figure; imagesc(maskWindow); colormap(gray(256)); axis image;
%         return
%     end
%
%
% end
% fprintf('Ch. %d background detections: %d/%d (%.2f %%)\n', 2, sum(ctrlStatus), npRef, 100*sum(ctrlStatus)/npRef);


%=================================
% Proportion of 'random' positives
%=================================

%=================================
% Classify tracks in slave channels
%=================================

% steps for required to reject H_0: binomial
% nBinSteps = ceil(log(alpha)/log(0.5));

% Loops through all the tracks
for k = 1:length(tracks);
    
    nt = numel(tracks(k).t);
    tracks(k).A_binary = NaN(nCh, nt);
    %tracks(k).signalSignificance = NaN(nCh,1);
    tracks(k).significantSignal = NaN(nCh,1);
    tracks(k).significantSignal_buffered = NaN(nCh,1);
    
    for c = 1:nCh % loop through all channels
        tracks(k).A_binary(c,:) = tracks(k).A(c,:) > sigmaT * tracks(k).sigma_r(c,:);
        % To do: apply clustering algorithm in place of global analysis
        % posLengths = find(diff([bd 0])==-1) - find(diff([0 bd 0])==1) + 1;
        
        if sum(tracks(k).A_binary(c,:)) > binoinv(0.95, nt, p);
            tracks(k).significantSignal(c) = 1;
        else
            tracks(k).significantSignal(c) = 0;
        end
        
        if isfield(tracks(k), 'startBuffer') && isfield(tracks(k), 'endBuffer') && ~isempty(tracks(k).startBuffer) && ~isempty(tracks(k).endBuffer)
            tracks(k).startBuffer.A_binary(c,:) = tracks(k).startBuffer.A(c,:) > sigmaT * tracks(k).startBuffer.sigma_r(c,:);
            tracks(k).endBuffer.A_binary(c,:) = tracks(k).endBuffer.A(c,:) > sigmaT * tracks(k).endBuffer.sigma_r(c,:);
            
            % test extended track
            A_binary_XT = [tracks(k).startBuffer.A_binary(c,:) tracks(k).A_binary(c,:) tracks(k).endBuffer.A_binary(c,:)];
            if sum(A_binary_XT) > binoinv(0.95, numel(A_binary_XT), p)
                tracks(k).significantSignal_buffered(c) = 1;
            else
                tracks(k).significantSignal_buffered(c) = 0;
            end
            
            % test each buffer individually
            if sum(tracks(k).startBuffer.A_binary(c,:)) > binoinv(0.95, size(tracks(k).startBuffer.A,2), p)
                tracks(k).startBuffer.significantSignal(c) = 1;
            else
                tracks(k).startBuffer.significantSignal(c) = 0;
            end
            if sum(tracks(k).endBuffer.A_binary(c,:)) > binoinv(0.95, size(tracks(k).endBuffer.A,2), p)
                tracks(k).endBuffer.significantSignal(c) = 1;
            else
                tracks(k).endBuffer.significantSignal(c) = 0;
            end
        end
    end
end

for c = slaveChannels
    nPos = sum(arrayfun(@(x) x.significantSignal(c), tracks));
    fprintf('Ch. %d positive tracks: %d/%d (%.2f %%)\n', c, nPos, length(tracks), 100*nPos/length(tracks));
end



% Cell mask: sum of first 10 frames
% frameInt = zeros(data.imagesize);
% for k = 1:10
%     frame = double(imread([data.source frameList(k).name]));
%     frameInt = frameInt + frame;
% end

% Sum all masks, generate cell mask
% maskInt = zeros(data.imagesize);
% for k = 1:data.movieLength
%     mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(k).name]));
%     maskInt = maskInt + mask;
% end
% maskInt(maskInt~=0) = 1;
% maskInt = imclose(maskInt, strel('disk', 2*w));
% maskInt = ~imclose(~maskInt, strel('disk', 2*w));
%
%
%