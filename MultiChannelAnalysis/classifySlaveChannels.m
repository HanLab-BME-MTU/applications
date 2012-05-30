% Francois Aguet, October 2010 (last modified: 05/29/2012)

% To do: apply clustering algorithm in place of global analysis


function tracks = classifySlaveChannels(data, tracks)

% Significance thresholds
alpha = 0.05;
kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background
% kLevel = norminv(1-alpha, 0, 1); % ~1.64 std above background


% load cell mask
cellmask = logical(getCellMask(data));

%=================================
% Determine master/slave channels
%=================================
nCh = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);
sCh = setdiff(1:nCh, mCh);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});
w = ceil(4*sigma);

%===============================================
% Classification
%===============================================

%-----------------------------------------------
% 1. Generate mask for background points
%-----------------------------------------------
frameIdx = 1:20:data.movieLength;
nf = numel(frameIdx);
pSlaveSignal = zeros(1,nf);
for i = 1:nf;
    k = frameIdx(i);
    % load mask and dilate
    dmask = double(imread(data.maskPaths{k}));
    dmask(dmask~=0) = 1;
    dmask = imdilate(dmask, strel('disk', 1*w));

    % % in first frame, determine proportion of positive detections in background
    %
    % % Mask excluding CCPs in current frame
    % mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(10).name]));
    % mask(mask~=0) = 1;
    % mask = maskInt-imdilate(mask, strel('disk', w));
    % mask(mask<0) = 0;
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
    
    % frame = double(imread(data.framePaths{2}{1}));
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
    
    
    % Load slave channel frame
    frame = double(imread(data.framePaths{2}{k}));
        
    %---------------------------------------------------------------------------------
    % Note: the following code is adapted from pointSourceDetection.m
    % Here, only a mask of significant pixels in the slave channel is computed
    %---------------------------------------------------------------------------------
    
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
    
    kLevel = norminv(1-alpha/2.0, 0, 1);
    
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;
    pval = tcdf(-T, df2);
    
    % mask of admissible positions for local maxima
    mask = pval < 0.05;
    
    % all local max
    allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);
    
    % local maxima above threshold in image domain
    imgLM = allMax .* mask;
    
    if sum(imgLM(:))~=0 % no local maxima found, likely a background image
        
        % -> set threshold in LoG domain
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;
        
        % combine masks
        mask = mask | logMask;
    end
    
    %=================================================================================
    % Ratio between significant pixels in slave channel and XXX
    %=================================================================================
    dmask = dmask & cellmask; % defines endocytically active zone (EAZ)
    mask = mask & dmask; % areas where significant slave signal was detected within EAZ
    
    pSlaveSignal(i) = sum(mask(:)) / sum(dmask(:));
end
pSlaveSignal = mean(pSlaveSignal);


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

%=================================================================================
% Classify tracks in slave channels
%=================================================================================

% steps for required to reject H_0: binomial
% nBinSteps = ceil(log(alpha)/log(0.5));

% Loops through all the tracks
for k = 1:length(tracks);
    
    nt = numel(tracks(k).t);
    tracks(k).A_binary = NaN(nCh, nt);
    tracks(k).significantSignal = NaN(nCh,1);
    %tracks(k).significantSignal_buffered = NaN(nCh,1);
    
    for c = 1:nCh % loop through all channels
        
        % significance test, binarization
        npx = round((tracks(k).sigma_r(c,:) ./ tracks(k).SE_sigma_r(c,:)).^2/2+1);
        A = tracks(k).A(c,:);
        sigma_A = tracks(k).A_pstd(c,:);
        sigma_r = tracks(k).sigma_r(c,:) * kLevel;
        SE_sigma_r = tracks(k).SE_sigma_r(c,:) * kLevel;
                
        df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
        scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
        T = (A - sigma_r) ./ scomb;
        pval = tcdf(-T, df2);
        tracks(k).A_binary(c,:) = pval < alpha;

        % test whether # significant points > 95th percentile of 'random' distribution
        tracks(k).significantSignal(c) = nansum(tracks(k).A_binary(c,:)) > binoinv(0.95, nt, pSlaveSignal);

        % Criterion based on length of significant regions
        % posLengths = find(diff([bd 0])==-1) - find(diff([0 bd 0])==1) + 1;
    end
end

for c = sCh
    nPos = sum(arrayfun(@(x) x.significantSignal(c), tracks));
    fprintf('Ch. %d positive tracks: %d/%d (%.2f %%)\n', c, nPos, length(tracks), 100*nPos/length(tracks));
end

%=================================================================================
% Determine whether disappearance of slave channel signal correlates with master
%=================================================================================
% lambda = 20;

% normalize master signal to [0..1]

% conditions:
% - signal in last 5 frames of the track must be significant (1 gap allowed)
% - binary signals: last 5 points of track, first 2 points of buffer correlate up to 1 point (= 1 gap allowed)
% - normalized correlation must be > 0.8
for k = 1:length(tracks);
    c = sCh(1);
    
    if tracks(k).catIdx==1 && numel(tracks(k).t)>4
        % binary classification (H==1 if significant signal)
        bc = [tracks(k).hval_Ar(mCh,end-4:end) == tracks(k).hval_Ar(c,end-4:end)...
            (tracks(k).endBuffer.pval_Ar(mCh,1:2)<0.05)==(tracks(k).endBuffer.pval_Ar(c,1:2)<0.05)];
        
        mEnd = [tracks(k).A(mCh,end-4:end) tracks(k).endBuffer.A(mCh,1:2)];
        sEnd = [tracks(k).A(sCh,end-4:end) tracks(k).endBuffer.A(sCh,1:2)];
        mEnd = mEnd/max(mEnd);
        sEnd = sEnd/max(sEnd);
        K = sum(mEnd.*sEnd)/sqrt(sum(mEnd.^2)*sum(sEnd.^2));
        
        tracks(k).corrDisappearance = (sum(bc)>=6) && K>0.8 &&...
            sum(tracks(k).A_binary(mCh,end-4:end))>=4 && sum(tracks(k).A_binary(sCh,end-4:end))>=4;
    end
end
