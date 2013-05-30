% Francois Aguet, October 2010 (last modified: 10/09/2012)

% To do: apply clustering algorithm in place of global analysis

function runSlaveChannelClassification(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.addParamValue('np', 10000);
ip.addParamValue('Mode', 'random', @(x) any(strcmpi(x, {'random', 'maskRatio'})));
ip.parse(data, varargin{:});

% reset random number generator to ensure reproducibility
rng('default');
fprintf('Random number generator set to defaults by ''runSlaveClassification()''.\n');

for i = 1:length(data)
    fprintf('Running slave ch. classification for %s\n', getShortPath(data(i)));
    main(data(i), ip.Results);
end


function main(data, opts)

kLevel = norminv(1-opts.Alpha/2.0, 0, 1); % ~2 std above background for 0.05

% load tracks (all)
ts = load([data.source 'Tracking' filesep opts.FileName]);

if isfield(ts.tracks, 'significantSignal') && ~opts.Overwrite
    fprintf('Classification has already been run for %s\n', getShortPath(data));
    return
end

% load cell mask
cellmask = logical(getCellMask(data));

% Determine master/slave channels
nc = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);
sCh = setdiff(1:nc, mCh);

load([data.source 'Detection' filesep 'detection_v2.mat']);
sigma = frameInfo(1).s;
w = max(ceil(4*sigma(sCh)));

frameIdx = round(linspace(1, data.movieLength, 16));
nf = numel(frameIdx);

nx = data.imagesize(2);
ny = data.imagesize(1);

bgA = cell(1,nf);
pSlaveSignal = cell(1,nf);

parfor f = 1:nf;
    k = frameIdx(f);
    
    %-----------------------------------------------
    % Generate masks
    %-----------------------------------------------
    % load CCP mask and dilate
    ccpMask = double(imread(data.maskPaths{k}));
    ccpMask(ccpMask~=0) = 1;
    ccpMask = imdilate(ccpMask, strel('disk', 1*w));
    
    
    %-----------------------------------------------
    % Probability of randomly occurring slave signal
    %-----------------------------------------------
    switch opts.Mode
        case 'random'
            %=================================================================================
            % Approach 1: fit at random positions outside CCPs, build distribution of 'A'
            %=================================================================================
            
            %-----------------------------------------
            % Distribution of points within cell
            %-----------------------------------------
            % generate candidate points
            xa = (nx-2*w-1)*rand(1,opts.np)+w+1;
            ya = (ny-2*w-1)*rand(1,opts.np)+w+1;
            xi = round(xa);
            yi = round(ya);
            
            % remove points outside of mask or within border
            linIdxA = sub2ind([ny nx], yi, xi);
            rmIdx = cellmask(linIdxA)==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
            xa(rmIdx) = [];
            ya(rmIdx) = [];
            linIdxA(rmIdx) = [];
            
            %-----------------------------------------
            % Distribution of points outside CCSs
            %-----------------------------------------
            mask = cellmask - ccpMask;
            
            % generate candidate points
            x = [];
            y = [];
            linIdx = [];
            while numel(x)<opts.np
                xcand = (nx-2*w-1)*rand(1,opts.np)+w+1;
                ycand = (ny-2*w-1)*rand(1,opts.np)+w+1;
                xi = round(xcand);
                yi = round(ycand);
                
                % remove points outside of mask or within border
                linIdxCand = sub2ind([ny nx], yi, xi);
                validIdx = mask(linIdxCand)==1 & xi>w & yi>w & xi<=nx-w & yi<=ny-w;
                x = [x xcand(validIdx)];
                y = [y ycand(validIdx)];
                linIdx = [linIdx linIdxCand(validIdx)];
            end
            x = x(1:opts.np);
            y = y(1:opts.np);
            linIdx = linIdx(1:opts.np);
            
            
            bgA{f} = NaN(nc,opts.np);
            pSlaveSignal{f} = NaN(nc,numel(xa));
            for c = sCh
                frame = double(imread(data.framePaths{c}{k}));
                
                % get local min & max for initial c and A
                ww = 2*ceil(4*sigma(c))+1;
                maxF = ordfilt2(frame, ww^2, true(ww));
                minF = ordfilt2(frame, 1, true(ww));
                
                pstruct = fitGaussians2D(frame, xa, ya, maxF(linIdxA), sigma(c), minF(linIdxA), 'Ac');
                pSlaveSignal{f}(c,:) = sum(pstruct.pval_Ar < 0.05) / numel(xa);
                
                % estimate background amplitude
                pstruct = fitGaussians2D(frame, x, y, maxF(linIdx), sigma(c), minF(linIdx), 'Ac');
                bgA{f}(c,:) = pstruct.A;
            end

%         case 'maskRatio'
%             %=================================================================================
%             % Approach 2: compute mask of significant pixels in slave channel (index 2)
%             %=================================================================================
%             % Note: the following code is adapted from pointSourceDetection.m
%             
%             % Gaussian kernel
%             x = -w:w;
%             g = exp(-x.^2/(2*sigma^2));
%             u = ones(1,length(x));
%             
%             % convolutions
%             imgXT = padarrayXT(frame, [w w], 'symmetric');
%             fg = conv2(g', g, imgXT, 'valid');
%             fu = conv2(u', u, imgXT, 'valid');
%             fu2 = conv2(u', u, imgXT.^2, 'valid');
%             
%             % Laplacian of Gaussian
%             gx2 = g.*x.^2;
%             imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
%             imgLoG = imgLoG / (2*pi*sigma^2);
%             
%             % 2-D kernel
%             g = g'*g;
%             n = numel(g);
%             gsum = sum(g(:));
%             g2sum = sum(g(:).^2);
%             
%             % solution to linear system
%             A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
%             c_est = (fu - A_est*gsum)/n;
%             
%             J = [g(:) ones(n,1)]; % g_dA g_dc
%             C = inv(J'*J);
%             
%             f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
%             RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
%             sigma_e2 = RSS/(n-3);
%             
%             sigma_A = sqrt(sigma_e2*C(1,1));
%             
%             % standard deviation of residuals
%             sigma_res = sqrt((RSS - (A_est*gsum+n*c_est - fu)/n)/(n-1));
%             
%             SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
%             df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
%             scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
%             T = (A_est - sigma_res*kLevel) ./ scomb;
%             pval = tcdf(-T, df2);
%             
%             % mask of admissible positions for local maxima
%             mask = pval < 0.05;
%             
%             % all local max
%             allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);
%             
%             % local maxima above threshold in image domain
%             imgLM = allMax .* mask;
%             
%             if sum(imgLM(:))~=0 % no local maxima found, likely a background image
%                 
%                 % -> set threshold in LoG domain
%                 logThreshold = min(imgLoG(imgLM~=0));
%                 logMask = imgLoG >= logThreshold;
%                 
%                 % combine masks
%                 mask = mask | logMask;
%             end
%             %---------------------------------------------------------------------------------
%             % Note: end of pointSourceDetection.m code
%             %---------------------------------------------------------------------------------
%             
%             %=================================================================================
%             % Ratio between significant pixels in slave channel and master channel
%             %=================================================================================
%             %mask = mask & ccpMask; % areas of significant slave signal within EAZ
%             %pSlaveSignal(i) = sum(mask(:)) / sum(ccpMask(:));
%             
%             pSlaveSignal(i) = sum(mask(:)) / sum(cellmask(:));
    end
end
pSlave = [pSlaveSignal{:}];
pSlave = mean(pSlave,2);
for c = sCh
    fprintf('P(random detection in ch. %d) = %.3f\n', c, pSlave(c));
end
%=================================================================================
% Classify tracks in slave channels
%=================================================================================

% steps for required to reject H_0: binomial
% nBinSteps = ceil(log(alpha)/log(0.5));

bg95 = prctile([bgA{:}], 95, 2);

% Loops through all the tracks
nt = numel(ts.tracks);
for k = 1:nt
    
    np = numel(ts.tracks(k).t); % # points/track
    ts.tracks(k).isDetected = NaN(nc, np);
    ts.tracks(k).significantMaster = NaN(nc,1);
    ts.tracks(k).significantVsBackground = NaN(nc,np);
    ts.tracks(k).significantSlave = NaN(nc,1);

    for c = sCh % loop through all slave channels
        
        % significance test, binarization
        npx = round((ts.tracks(k).sigma_r(c,:) ./ ts.tracks(k).SE_sigma_r(c,:)).^2/2+1);
        A = ts.tracks(k).A(c,:);
        sigma_A = ts.tracks(k).A_pstd(c,:);
        
        % significance test for independent detecion
        sigma_r = ts.tracks(k).sigma_r(c,:) * kLevel;
        SE_sigma_r = ts.tracks(k).SE_sigma_r(c,:) * kLevel;
        df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
        scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
        T = (A - sigma_r) ./ scomb;
        pval = tcdf(-T, df2);
        ts.tracks(k).isDetected(c,:) = pval < 0.05;
        
        % test whether # significant points > 95th percentile of 'random' distribution
        ts.tracks(k).significantMaster(c) = nansum(ts.tracks(k).isDetected(c,:)) > binoinv(0.95, np, pSlave(c));
        
        % significance test relative to background slave signal
        sigma_r = bg95(c);
        SE_sigma_r = sigma_r ./ sqrt(2*npx-1);
        df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
        scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
        T = (A - sigma_r) ./ scomb;
        pval = tcdf(-T, df2);
        ts.tracks(k).significantVsBackground(c,:) = pval < 0.05;
        ts.tracks(k).significantSlave(c) = nansum(ts.tracks(k).significantVsBackground(c,:)) > binoinv(0.95, np, 0.05); %%%%%%%%%%%%%%

        % Criterion based on length of significant regions
        % posLengths = find(diff([bd 0])==-1) - find(diff([0 bd 0])==1) + 1;
    end
end

idx = [ts.tracks.catIdx]==1 & [ts.tracks.lifetime_s]>=data.framerate*opts.Cutoff_f;
nPosM = sum([ts.tracks(idx).significantMaster],2);
nPosS = sum([ts.tracks(idx).significantSlave],2);
for c = sCh
    fprintf('Ch. %d positive tracks as master: %.2f %% (%d/%d valid, %d total)\n', c, 100*nPosM(c)/sum(idx), nPosM(c), sum(idx), nt);
    fprintf('Ch. %d positive tracks as slave:  %.2f %% (%d/%d valid, %d total)\n', c, 100*nPosS(c)/sum(idx), nPosS(c), sum(idx), nt);
end

%=================================================================================
% Determine whether disappearance of slave channel signal correlates with master
%=================================================================================
% conditions:
% - signal in last 5 frames of the track must be significant (1 gap allowed)
% - binary signals: last 5 points of track, first 2 points of buffer correlate up to 1 point (= 1 gap allowed)
% - normalized correlation must be > 0.8

% for k = 1:nt
%     c = sCh(1);
%     
%     if ts.tracks(k).catIdx==1 && numel(ts.tracks(k).t)>=5
%         % binary classification (hval_Ar==1, pval_Ar<0.05 if significant signal)
%         bc = [ts.tracks(k).hval_Ar(mCh,end-4:end) == ts.tracks(k).hval_Ar(c,end-4:end)...
%             (ts.tracks(k).endBuffer.pval_Ar(mCh,1:2)<0.05)==(ts.tracks(k).endBuffer.pval_Ar(c,1:2)<0.05)];
%         
%         mEnd = [ts.tracks(k).A(mCh,end-4:end) ts.tracks(k).endBuffer.A(mCh,1:2)];
%         sEnd = [ts.tracks(k).A(sCh,end-4:end) ts.tracks(k).endBuffer.A(sCh,1:2)];
%         mEnd = mEnd/max(mEnd);
%         sEnd = sEnd/max(sEnd);
%         K = sum(mEnd.*sEnd)/sqrt(sum(mEnd.^2)*sum(sEnd.^2));
%         
%         ts.tracks(k).corrDisappearance = (sum(bc)>=6) && K>0.8 &&...
%             sum(ts.tracks(k).isDetected(mCh,end-4:end))>=4 && sum(ts.tracks(k).isDetected(sCh,end-4:end))>=4;
%     end
% end

save([data.source 'Tracking' filesep opts.FileName], '-struct', 'ts');

% update lifetime data structures
getLifetimeData(data, 'Overwrite', true);
