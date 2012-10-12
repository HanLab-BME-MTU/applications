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
ip.addParamValue('Mode', 'maskRatio', @(x) any(strcmpi(x, {'random', 'maskRatio'})));
ip.parse(data, varargin{:});

% reset random number generator to assure reproducibility
rng('default');
fprintf('Random number generator set to defaults by ''runSlaveClassification()''.\n');

for i = 1:length(data)
    fprintf('Running slave ch. classification for %s\n', getShortPath(data(i)));
    main(data(i), ip.Results.Alpha, ip.Results.Overwrite, ip.Results.Mode, ip.Results.FileName, ip.Results.Cutoff_f);
end
    

function main(data, alpha, overwrite, mode, fileName, cutoff_f)

kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background
% kLevel = norminv(1-alpha, 0, 1); % ~1.64 std above background

% load tracks (all)
ts = load([data.source 'Tracking' filesep fileName]);

if ~isfield(ts.tracks, 'significantSignal') || overwrite
    
    % load cell mask
    cellmask = logical(getCellMask(data));
    
    % Determine master/slave channels
    nCh = length(data.channels); % number of channels
    mCh = strcmp(data.source, data.channels);
    sCh = setdiff(1:nCh, mCh);
    
    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{sCh});
    w = ceil(4*sigma);

    frameIdx = round(linspace(1, data.movieLength, 16));
    nf = numel(frameIdx);
    pSlaveSignal = zeros(1,nf);
    
    parfor i = 1:nf;
        k = frameIdx(i);
        
        %-----------------------------------------------
        % Generate masks
        %-----------------------------------------------
        % load CCP mask and dilate
        ccpMask = double(imread(data.maskPaths{k}));
        ccpMask(ccpMask~=0) = 1;
        ccpMask = imdilate(ccpMask, strel('disk', 1*w));
        % mask = maskInt-imdilate(mask, strel('disk', w));
        ccpMask = ccpMask & cellmask; % endocytically active zone (EAZ)
        
        frame = double(imread(data.framePaths{2}{k}));
        [ny,nx] = size(frame);
        %-----------------------------------------------
        % Probability of randomly occurring slave signal
        %-----------------------------------------------
        switch mode
            case 'random'
                %=================================================================================
                % Approach 1: fit at random positions outside CCPs, build distribution of 'A'
                %=================================================================================
                
                % Mask excluding CCPs in current frame
                %xmask = cellMask - ccpMask;
                
                npRef = 50000;
                % generate candidate points
                x = (nx-2*w-1)*rand(1,npRef)+w+1;
                y = (ny-2*w-1)*rand(1,npRef)+w+1;
                xi = round(x);
                yi = round(y);
                
                % remove points outside of mask or within border
                linIdx = sub2ind([ny nx], yi, xi);
                rmIdx = cellmask(linIdx)==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
                x(rmIdx) = [];
                y(rmIdx) = [];
                linIdx(rmIdx) = [];
                npRef = numel(x);
                
                % get local min & max for initial c and A
                ww = 2*w+1;
                maxF = ordfilt2(frame, ww^2, true(ww));
                minF = ordfilt2(frame, 1, true(ww));
                
                pStruct = fitGaussians2D(frame, x, y, maxF(linIdx), sigma, minF(linIdx), 'Ac');
                pSlaveSignal(i) = sum(pStruct.pval_Ar < 0.05) / npRef;
           
            case 'maskRatio'
                %=================================================================================
                % Approach 2: compute mask of significant pixels in slave channel (index 2)
                %=================================================================================
                % Note: the following code is adapted from pointSourceDetection.m
                
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
                %---------------------------------------------------------------------------------
                % Note: end of pointSourceDetection.m code
                %---------------------------------------------------------------------------------
                
                %=================================================================================
                % Ratio between significant pixels in slave channel and master channel
                %=================================================================================
                %mask = mask & ccpMask; % areas of significant slave signal within EAZ
                %pSlaveSignal(i) = sum(mask(:)) / sum(ccpMask(:));
                
                pSlaveSignal(i) = sum(mask(:)) / sum(cellmask(:));
        end
    end
    pSlaveSignal = mean(pSlaveSignal);
    fprintf('P(random slave detection) = %.3f\n', pSlaveSignal);
    
    %=================================================================================
    % Classify tracks in slave channels
    %=================================================================================
    
    % steps for required to reject H_0: binomial
    % nBinSteps = ceil(log(alpha)/log(0.5));
    
    % Loops through all the tracks
    nt = numel(ts.tracks);
    for k = 1:nt
        
        np = numel(ts.tracks(k).t); % # points/track
        ts.tracks(k).A_binary = NaN(nCh, np);
        ts.tracks(k).significantSignal = NaN(nCh,1);
        
        for c = 1:nCh % loop through all channels
            
            % significance test, binarization
            npx = round((ts.tracks(k).sigma_r(c,:) ./ ts.tracks(k).SE_sigma_r(c,:)).^2/2+1);
            A = ts.tracks(k).A(c,:);
            sigma_A = ts.tracks(k).A_pstd(c,:);
            sigma_r = ts.tracks(k).sigma_r(c,:) * kLevel;
            SE_sigma_r = ts.tracks(k).SE_sigma_r(c,:) * kLevel;
            
            df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
            scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
            T = (A - sigma_r) ./ scomb;
            pval = tcdf(-T, df2);
            ts.tracks(k).A_binary(c,:) = pval < alpha;
            
            % test whether # significant points > 95th percentile of 'random' distribution
            ts.tracks(k).significantSignal(c) = nansum(ts.tracks(k).A_binary(c,:)) > binoinv(0.95, np, pSlaveSignal);
            
            % Criterion based on length of significant regions
            % posLengths = find(diff([bd 0])==-1) - find(diff([0 bd 0])==1) + 1;
        end
    end
    
    idx = [ts.tracks.catIdx]==1 & [ts.tracks.lifetime_s]>=data.framerate*cutoff_f;
    nPos = sum([ts.tracks(idx).significantSignal],2);
    for c = sCh
        %fprintf('Ch. %d positive tracks: %d/%d (%.2f %%)\n', c, nPos(c), nt, 100*nPos(c)/nt);
        fprintf('Ch. %d positive tracks: %.2f %% (%d/%d valid, %d total)\n', c, 100*nPos(c)/nPos(mCh), nPos(c), nPos(mCh), nt);
    end
    
    %=================================================================================
    % Determine whether disappearance of slave channel signal correlates with master
    %=================================================================================
    % conditions:
    % - signal in last 5 frames of the track must be significant (1 gap allowed)
    % - binary signals: last 5 points of track, first 2 points of buffer correlate up to 1 point (= 1 gap allowed)
    % - normalized correlation must be > 0.8
    for k = 1:nt
        c = sCh(1);
        
        if ts.tracks(k).catIdx==1 && numel(ts.tracks(k).t)>=5
            % binary classification (hval_Ar==1, pval_Ar<0.05 if significant signal)
            bc = [ts.tracks(k).hval_Ar(mCh,end-4:end) == ts.tracks(k).hval_Ar(c,end-4:end)...
                (ts.tracks(k).endBuffer.pval_Ar(mCh,1:2)<0.05)==(ts.tracks(k).endBuffer.pval_Ar(c,1:2)<0.05)];
            
            mEnd = [ts.tracks(k).A(mCh,end-4:end) ts.tracks(k).endBuffer.A(mCh,1:2)];
            sEnd = [ts.tracks(k).A(sCh,end-4:end) ts.tracks(k).endBuffer.A(sCh,1:2)];
            mEnd = mEnd/max(mEnd);
            sEnd = sEnd/max(sEnd);
            K = sum(mEnd.*sEnd)/sqrt(sum(mEnd.^2)*sum(sEnd.^2));
            
            ts.tracks(k).corrDisappearance = (sum(bc)>=6) && K>0.8 &&...
                sum(ts.tracks(k).A_binary(mCh,end-4:end))>=4 && sum(ts.tracks(k).A_binary(sCh,end-4:end))>=4;
        end
    end
    
    save([data.source 'Tracking' filesep fileName], '-struct', 'ts');
else
    fprintf('Classification has already been run for %s\n', getShortPath(data));
end
