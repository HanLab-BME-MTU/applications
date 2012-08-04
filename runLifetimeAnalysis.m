function [lftRes fitRes] = runLifetimeAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addOptional('lb', [3:10 11 16 21 41 61 81 101 141]);
ip.addOptional('ub', [3:10 15 20 40 60 80 100 140 200]);
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('FileName', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('Type', 'all', @ischar);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.addParamValue('Print', false, @islogical);
ip.addParamValue('Buffer', 5);
ip.addParamValue('MaxIntensityThreshold', []);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ClassificationSelector', 'significantSignal');
ip.addParamValue('ShowThresholdRange', false, @islogical);
ip.addParamValue('MaxP', 3);
ip.addParamValue('YLim', []);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.parse(data, varargin{:});
lb = ip.Results.lb;
ub = ip.Results.ub;
nd = length(data);
nc = numel(lb); % # cohorts
mCh = find(strcmp(data(1).source, data(1).channels));

% median absolute deviation -> standard deviation
madFactor = 1/norminv(0.75, 0, 1);

% Extend all to max. movie length, in case of mismatch
Nmax = max([data.movieLength])-2;
buffer = ip.Results.Buffer;
cutoff_f = ip.Results.Cutoff_f;

framerate = data(1).framerate;

firstN = 3:20;

lftData = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite, 'FileName', ip.Results.FileName);

% loop through data sets, load tracks, store max. intensities and lifetimes
res = struct([]);
% lftRes = struct([]);

%==============================================================
% Outlier detection (based on max. intensity distribution)
%==============================================================
maxA_all = arrayfun(@(i) nanmax(i.A(:,:,mCh),[],2)', lftData, 'UniformOutput', false);

% Rescale EDFs (correction for FP-fusion expression level)
[a offset refIdx] = rescaleEDFs(maxA_all, 'Display', true);

% apply intensity scaling
% for i = 1:nd
%     maxA_all{i} = a(i) * maxA_all{i};
%     lftData(i).A(:,:,mCh) = a(i) * lftData(i).A(:,:,mCh);
%     lftData(i).sbA(:,:,mCh) = a(i) * lftData(i).sbA(:,:,mCh);
%     lftData(i).ebA(:,:,mCh) = a(i) * lftData(i).ebA(:,:,mCh);
%     lftData(i).sigma_r(:,:,mCh) = a(i) * lftData(i).sigma_r(:,:,mCh);
% end

% outlierIdx = detectEDFOutliers(maxA_all, offset, refIdx);
% fprintf('Outlier data sets:\n');
% for i = 1:numel(outlierIdx)
%     fprintf('%s\n', getShortPath(data(outlierIdx(i))));
% end
% 
% lftData(outlierIdx) = [];
% data(outlierIdx) = [];
% nd = numel(data);
% clear a outlierIdx maxA_all;

% apply lifetime scaling (based on intensity scaling)
% reference: lowest-intensity scale
relativeScale = min(a)./a;
alpha = 0.05;
kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background

% for each dataset, loop through tracks, divide amplitude by relative scale, calc. new lifetime
% and #frames lost at beginning and end => statistics
for i = 1:nd
    % loop through 'Ia' tracks only
    aMat = relativeScale(i) * lftData(i).A(:,:,mCh);
    aStdMat = relativeScale(i) * lftData(i).A_pstd(:,:,mCh);
    % background level remains the same!
    sMat = kLevel * lftData(i).sigma_r(:,:,mCh);
    sStdMat = kLevel * lftData(i).SE_sigma_r(:,:,mCh);
    lifetime_s = lftData(i).lifetime_s(lftData(i).catIdx==1);
    trackLengths = lftData(i).trackLengths(lftData(i).catIdx==1);
    nt = numel(trackLengths);
    deltaS = NaN(1,nt);
    deltaE = NaN(1,nt);
    lifetimeScaled = NaN(1,nt);
    for k = 1:nt
        A = aMat(k,1:trackLengths(k));
        sigma_r = sMat(k,1:trackLengths(k));
        A_pstd = aStdMat(k,1:trackLengths(k));
        SE_sigma_r = sStdMat(k,1:trackLengths(k));
        npx = round((sigma_r./SE_sigma_r).^2/2+1);
        df2 = (npx-1) .* (A_pstd.^2 + SE_sigma_r.^2).^2 ./ (A_pstd.^4 + SE_sigma_r.^4);
        scomb = sqrt((A_pstd.^2 + SE_sigma_r.^2)./npx);
        T = (A - sigma_r) ./ scomb;
        pval_Ar_scaled = tcdf(-T, df2);
        %pval_Ar_scaled
        %lftData(i).pvalMat(k,1:trackLengths(k))
        %hvalRef = lftData(i).pvalMat(k,1:trackLengths(k))<0.05;

        % binary mask of the track
        hval_Ar_scaled = pval_Ar_scaled<0.05;
        if sum(hval_Ar_scaled)~=0
            deltaS(k) = find(hval_Ar_scaled==1, 1, 'first')-1;
            deltaE(k) = trackLengths(k) - find(hval_Ar_scaled==1, 1, 'last');
            lifetimeScaled(k) = lifetime_s(k) - (deltaS(k)+deltaE(k))*framerate;
        end
        % rough estimate: first and last points detected are limits
    end
    lftData(i).lifetimeScaled = lifetimeScaled;
    lftData(i).deltaS = deltaS;
    lftData(i).deltaE = deltaE;
end


lftFields = {'lifetime_s', 'trackLengths', 'start', 'catIdx'};
fprintf('=================================================\n');
fprintf('Lifetime analysis - processing:   0%%');
lftRes.cellArea = zeros(nd,1);
for i = 1:nd
    
    % apply frames cutoff for short tracks
    lftData(i).A(lftData(i).trackLengths(lftData(i).catIdx==1)<cutoff_f,:) = [];
    idx = lftData(i).trackLengths(lftData(i).catIdx==1) < cutoff_f;
    for f = 1:numel(lftFields)
        lftData(i).(lftFields{f}) = lftData(i).(lftFields{f})(lftData(i).catIdx==1);
        lftData(i).(lftFields{f})(idx) = [];
    end
    lifetime_s = lftData(i).lifetime_s;
    lftData(i).lifetimeScaled = lftData(i).lifetimeScaled(lftData(i).lifetimeScaled>=cutoff_f);
    
    % Category statistics
    idx_Ia = [lftData(i).catIdx]==1;
    idx_Ib = [lftData(i).catIdx]==2;
    idx_IIa = [lftData(i).catIdx]==5;
    v = hist([lftData(i).catIdx], 1:8);
    v = v/numel(lifetime_s);
    lftRes.trackClassStats(i,:) = v;
    
    % raw histograms
    N = data(i).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    lftHist_Ia = hist(lifetime_s(idx_Ia), t);
    lftHist_Ib = hist(lifetime_s(idx_Ib), t);
    lftHist_IIa = hist(lifetime_s(idx_IIa), t);
    lftHist_scaled = hist(lftData(i).lifetimeScaled, t);
    
    % apply correction
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    w = N./(N-cutoff_f+1:-1:1);
    pad0 = zeros(1,Nmax-N);
    lftHist_Ia =  [lftHist_Ia.*w  pad0];
    lftHist_Ib =  [lftHist_Ib.*w  pad0];
    lftHist_IIa = [lftHist_IIa.*w pad0];
    lftHist_scaled = [lftHist_scaled.*w pad0];
    res(i).t = (cutoff_f:Nmax)*framerate;
    
    % Normalization
    res(i).lftHist_Ia = lftHist_Ia / sum(lftHist_Ia) / framerate;
    res(i).lftHist_Ib = lftHist_Ib / sum(lftHist_Ib) / framerate;
    res(i).lftHist_IIa = lftHist_IIa / sum(lftHist_IIa) / framerate;
    res(i).lftHist_scaled = lftHist_scaled / sum(lftHist_scaled) / framerate;
    res(i).nSamples_Ia = sum(idx_Ia);
    
    % birth/death statistics
    startsPerFrame_all = hist(lftData(i).start, 1:data(i).movieLength);
    startsPerFrame_all = startsPerFrame_all(6:end-2);
    startsPerFrame_Ia = hist(lftData(i).start(idx_Ia), 1:data(i).movieLength);
    startsPerFrame_Ia = startsPerFrame_Ia(6:end-2);
    
    %-------------------------------------------------------------
    % Initiation density
    %-------------------------------------------------------------
    
    % Cell area
    px = data(i).pixelSize / data(i).M; % pixels size in object space
    mpath = [data(i).source 'Detection' filesep 'cellmask.tif'];
    mask = logical(imread(mpath));
    lftRes.cellArea(i) = sum(mask(:)) * px^2 / 1e-12; % in µm^2
    
    % in µm^-2 min^-1
    lftRes.initDensity_all(i,:) = [median(startsPerFrame_all); madFactor*mad(startsPerFrame_all, 1)]/data(i).framerate*60/lftRes.cellArea(i);
    lftRes.initDensity_Ia(i,:) = [median(startsPerFrame_Ia); madFactor*mad(startsPerFrame_Ia, 1)]/data(i).framerate*60/lftRes.cellArea(i);
    
    %-------------------------------------------------------------
    % Max. intensity distribution for cat. Ia CCP tracks
    %-------------------------------------------------------------
    lifetime_s = lifetime_s(idx_Ia);
    
    for k = 1:nc
        % indexes within cohorts
        cidx = lb(k)<=lifetime_s & lifetime_s<=ub(k);
        res(i).maxA{k} = nanmax(lftData(i).A(cidx,:),[],2);
        for n = firstN
           res(i).(['maxA_f' num2str(n)]){k} = nanmax(lftData(i).A(cidx,1:n),[],2);
        end
        
        % lifetimes for given cohort
        res(i).lft{k} = lifetime_s(cidx);
    end
    
    res(i).lft_all = lifetime_s;
    res(i).maxA_all = nanmax(lftData(i).A,[],2)';
    if isfield(lftData, 'significantSignal')
        res(i).significantSignal = lftData(i).significantSignal(:,idx_Ia);
    end
    res(i).firstN = firstN;
    
fprintf('\b\b\b\b%3d%%', round(100*i/nd));
end
fprintf('\n');

% Plot raw and scaled lifetime histograms
% figure;
% hold on;
% for i = 1:nd
%     plot(res(i).t, res(i).lftHist_Ia, 'k');    
%     plot(res(i).t, res(i).lftHist_scaled, 'r--');
% end

% tmp = mean(vertcat(res.lftHist_scaled),1);
% medIdx = find(a==1);
% D1 = hist(lftData(medIdx).deltaS, 0:25);
% D1 = D1/sum(D1);
% D2 = hist(lftData(medIdx).deltaE, 0:25);
% D2 = D2/sum(D2);
% D = conv(D1,D2(end:-1:1));
% D = D(numel(D2):end);
% D = D/sum(D);
% 
% figure;
% hold on;
% % plot(res(i).t, mean(vertcat(res.lftHist_Ia),1), 'k');
% % plot(res(i).t, tmp, 'r--');
% plot(mean(vertcat(res.lftHist_Ia),1), 'k');
% plot(tmp, 'r--');
% tt = conv(tmp, D(end:-1:1));
% tt = tt(numel(D):end);
% plot(tt, 'c');

%====================
% Initiation density
%====================
fprintf(2, 'Initiation density, all tracks:\n');
for i = 1:nd
    fprintf('%s: %.3f ± %.3f [µm^-2 min^-1]\n', getCellDir(data(i)), lftRes.initDensity_all(i,1), lftRes.initDensity_all(i,2));
end
% fprintf('Initiation density, SEM: %f ± %f [µm^-2 min^-1]\n', mean(D(1,:)), std(D(1,:))/sqrt(nd));
fprintf(2, 'Average: %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensity_all(:,1)), std(lftRes.initDensity_all(:,1)));
fprintf('-------------------------------------------------\n');
fprintf(2, 'Initiation density, valid tracks only:\n');
for i = 1:nd
    fprintf('%s: %.3f ± %.3f [µm^-2 min^-1]\n', getCellDir(data(i)), lftRes.initDensity_Ia(i,1), lftRes.initDensity_Ia(i,2));
end
fprintf(2, 'Average: %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensity_Ia(:,1)), std(lftRes.initDensity_Ia(:,1)));
fprintf('-------------------------------------------------\n');



%====================
% Threshold
%====================
if isempty(ip.Results.MaxIntensityThreshold)
    % lifetime cohort: [5..10] seconds
    % combine first 5 frames from all cohorts
    
    maxIntDistCat_f5 = horzcat(res.maxA_f5);
    maxIntDistCat_f5 = vertcat(maxIntDistCat_f5{:});
    
    %[mu_g sigma_g] = fitGaussianModeToPDF(maxIntDistCat_f5);
    [mu_g sigma_g] = fitGaussianModeToCDF(maxIntDistCat_f5);
    T = norminv(0.99, mu_g, sigma_g);
    fprintf('Max. intensity threshold: %.2f\n', T);
else
    T = ip.Results.MaxIntensityThreshold;
end

% intMat_Ia_all = vertcat(lftData.intMat_Ia);
minLength = min([data.movieLength]);
intMat_Ia_all = arrayfun(@(i) i.A(:,1:minLength), lftData, 'UniformOutput', false);
intMat_Ia_all = vertcat(intMat_Ia_all{:});
lifetime_s_all = [res.lft_all];

tx = 30;

% 95th percentile of the reference (above threshold) distribution
pRef = prctile(intMat_Ia_all(lifetime_s_all>=tx,2:4), 95, 1);

% loop through data sets, apply max. intensity threshold
for i = 1:nd
    
    % Selection indexes for each data set
    % 1) Intensity threshold based on maximum intensity distribution
    idxMI = res(i).maxA_all >= T;
    
    % 2) Lifetime threshold for objects with a faster-than-tolerated* growth rate
    if ip.Results.ExcludeVisitors
        idxLft = sum(lftData(i).A(:,2:4)>repmat(pRef, [size(lftData(i).A,1) 1]),2)>0 & res(i).lft_all'<tx;
        % combined index
        idxMI = idxMI & ~idxLft';
    end
    
    res(i).lftAboveT = res(i).lft_all(idxMI);
    res(i).lftBelowT = res(i).lft_all(~idxMI);
    lftRes.pctAbove(i) = sum(idxMI)/numel(idxMI);
    
    N = data(i).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    %t = (1:N)*framerate;
    lftHist_A = hist(res(i).lftAboveT, t);
    lftHist_B = hist(res(i).lftBelowT, t);
    
    %lftHist_A(1:cutoff_f-1) = [];
    %lftHist_B(1:cutoff_f-1) = [];
    
    % apply correction
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    w = N./(N-cutoff_f+1:-1:1);
    pad0 = zeros(1,Nmax-N);
    lftHist_A = [lftHist_A.*w pad0];
    lftHist_B = [lftHist_B.*w pad0];
    
    % Normalization
    %normA = sum(lftHist_A);
    %normB = sum(lftHist_B);
    lftRes.lftHist_A(i,:) = lftHist_A / sum(lftHist_A) / framerate;
    lftRes.lftHist_B(i,:) = lftHist_B / sum(lftHist_B) / framerate;
    
    % Multi-channel data
    if isfield(res, 'significantSignal')
        lftHist_Apos = hist(res(i).lft_all(idxMI & res(i).significantSignal(2,:)), t);
        lftHist_Aneg = hist(res(i).lft_all(idxMI & ~res(i).significantSignal(2,:)), t);
        lftHist_Bpos = hist(res(i).lft_all(~idxMI & res(i).significantSignal(2,:)), t);
        lftHist_Bneg = hist(res(i).lft_all(~idxMI & ~res(i).significantSignal(2,:)), t);
        lftHist_Apos =  [lftHist_Apos.*w  pad0];
        lftHist_Aneg =  [lftHist_Aneg.*w  pad0];
        lftHist_Bpos =  [lftHist_Bpos.*w  pad0];
        lftHist_Bneg =  [lftHist_Bneg.*w  pad0];
        lftRes.lftHist_Apos(i,:) = lftHist_Apos / sum(lftHist_Apos) / framerate;
        lftRes.lftHist_Aneg(i,:) = lftHist_Aneg / sum(lftHist_Aneg) / framerate;
        lftRes.lftHist_Bpos(i,:) = lftHist_Bpos / sum(lftHist_Bpos) / framerate;
        lftRes.lftHist_Bneg(i,:) = lftHist_Bneg / sum(lftHist_Bneg) / framerate;
        %lftRes.lftHist_Apos(i,:) = lftHist_Apos / normA / framerate;
        %lftRes.lftHist_Aneg(i,:) = lftHist_Aneg / normA / framerate;
        %lftRes.lftHist_Bpos(i,:) = lftHist_Bpos / normB / framerate;
        %lftRes.lftHist_Bneg(i,:) = lftHist_Bneg / normB / framerate;
        
        lftRes.pctAboveSignificant(i) = sum(idxMI & res(i).significantSignal(2,:))/numel(idxMI);
        lftRes.pctAboveNS(i) = sum(idxMI & ~res(i).significantSignal(2,:))/numel(idxMI);
        lftRes.pctBelowSignificant(i) = sum(~idxMI & res(i).significantSignal(2,:))/numel(idxMI);
    end
end

lftRes.t = (cutoff_f:Nmax)*framerate;
lftRes.meanLftHist_A = mean(lftRes.lftHist_A,1);
lftRes.meanLftHist_B = mean(lftRes.lftHist_B,1);
% lftRes.meanLftHist_Ia = mean(vertcat(res.lftHist_Ia),1);
lftRes.lftHist_Ia = vertcat(res.lftHist_Ia);
lftRes.nSamples_Ia = [res.nSamples_Ia];
lftRes.data = data;




%%
%---------------------------------------
% Raw lifetime distributions + average
%---------------------------------------
colorV = hsv(nd);
[~,idxa] = sort(a);
% [~,idxa] = sort(idxa);
colorV = colorV(idxa,:);
fset = loadFigureSettings('print');

figure;
axes(fset.axOpts{:});
hold on;
for i = 1:nd
    plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
end
plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
ya = 0:0.02:0.1;
axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
set(gca, fset.axOpts{:}, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});

% Inset with zoom
zf = 0.6;
aw = fset.axPos(3);
ah = fset.axPos(4);
axes(fset.axOpts{:}, 'Units', fset.units, 'Position', [fset.axPos(1)+(1-zf)*aw fset.axPos(2)+(1-zf)*ah zf*aw zf*ah]);
hold on;
for i = 1:nd
    plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
end
hp = plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
axis([0 60 0 0.035]);
ya = 0:0.01:0.04;
set(gca, 'TickLength', fset.TickLength/zf, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);

% print('-depsc2', ['LftRaw_dataOX_10_cut' num2str(cutoff_f) '.eps']);

%%
%---------------------------------------
% CDF plot of the raw lifetimes
%---------------------------------------
lftCDF = cumsum(mean(vertcat(lftRes.lftHist_Ia),1));
[uCDF idx] = unique(lftCDF);
lft50 = interp1(uCDF, lftRes.t(idx), 0.5);

figure;
axes(fset.axOpts{:});
hold on;
idx = find(lftRes.t==round(lft50/framerate));
fill([lftRes.t(1:idx) lftRes.t(idx:-1:1)], [lftCDF(1:idx) zeros(1,idx)], fset.ceB, 'EdgeColor', 'none');
plot(lftRes.t, cumsum(mean(vertcat(lftRes.lftHist_Ia), 1)), 'k', 'LineWidth', 2);

ya = 0:0.2:1;
axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
set(gca, fset.axOpts{:}, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Cumulative frequency', fset.lfont{:});
% print('-depsc2', ['LftRawCDF_dataOX_10_cut' num2str(cutoff_f) '.eps']);

%%


plotLifetimes(lftRes);


%=====================================================================
% Display histogram for range of thresholds
%=====================================================================
if ip.Results.ShowThresholdRange
    Trange = 40:10:200;
    for ti = 1:numel(Trange)
        T = Trange(ti);
        
        for i = 1:nd
            idx = res(i).maxA_all >= T;
            lftAboveT = res(i).lft_all(idx);
            lftBelowT = res(i).lft_all(~idx);
            tmp.pctAbove(i) = sum(idx)/numel(idx);
            
            N = data(i).movieLength-2*buffer;
            t = (cutoff_f:N)*framerate;
            lftHist_A = hist(lftAboveT, t);
            lftHist_B = hist(lftBelowT, t);
            
            % apply correction
            w = N./(N-cutoff_f+1:-1:1);
            pad0 = zeros(1,Nmax-N);
            lftHist_A =  [lftHist_A.*w  pad0];
            lftHist_B =  [lftHist_B.*w  pad0];
            
            % Normalization
            tmp.lftHist_A(i,:) = lftHist_A / sum(lftHist_A) / framerate;
            tmp.lftHist_B(i,:) = lftHist_B / sum(lftHist_B) / framerate;
            
        end
        tComp(ti).meanLftHist_A = mean(tmp.lftHist_A,1);
        tComp(ti).meanLftHist_B = mean(tmp.lftHist_B,1);
        tComp(ti).t_hist = (cutoff_f:Nmax)*framerate;
        tComp(ti).pctAbove = mean(tmp.pctAbove);
    end

    opts = {'.-', 'LineWidth', 2, 'MarkerSize', 16};
    fset = loadFigureSettings();
    cmap = jet(numel(Trange));
    cv = rgb2hsv(cmap);
    cv(:,2) = 0.2;
    cv = hsv2rgb(cv);
    
    figure;
    hold on;
    for ti = 1:numel(Trange)
        
        plot(tComp(ti).t_hist, tComp(ti).meanLftHist_B, opts{:}, 'Color', cv(ti,:));
        hp(ti) = plot(tComp(ti).t_hist, tComp(ti).meanLftHist_A, opts{:}, 'Color', cmap(ti,:));
        
        %legendText = {['Above threshold (' num2str(mean(lftRes.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100,'%.1f') ' %)'],...
        %    ['Below threshold (' num2str(mean(1-lftRes.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100,'%.1f') ' %)']};
    end
    axis([0 min(120, lftRes.t(end)) 0 0.05]);
    set(gca, 'LineWidth', 2, fset.sfont{:}, fset.axOpts{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    hl = legend(hp, arrayfun(@(x) num2str(x, '%.2f'), [tComp.pctAbove], 'UniformOutput', false), 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.ifont{:});
    
end






%=====================================================================
% Fit lifetime histogram with Weibull-distributed populations
%=====================================================================
% fitResCDF = fitLifetimeDistWeibullModel(lftRes, 'Mode', 'CDF');
% plotLifetimeDistModel(lftRes, fitResCDF);

fitResPDF = fitLifetimeDistWeibullModel(lftRes, 'Mode', 'PDF', 'MaxP', ip.Results.MaxP);
plotLifetimeDistModel(lftRes, fitResPDF, 'YLim', ip.Results.YLim);

fitRes = fitResPDF;

% fitResCDF = fitLifetimeDistGammaModel(lftRes, 'Mode', 'CDF');
% plotLifetimeDistModel(lftRes, fitResCDF);

% fitResPDF = fitLifetimeDistGammaModel(lftRes, 'Mode', 'PDF', 'MaxP', ip.Results.MaxP);
% plotLifetimeDistModel(lftRes, fitResPDF);



return
  
    %====================
    % Gap statistics
    %==================== 
%     binEdges = [0:20:120 data(k).movieLength-data(k).framerate];
%     nb = length(binEdges)-1;
%     gapsPerTrack_Ia = zeros(1,nb);
%     gapsPerTrack_Ib = zeros(1,nb);
%     gapsPerTrack_IIa = zeros(1,nb);
%     for b = 1:nb
%         tidx = binEdges(b)<=lifetimes_s & lifetimes_s<binEdges(b+1);
%         gapsPerTrack_Ia(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_Ia & tidx)));
%         gapsPerTrack_Ib(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_Ib & tidx)));
%         gapsPerTrack_IIa(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_IIa & tidx)));
%     end
%     res.gapsPerTrack_Ia{k} = gapsPerTrack_Ia;
%     res.gapsPerTrack_Ib{k} = gapsPerTrack_Ib;
%     res.gapsPerTrack_IIa{k} = gapsPerTrack_IIa;


if strcmpi(ip.Results.Display, 'on')
    
    % gap statistics
%     ce = fset.ceTrackClasses([1 2 5],:);
%     cf = fset.cfTrackClasses([1 2 5],:);
%     xlabels = arrayfun(@(b) [num2str(binEdges(b)) '-' num2str(binEdges(b+1)) ' s'], 1:numel(binEdges)-1, 'UniformOutput', false);
    
    
%     M_Ia = vertcat(res.gapsPerTrack_Ia{:});
%     M_Ib = vertcat(res.gapsPerTrack_Ib{:});
%     M_IIa = vertcat(res.gapsPerTrack_IIa{:});
%     
%     M = [mean(M_Ia,1); mean(M_Ib,1); mean(M_IIa,1)]';
%     S = [std(M_Ia,[],1); std(M_Ib,[],1); std(M_IIa,[],1)]';
%     
%     hf(3) = figure; barplot2(M, S, 'FaceColor', cf, 'EdgeColor', ce,...
%         'XLabels', xlabels, 'XLabel', 'Lifetime cohort', 'YLabel', 'gaps/track');
    
end

if ip.Results.Print
    
    fpath = cell(1,nd);
    for k = 1:nd
        [~,fpath{k}] = getCellDir(data(k));
    end
    fpath = unique(fpath);
    if numel(fpath)>1
        fprintf('Figures could not be printed.');
    else
        fpath = [fpath{1} 'Figures' filesep];
        [~,~] = mkdir(fpath);
        print(hf(1), '-depsc2', [fpath 'trackClassDistribution.eps']);
        print(hf(2), '-depsc2', [fpath 'meanLftHist_classes.eps']);
        print(hf(3), '-depsc2', [fpath 'gapStatistics.eps']);
    end
end
