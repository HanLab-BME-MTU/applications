function [lftRes, res] = runLifetimeAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addOptional('lb', [1  11 16 21 41 61]);
ip.addOptional('ub', [10 15 20 40 60 120]);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off', 'all'})));
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('LifetimeData', 'lifetimeData.mat', @ischar);
ip.addParamValue('Type', 'all', @ischar);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.addParamValue('Print', false, @islogical);
ip.addParamValue('Buffer', 5);
ip.addParamValue('MaxIntensityThreshold', []);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ClassificationSelector', 'significantMaster');
ip.addParamValue('ShowThresholdRange', false, @islogical);
ip.addParamValue('MaxP', 3);
ip.addParamValue('YLim', []);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('RemoveOutliers', true, @islogical);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('FirstNFrames', []);
ip.addParamValue('PoolDatasets', false, @islogical);
%ip.addParamValue('TrackIndex', []);
ip.addParamValue('ShowStatistics', false, @islogical);
ip.addParamValue('SelectIndex', [], @iscell);
ip.parse(data, varargin{:});
lb = ip.Results.lb;
ub = ip.Results.ub;
nc = numel(lb); % # cohorts
mCh = find(strcmp(data(1).source, data(1).channels));
FirstNFrames = ip.Results.FirstNFrames;
selIdx = ip.Results.SelectIndex;

printPath = [getExpDir(data) 'Figures' filesep];
[~,~] = mkdir(printPath);

% median absolute deviation -> standard deviation
% madFactor = 1/norminv(0.75, 0, 1);

% Extend all to max. movie length, in case of mismatch
Nmax = max([data.movieLength])-2;
buffer = ip.Results.Buffer;
cutoff_f = ip.Results.Cutoff_f;

% movieLength = min([data.movieLength]);
framerate = data(1).framerate;

firstN = 3:20;

% loop through data sets, load tracks, store max. intensities and lifetimes
res = struct([]);

[lftData, outlierIdx] = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite,...
    'ReturnValidOnly', false, 'ExcludeVisitors', ip.Results.ExcludeVisitors, 'Cutoff_f', cutoff_f,...
    'Scale', ip.Results.Rescale, 'DisplayScaling', strcmpi(ip.Results.Display, 'on'),...
    'RemoveOutliers', ip.Results.RemoveOutliers,...
    'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData);
if ~isempty(selIdx)
    selIdx(outlierIdx) = [];
end

if ip.Results.PoolDatasets
    pnames = {'lifetime_s', 'start', 'catIdx', 'A', 'lifetime_s_all', 'start_all', 'catIdx_all'};
    if isfield(lftData, 'significantMaster')
        pnames = [pnames 'significantMaster'];
    end
    for k = 1:numel(pnames)
        tmp.(pnames{k}) = vertcat(lftData.(pnames{k}));
    end
    tmp2 = arrayfun(@(i) i.visitors.lifetime_s, lftData, 'unif', 0);
    tmp.visitors.lifetime_s = vertcat(tmp2{:});
    tmp.a = 1;
    lftData = tmp;    
end
fprintf('=================================================\n');
fprintf('Lifetime analysis - processing:   0%%');
nd = numel(lftData);

lftRes.cellArea = zeros(nd,1);
for i = 1:nd
    
    % Category statistics
    idx_Ia = [lftData(i).catIdx_all]==1;
    idx_Ib = [lftData(i).catIdx_all]==2;
    idx_IIa = [lftData(i).catIdx_all]==5;
    
    % raw histograms
    N = data(i).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    lftHist_Ia =  hist(lftData(i).lifetime_s_all(idx_Ia), t);
    lftHist_Ib =  hist(lftData(i).lifetime_s_all(idx_Ib), t);
    lftHist_IIa = hist(lftData(i).lifetime_s_all(idx_IIa), t);
    
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
    res(i).t = (cutoff_f:Nmax)*framerate;
    
    % Normalization
    res(i).lftHist_Ia = lftHist_Ia / sum(lftHist_Ia) / framerate;
    res(i).lftHist_Ib = lftHist_Ib / sum(lftHist_Ib) / framerate;
    res(i).lftHist_IIa = lftHist_IIa / sum(lftHist_IIa) / framerate;
    res(i).nSamples_Ia = sum(idx_Ia);
    
    %-------------------------------------------------------------
    % Max. intensity distribution for cat. Ia CCP tracks
    %-------------------------------------------------------------
    
    for k = 1:nc
        % indexes within cohorts
        cidx = lb(k)<=lftData(i).lifetime_s & lftData(i).lifetime_s<=ub(k);
        res(i).maxA{k} = nanmax(lftData(i).A(cidx,:,mCh),[],2);
        %for n = firstN
        %   res(i).(['maxA_f' num2str(n)]){k} = nanmax(lftData(i).A(cidx,1:n,mCh),[],2);
        %end
        
        % lifetimes for given cohort
        res(i).lft{k} = lftData(i).lifetime_s(cidx);
    end
    
    res(i).maxA_all = nanmax(lftData(i).A(:,:,mCh),[],2);
    if isfield(lftData, 'significantMaster')
       res(i).significantMaster = lftData(i).significantMaster;
    end
    res(i).firstN = firstN;
        
fprintf('\b\b\b\b%3d%%', round(100*i/nd));
end
fprintf('\n');


%====================
% Threshold
%====================
if isempty(ip.Results.MaxIntensityThreshold)
    A = arrayfun(@(i) i.A(:,:,mCh), lftData, 'UniformOutput', false);
    A = vertcat(A{:});
    lft = vertcat(lftData.lifetime_s);
   
    if isempty(FirstNFrames)
        frameRange = 3:12;
        hval = zeros(1,frameRange(end));
        for ni = 1:numel(frameRange)
            M = max(A(:,1:frameRange(ni)),[],2);
            
            muC = zeros(1,nc);
            sC = zeros(1,nc);
            for c = 1:nc
                cidx = lb(c)<=lft & lft<=ub(c);
                [muC(c), sC(c)] = fitGaussianModeToCDF(M(cidx,:));
            end
            hval(frameRange(ni)) = adtest1(muC(2:end), 'mu', muC(1), 'sigma', sC(1)/sqrt(nc));
        end
        FirstNFrames = find(hval==1, 1, 'first')-1;
    end
    
    M = nanmax(A(:,1:FirstNFrames,mCh),[],2);
    
    [mu_g, sigma_g] = fitGaussianModeToCDF(M);
    %[mu_g sigma_g] = fitGaussianModeToPDF(M);
    T = norminv(0.99, mu_g, sigma_g);

    % 95th percentile of first frame intensity distribution
    T95 = prctile(A(:,1,mCh), 95);
    
    fprintf('Max. intensity threshold on first %d frames: %.2f\n', FirstNFrames, T);
    fprintf('95th percentile of 1st frame distribution: %.2f\n', T95);
else
    T = ip.Results.MaxIntensityThreshold;
end

% loop through data sets, apply max. intensity threshold
for i = 1:nd
    
    % Selection indexes for each data set
    % 1) Intensity threshold based on maximum intensity distribution
    idxMI = res(i).maxA_all >= T;
    %idxMI = nanmax(lftData(i).A(:,FirstNFrames:end),[],2)' >= T;
    
    if ~isempty(selIdx)
        idxMI = idxMI & selIdx{i};
    end
    
    % 2) Remove non-endocytic structures (i.e., endosomal CCSs)
    if ip.Results.ExcludeVisitors
        res(i).lftVisitors = lftData(i).visitors.lifetime_s;
        nCS = numel(lftData(i).lifetime_s) + numel(lftData(i).visitors.lifetime_s);
        lftRes.pctVisit(i) = numel(res(i).lftVisitors) / nCS;
    else
        nCS = numel(lftData(i).lifetime_s);
    end
    res(i).lftAboveT = lftData(i).lifetime_s(idxMI);
    res(i).lftBelowT = lftData(i).lifetime_s(~idxMI);
    res(i).nSamplesAboveT = sum(idxMI);
    
    %res(i).maxAAboveT = res(i).maxA_all(idxMI);
    %res(i).AaboveT = lftData(i).A(idxMI,:,:);
    lftRes.pctAbove(i) = numel(res(i).lftAboveT)/nCS;
    lftRes.pctBelow(i) = numel(res(i).lftBelowT)/nCS;
            
    N = data(i).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    %t = (1:N)*framerate;
    lftHist_A = hist(res(i).lftAboveT, t);
    lftHist_B = hist(res(i).lftBelowT, t);

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
    if ip.Results.ExcludeVisitors
        lftHist_V = hist(res(i).lftVisitors, t);
        lftHist_V = [lftHist_V.*w pad0];
        lftRes.lftHist_V(i,:) = lftHist_V / sum(lftHist_V) / framerate;
    end
    
    % Multi-channel data
    if isfield(res, 'significantMaster')
        
        lftHist_pos = hist(lftData(i).lifetime_s(res(i).significantMaster(:,2)==1), t);
        lftHist_neg = hist(lftData(i).lifetime_s(res(i).significantMaster(:,2)==0), t);
        lftHist_pos =  [lftHist_pos.*w  pad0];
        lftHist_neg =  [lftHist_neg.*w  pad0];
        lftRes.lftHist_pos(i,:) = lftHist_pos / sum(lftHist_pos) / framerate;
        lftRes.lftHist_neg(i,:) = lftHist_neg / sum(lftHist_neg) / framerate;
        
        lftHist_Apos = hist(lftData(i).lifetime_s(idxMI & res(i).significantMaster(:,2)), t);
        lftHist_Aneg = hist(lftData(i).lifetime_s(idxMI & ~res(i).significantMaster(:,2)), t);
        lftHist_Bpos = hist(lftData(i).lifetime_s(~idxMI & res(i).significantMaster(:,2)), t);
        lftHist_Bneg = hist(lftData(i).lifetime_s(~idxMI & ~res(i).significantMaster(:,2)), t);
        lftHist_Apos =  [lftHist_Apos.*w  pad0];
        lftHist_Aneg =  [lftHist_Aneg.*w  pad0];
        lftHist_Bpos =  [lftHist_Bpos.*w  pad0];
        lftHist_Bneg =  [lftHist_Bneg.*w  pad0];
        lftRes.lftHist_Apos(i,:) = lftHist_Apos / sum(lftHist_Apos) / framerate;
        lftRes.lftHist_Aneg(i,:) = lftHist_Aneg / sum(lftHist_Aneg) / framerate;
        lftRes.lftHist_Bpos(i,:) = lftHist_Bpos / sum(lftHist_Bpos) / framerate;
        lftRes.lftHist_Bneg(i,:) = lftHist_Bneg / sum(lftHist_Bneg) / framerate;
        
        %lftRes.lftHist_Apos(i,:) = lftHist_Apos / normA;
        %lftRes.lftHist_Aneg(i,:) = lftHist_Aneg / normA;
        %lftRes.lftHist_Bpos(i,:) = lftHist_Bpos / normB;
        %lftRes.lftHist_Bneg(i,:) = lftHist_Bneg / normB;
        
        lftRes.pctAboveSignificant(i) =    sum(idxMI & res(i).significantMaster(:,2))/numel(idxMI);
        lftRes.pctAboveNotSignificant(i) = sum(idxMI & ~res(i).significantMaster(:,2))/numel(idxMI);
        lftRes.pctBelowSignificant(i) =    sum(~idxMI & res(i).significantMaster(:,2))/numel(idxMI);
        lftRes.pctBelowNotSignificant(i) = sum(~idxMI & ~res(i).significantMaster(:,2))/numel(idxMI);
    end
    %-----------------------------------
    % Initiation density
    %-----------------------------------
    % birth/death statistics
    startsPerFrame_all = hist(lftData(i).start_all, 1:data(i).movieLength);
    startsPerFrame_all = startsPerFrame_all(6:end-2);
    startsPerFrame_Ia = hist(lftData(i).start(lftData(i).catIdx==1), 1:data(i).movieLength);
    startsPerFrame_Ia = startsPerFrame_Ia(6:end-2);
    startsPerFrameAbove = hist(lftData(i).start(idxMI), 1:data(i).movieLength);
    startsPerFrameAbove = startsPerFrameAbove(6:end-2);
    
    % Cell area
    px = data(i).pixelSize / data(i).M; % pixels size in object space
    mpath = [data(i).source 'Detection' filesep 'cellmask.tif'];
    mask = logical(imread(mpath));
    lftRes.cellArea(i) = sum(mask(:)) * px^2 / 1e-12; % in µm^2
    
    % in µm^-2 min^-1
    %lftRes.initDensity_all(i,:) = [median(startsPerFrame_all); madFactor*mad(startsPerFrame_all, 1)]/data(i).framerate*60/lftRes.cellArea(i);
    %lftRes.initDensity_Ia(i,:) = [median(startsPerFrame_Ia); madFactor*mad(startsPerFrame_Ia, 1)]/data(i).framerate*60/lftRes.cellArea(i);
    %lftRes.initDensity_above(i,:) = [median(startsPerFrameAbove); madFactor*mad(startsPerFrameAbove,1)]/data(i).framerate*60/lftRes.cellArea(i);
    lftRes.initDensity_all(i,:) = [mean(startsPerFrame_all); std(startsPerFrame_all)]/data(i).framerate*60/lftRes.cellArea(i);
    lftRes.initDensity_Ia(i,:) = [mean(startsPerFrame_Ia); std(startsPerFrame_Ia)]/data(i).framerate*60/lftRes.cellArea(i);
    lftRes.initDensity_above(i,:) = [mean(startsPerFrameAbove); std(startsPerFrameAbove)]/data(i).framerate*60/lftRes.cellArea(i);
end
%====================
% Initiation density
%====================
fprintf(2, 'Initiation density, all tracks:\n');
% for i = 1:nd
%     fprintf('%s: %.3f ± %.3f [µm^-2 min^-1]\n', getCellDir(data(i)), lftRes.initDensity_all(i,1), lftRes.initDensity_all(i,2));
% end
% fprintf('Initiation density, SEM: %f ± %f [µm^-2 min^-1]\n', mean(D(1,:)), std(D(1,:))/sqrt(nd));
fprintf(2, 'Average: %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensity_all(:,1)), std(lftRes.initDensity_all(:,1)));
fprintf('-------------------------------------------------\n');
fprintf(2, 'Initiation density, valid tracks only:\n');
% for i = 1:nd
%     fprintf('%s: %.3f ± %.3f [µm^-2 min^-1]\n', getCellDir(data(i)), lftRes.initDensity_Ia(i,1), lftRes.initDensity_Ia(i,2));
% end
fprintf(2, 'Average: %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensity_Ia(:,1)), std(lftRes.initDensity_Ia(:,1)));
fprintf('-------------------------------------------------\n');


lftRes.t = (cutoff_f:Nmax)*framerate;
lftRes.meanLftHist_A = nanmean(lftRes.lftHist_A,1);
lftRes.meanLftHist_B = nanmean(lftRes.lftHist_B,1);
% lftRes.meanLftHist_Ia = mean(vertcat(res.lftHist_Ia),1);
lftRes.lftHist_Ia = vertcat(res.lftHist_Ia);
lftRes.nSamples_Ia = [res.nSamples_Ia];
lftRes.nSamplesAboveT = [res.nSamplesAboveT];
lftRes.data = data;

if ip.Results.ExcludeVisitors
    lftRes.meanLftHist_V = mean(lftRes.lftHist_V(i,:),1);
end

pctV = [50 25 75 5 95]';
pAll = arrayfun(@(i) prctile(i.lifetime_s,pctV), lftData, 'UniformOutput', false);
pA = arrayfun(@(i) prctile(i.lftAboveT,pctV), res, 'UniformOutput', false);
pB = arrayfun(@(i) prctile(i.lftBelowT,pctV), res, 'UniformOutput', false);
lftRes.stats = [mean([pAll{:}],2) mean([pB{:}],2) mean([pA{:}],2)];
%lftRes.stats = [prctile([lftData.lifetime_s],pctV) prctile([res.lftBelowT],pctV) prctile([res.lftAboveT],pctV)]

% stats = [prctile([lftData.lifetime_s], [50 25 75])' prctile([res.lftBelowT], [50 25 75])' prctile([res.lftAboveT], [50 25 75])' ];
% w = 1.5*(stats(3,:)-stats(2,:));
% lftRes.stats = [stats; stats(2,:)-w; stats(3,:)+w];


%---------------------------------------
% Raw lifetime distributions + average
%---------------------------------------
if strcmpi(ip.Results.Display, 'all')
    a = [lftData.a];
    colorV = hsv(nd);
    [~,idxa] = sort(a);
    [~,idxa] = sort(idxa);
    colorV = colorV(idxa,:);
    fset = loadFigureSettings('print');
    
    figure(fset.fOpts{:}, 'Name', 'Raw lifetime distribution');
    axes(fset.axOpts{:});
    hold on;
    for i = nd:-1:1
        plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
    end
    plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
    ya = 0:0.02:0.1;
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    % Inset with zoom
    zf = 0.6;
    aw = fset.axPos(3);
    ah = fset.axPos(4);
    axes(fset.axOpts{:}, 'Units', fset.units, 'Position', [fset.axPos(1)+(1-zf)*aw fset.axPos(2)+(1-zf)*ah zf*aw zf*ah]);
    hold on;
    for i = nd:-1:1
        plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
    end
    plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
    axis([0 60 0 0.035]);
    ya = 0:0.01:0.04;
    set(gca, 'FontSize', 7, 'TickLength', fset.TickLength/zf, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
    
    % print('-depsc2', '-loose', ['LftRaw_dataOX_10_cut' num2str(cutoff_f) '.eps']);
    
    lftCDF = cumsum(mean(vertcat(lftRes.lftHist_Ia),1))*framerate;
    [uCDF, idx] = unique(lftCDF);
    lft50 = interp1(uCDF, lftRes.t(idx), 0.5);
    
    
    figure(fset.fOpts{:}, 'Name', 'Raw lifetime distribution');
    axes(fset.axOpts{:});
    hold on;
    meanHist = mean(vertcat(lftRes.lftHist_Ia), 1);
    plot(lftRes.t, meanHist, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    ya = 0:0.02:0.1;
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, fset.axOpts{:}, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    [mu,~,Aexp] = fitExpToHist(lftRes.t, meanHist);
    % plot(lftRes.t, expF, 'r--', 'LineWidth', 1);
    ti = 0:0.1:120;
    %plot(ti, Aexp/mu*exp(-1/mu*ti), 'r--', 'LineWidth', 1);
    plot(lftRes.t, Aexp/mu*exp(-1/mu*lftRes.t), 'r-', 'LineWidth', 1);
    %hl = legend(' Best exponential fit', 'Location', 'SouthEast');
    %set(hl, 'Box', 'off', 'Position', [4.5 1.5 2 1]);
    
    % Inset with zoom
    zf = 0.6;
    aw = fset.axPos(3);
    ah = fset.axPos(4);
    axes(fset.axOpts{:}, 'Units', fset.units, 'Position', [fset.axPos(1)+(1-zf)*aw fset.axPos(2)+(1-zf)*ah zf*aw zf*ah]);
    hold on;
    idx = find(lftRes.t==round(lft50/framerate)*framerate);
    fill([lftRes.t(1:idx) lftRes.t(idx:-1:1)], [lftCDF(1:idx) zeros(1,idx)], fset.ceB, 'EdgeColor', 'none');
    plot(lftRes.t, lftCDF, 'k', 'LineWidth', 1.5);
    plot([0 lft50], [0.5 0.5], 'k--', 'LineWidth', 1);
    ya = 0:0.25:1;
    set(gca, 'FontSize', 7, 'TickLength', fset.TickLength/zf, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    % xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Cumulative freq.', fset.sfont{:});
    
    % print('-depsc2', '-loose', ['LftMean+CDF_dataOX_10_cut' num2str(cutoff_f) '_mu=' num2str(mu, '%.2f') '.eps']);
end

%---------------------------------------
% CDF plot of the raw lifetimes
%---------------------------------------

% figure(fset.fOpts{:}, 'Name', 'Cumulative lifetime distribution');
% axes(fset.axOpts{:});
% hold on;
% idx = find(lftRes.t==round(lft50/framerate)*framerate);
% fill([lftRes.t(1:idx) lftRes.t(idx:-1:1)], [lftCDF(1:idx) zeros(1,idx)], fset.ceB, 'EdgeColor', 'none');
% plot(lftRes.t, lftCDF, 'k', 'LineWidth', 2);
% plot([0 lft50], [0.5 0.5], 'k--');
% 
% ya = 0:0.25:1;
% axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
% set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
% xlabel('Lifetime (s)', fset.lfont{:});
% ylabel('Cumulative frequency', fset.lfont{:});
% % print('-depsc2', '-loose', ['LftRawCDF_dataOX_10_cut' num2str(cutoff_f) '.eps']);


if strcmpi(ip.Results.Display, 'on')
    h = plotLifetimes(lftRes, 'ShowStatistics', ip.Results.ShowStatistics);
    print(h(1), '-depsc2', '-loose', [printPath 'lifetimeDistributions.eps']);
    if ip.Results.ShowStatistics
        print(h(2), '-depsc2', '-loose', [printPath 'lifetimeDistributionsStats.eps']);
    end
end
return


%=====================================================================
% Display histogram for range of thresholds
%=====================================================================
if ip.Results.ShowThresholdRange
    Trange = 40:10:200;
    for ti = 1:numel(Trange)
        T = Trange(ti);
        
        for i = 1:nd
            idx = res(i).maxA_all >= T;
            lftAboveT = lftData(i).lifetime_s(idx);
            lftBelowT = lftData(i).lifetime_s(~idx);
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
            tmp.lftHist_A(i,:) = lftHist_A / sum(lftHist_A);
            tmp.lftHist_B(i,:) = lftHist_B / sum(lftHist_B);
            
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
