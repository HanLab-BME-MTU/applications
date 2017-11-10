
function tracks = runTrackPostProcessing(data, varargin)

%----------------------------------------------------------------------------
% Figure settings
%----------------------------------------------------------------------------
fset = loadFigureSettings();

%----------------------------------------------------------------------------
% Parse inputs
%----------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Cutoff_f', 2, @isscalar);
ip.addParamValue('Display', false, @islogical);
% Default bounds: [cut-10) [10-20) [20-40) [60-80) [80-100) [100 125) [125-150) [150-end)
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 125 150]);
ip.parse(varargin{:});

cutoff_f = ip.Results.Cutoff_f;

minLft = cutoff_f*data.framerate;
cohortBounds = ip.Results.CohortBounds_s;
cohortBounds(cohortBounds<=minLft) = [];
cohortBounds = [minLft cohortBounds data.movieLength*data.framerate];

% load tracks
load([data.source 'Tracking' filesep 'trackAnalysis.mat']);

%----------------------------------------------------------------------------
% Initialize
%----------------------------------------------------------------------------
% keep only tracks with lifetime above cutoff
tracks = tracks([tracks.end]-[tracks.start]+1 >= cutoff_f);

% master channel
mCh = find(strcmp(data.source, data.channels));

kLevel = norminv(1-0.05/2, 0, 1);

%============================================================================
% I. Assign category to each track
%============================================================================
% Categories:
% Ia)  Single tracks with valid gaps
% Ib)  Single tracks with invalid gaps
% Ic)  Single tracks cut at beginning or end
% Id)  Single tracks, persistent
% IIa) Compound tracks with valid gaps
% IIb) Compound tracks with invalid gaps
% IIc) Compound tracks cut at beginning or end
% IId) Compound tracks, persistent

% The categories correspond to index 1-8, in the above order

validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
singleIdx = [tracks.nSeg]==1;
vis = [tracks.visibility];

mask_Ia = singleIdx & validGaps & vis==1;
mask_Ib = singleIdx & ~validGaps & vis==1;

C = [mask_Ia;
     2*mask_Ib;
     3*(singleIdx & vis==2);
     4*(singleIdx & vis==3);
     5*(~singleIdx & validGaps & vis==1);
     6*(~singleIdx & ~validGaps & vis==1);
     7*(~singleIdx & vis==2);
     8*(~singleIdx & vis==3)];

C = num2cell(sum(C,1));
% assign category
[tracks.catIdx] = deal(C{:});

%============================================================================
% II. Identify diffraction-limited tracks (CCPs)
%============================================================================
% Criterion: if all detected points pass AD-test, then track is a CCP.
% (gaps are not excluded)

% # diffraction-limited points per track (can be different from track length!)
nPl = arrayfun(@(i) nansum(i.hval_AD(mCh,:) .* ~i.gapVect), tracks);

isCCP = num2cell(nPl==0);
[tracks.isCCP] = deal(isCCP{:});
isCCP = [isCCP{:}];

% average mask area per track
meanMaskAreaCCP = arrayfun(@(i) nanmean(i.maskN), tracks(isCCP));
meanMaskAreaNotCCP = arrayfun(@(i) nanmean(i.maskN), tracks(~isCCP));

if ip.Results.Display
    %------------------------------------------------
    % Histograms of mask area after AD test
    %------------------------------------------------
    
    % individual points: separation not as clear
    % meanMaskAreaCCP = [tracks([tracks.isCCP]==1).maskN];
    % meanMaskAreaNotCCP = [tracks([tracks.isCCP]==0).maskN];
    
    % at intensity maximum
    % meanMaskAreaCCP = arrayfun(@(i) i.maskN(i.A(mCh,:)==max(i.A(mCh,:))), tracks([tracks.isCCP]==1));
    % meanMaskAreaNotCCP = arrayfun(@(i) i.maskN(i.A(mCh,:)==max(i.A(mCh,:))), tracks([tracks.isCCP]==0));
    
    da = 1;
    ta = 0:da:da*ceil(max(meanMaskAreaNotCCP)/da);
    histCCP = hist(meanMaskAreaCCP, ta);
    histCCP = histCCP/sum(histCCP)/da;
    histNotCCP = hist(meanMaskAreaNotCCP, ta);
    histNotCCP = histNotCCP/sum(histNotCCP)/da;
    
    figure;
    hold on;
    bar(ta, [histCCP; histNotCCP]', 'BarWidth', 1);
    set(gca, 'XLim', [0 50], fset.axOpts{:}, fset.sfont{:});
    title('Mean mask area / track (pixels)');
    xlabel('# mask pixels');
    hl = legend('CCPs (A-D test)', 'Larger assemblies');
    set(hl, 'Box', 'off');
    % print('-depsc2', 'tracks_ADtest.eps');
    
    %------------------------------------------------
    % Examples for extreme and median percentiles
    %------------------------------------------------
    
%     % Median CCP
%     nm = prctile(meanMaskAreaCCP, 50);
%     % sort as a fct of distance to nm
%     D = (meanMaskAreaCCP-nm).^2;
%     [~,idx] = sort(D);
%     idxCCP = find(isCCP);
%     idxCCP = idxCCP(idx);
%     k = 1;
%     [stack, xa, ya] = getTrackStack(data, tracks(idxCCP(k)));
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
%     print('-depsc', 'track_CCP_p50.eps');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', true);
%     print('-depsc', 'track_CCP_p50_det.eps');
%     plotTrack(data, tracks(idxCCP(k)), 'DisplayMode', 'print', 'YTick', [-20 0:50:350]);
%     print('-depsc', 'track_CCP_p50_track.eps');
%     
%     % 99% CCP
%     nm = prctile(meanMaskAreaCCP, 99);
%     D = (meanMaskAreaCCP-nm).^2;
%     [~,idx] = sort(D);
%     idxCCP = find(isCCP);
%     idxCCP = idxCCP(idx);
%     k = 6;
%     [stack, xa, ya] = getTrackStack(data, tracks(idxCCP(k)), 'Reference', 'frame');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
%     print('-depsc', 'track_CCP_p99.eps');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', true);
%     print('-depsc', 'track_CCP_p99_det.eps');
%     plotTrack(data, tracks(idxCCP(k)), 'DisplayMode', 'print', 'YTick', [-20 0:50:350]);
%     print('-depsc', 'track_CCP_p99_track.eps');
%     
%     % Median other
%     nm = prctile(meanMaskAreaNotCCP, 50);
%     D = (meanMaskAreaNotCCP-nm).^2;
%     [~,idx] = sort(D);
%     idxCCP = find(~isCCP);
%     idxCCP = idxCCP(idx);
%     % k = 5;
%     k = 13;
%     [stack, xa, ya] = getTrackStack(data, tracks(idxCCP(k)), 'Reference', 'frame');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
%     print('-depsc', 'track_notCCP_p50.eps');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', true);
%     print('-depsc', 'track_notCCP_p50_det.eps');
%     plotTrack(data, tracks(idxCCP(k)), 'DisplayMode', 'print', 'YTick', [-20 0:50:350]);
%     print('-depsc', 'track_notCCP_p50_track.eps');
%     
%     % 1% other
%     nm = prctile(meanMaskAreaNotCCP, 1);
%     D = (meanMaskAreaNotCCP-nm).^2;
%     [~,idx] = sort(D);
%     idxCCP = find(~isCCP);
%     idxCCP = idxCCP(idx);
%     k = 4;
%     [stack, xa, ya] = getTrackStack(data, tracks(idxCCP(k)), 'Reference', 'frame');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
%     print('-depsc', 'track_notCCP_p1.eps');
%     plotTrackMontage(tracks(idxCCP(k)), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', true);
%     print('-depsc', 'track_notCCP_p1_det.eps');
%     plotTrack(data, tracks(idxCCP(k)), 'DisplayMode', 'print', 'YTick', [-20 0:50:350]);
%     print('-depsc', 'track_notCCP_p1_track.eps');
    
    
end

%%

v = hist([tracks.catIdx], 1:8);
v = v/numel(tracks)*100;
plotTrackClasses(v');

% % CCPs
% v = hist([tracks([tracks.isCCP]).catIdx], 1:8);
% v = v/numel(tracks)*100;
% plotTrackClasses(v');
% % print('-depsc2', 'trackStats_CCPs.eps');
% 
% % Plaques
% v = hist([tracks(~[tracks.isCCP]).catIdx], 1:8);
% v = v/numel(tracks)*100;
% plotTrackClasses(v');
% % print('-depsc2', 'trackStats_plaques.eps');

%%
% background variance comparison (higher for non-CCPs due to worse fit)

% meanSigmaRCPP = arrayfun(@(i) nanmean(i.A./i.sigma_r), tracks([tracks.isCCP]==1));
% meanSigmaRNotCPP = arrayfun(@(i) nanmean(i.A./i.sigma_r), tracks([tracks.isCCP]==0));
% 
% da = 0.1;
% ta = 0:da:10;
% histCCP = hist(meanSigmaRCPP, ta);
% histCCP = histCCP/sum(histCCP)/da;
% histNotCCP = hist(meanSigmaRNotCPP, ta);
% histNotCCP = histNotCCP/sum(histNotCCP)/da;
% figure;
% hold on;
% bar(ta, [histCCP; histNotCCP]');


%%
% Correlation between A and sigma_r for diffraction-limited structures
% Essentially no correlation
% No change if split along CCP/plaques

% av = [tracks.A];
% sv = [tracks.sigma_r];
% 
% ta = -100:400;
% ts = 0:60;
% n3 = hist3([sv' av' ], {ts, ta});
% figure; 
% % plot(av, sv, 'k.');
% imagesc(ta, ts, n3); axis xy;
% xlabel('A');
% ylabel('\sigma_r');



%%
%----------------------------------------------------------------------------
% Correlation between A and maskN for diffraction-limited structures
%----------------------------------------------------------------------------
% % MC simulation of relationship between maskN and A
% % All tracks
% av = [tracks.A];
% nv = [tracks.maskN];
% % CCP tracks (A-D test)
% % av = [tracks(isCCP).A];
% % nv = [tracks(isCCP).maskN];
% % 
% % av = [tracks(~isCCP).A];
% % nv = [tracks(~isCCP).maskN];
% 
% % av = arrayfun(@(i) nanmean(i.A), tracks(isCCP));
% % nv = arrayfun(@(i) nanmean(i.maskN), tracks(isCCP));
% % av = arrayfun(@(i) nanmean(i.A), tracks(~isCCP));
% % nv = arrayfun(@(i) nanmean(i.maskN), tracks(~isCCP));
% 
% % threshold
% T_mask = prctile(meanMaskAreaCCP, 99);
% 
% 
% % sigma_r = 23;
% Av = 0:10:500;
% % N = 100;
% % sigma_PSF = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, name2wavelength(data.markers{mCh}));
% % % sigma_PSF = 1.5;
% % g0 = simGaussianSpots(13, 13, sigma_PSF, 'x', 7, 'y', 7, 'A', 1);
% % ni = zeros(1,N);
% % ct = zeros(1,numel(Av));
% % ctstd = zeros(1,numel(Av));
% % 
% % for i = 1:numel(Av)
% %     %g = (Av(i) + kLevel*sigma_r) * g0;
% %     g = Av(i) * g0;
% %     for k = 1:N
% %         gn = g + sigma_r*randn(13,13);
% %         [pstruct, mask] = pointSourceDetection(gn, sigma_PSF);
% %         %mask = gn >= kLevel*sigma_r;
% %         CC = bwconncomp(mask);
% %         if ~isempty(pstruct) && ~isempty(CC.PixelIdxList)
% %             ni(k) = max(cellfun(@(i) numel(i), CC.PixelIdxList));
% %         else
% %             ni(k) = 0;
% %         end
% %     end
% %     ct(i) = mean(ni);
% %     ctstd(i) = std(ni);
% % end
% 
% 
% ta = 0:400;
% tn = 0:60;
% n3 = hist3([nv' av' ], {tn, ta});
% % upsample
% n3 = interp1(tn, n3, 1:0.1:60);
% % smooth
% s = 1;
% w = ceil(4*s);
% k = exp(-(-w:w).^2/(2*s^2)); k = k/sum(k);
% n3 = conv2(k, k, padarrayXT(n3, [w w], 'symmetric'), 'valid');
% 
% % n3 = n3/sum(n3(:));
% n3 = n3/max(n3(:));
% 
% % n3 = hist3([nv([tracks.isCCP])' av([tracks.isCCP])' ], {tn, ta});
% % n3 = hist3([nv' (av-kLevel*sv)' ], {tn, ta});
% 
% figure;
% % plot(av, nv, 'k.');
% imagesc(ta, tn, n3); axis xy square;
% axis([0 400 0 40]);
% 
% pc = 0.5;
% nc = 10000;
% cidx = (0:nc-1)/(nc-1);
% cmap = jet(nc);
% cmap = interp1(cidx, cmap, cidx.^pc);
% 
% % caxis([0:0.01:1].^0.5);
% xlabel('Amplitude (A.U.)', fset.lfont{:});
% ylabel('# mask pixels', fset.lfont{:});
% hold on;
% % plot(Av, ct, 'r');
% % plot(Av, ct+ctstd, 'm');
% % plot(Av, ct-ctstd, 'm');
% % plot(Av, T_mask*ones(size(Av)), 'k--', 'LineWidth', 2);
% colormap(cmap);
% hc = colorbar;
% set(gca, fset.sfont{:});%, 'LineWidth', 0);
% set(hc, fset.tfont{:}, 'TickLength', [0 0]);
% set(get(hc, 'YLabel'), 'String', 'Relative frequency', fset.sfont{:});

% print('-depsc', 'A_vs_maskN.eps');

%%
% Correlation between max(A) and maskN at max(A) for diffraction-limited structures
% idx = find([tracks.isCCP]==1);
% N = numel(idx);
% maxA = zeros(1,N);
% maskNatMaxA = zeros(1,N);
% for k = 1:N
%     i = idx(k);
%     iA = tracks(i).A(mCh,:);
%     maxIdx = find(iA==max(iA), 1, 'first');
%     maxA(k) = iA(maxIdx);
%     maskNatMaxA(k) = tracks(i).maskN(maxIdx);
% end
% 
% n3 = hist3([maskNatMaxA' maxA' ], {tn, ta});
% figure;
% imagesc(ta, tn, n3.^0.25); axis xy square;
% xlabel('A');
% ylabel('# mask pixels');
% % plot(maxA, maskNatMaxA, 'k.');




%%
%============================================================================
% III. Map gaps
%============================================================================

%----------------------------------------------------------------------------
% Track classes
%----------------------------------------------------------------------------

% Reference distribution: class Ia tracks
idx_Ia = find(mask_Ia);
idx_Ib = find(mask_Ib);

trackLengths = [tracks.end]-[tracks.start]+1;

%----------------------------------------------------------------------------
% Reference distributions for each cohort from class Ia tracks
%----------------------------------------------------------------------------
% Determine critical max. intensity values from class Ia tracks

% # cohorts
nc = numel(cohortBounds)-1;

% max intensities of all 'Ia' tracks
maxInt = arrayfun(@(i) max(i.A(mCh,:)), tracks(idx_Ia));
maxIntDistr = cell(1,nc);
criticalInt = zeros(nc,2);
lft_Ia = [tracks(idx_Ia).lifetime_s];
for i = 1:nc
    maxIntDistr{i} = maxInt(cohortBounds(i)<=lft_Ia & lft_Ia<cohortBounds(i+1));
    % percentiles (critical values for test)
    criticalInt(i,:) = prctile(maxIntDistr{i}, [2.5 97.5]);
end
% Plot lower/upper intensity bounds
% figure; hold on; plot(criticalInt(:,1)); plot(criticalInt(:,2),'r')

if ip.Results.Display
    %----------------------------------------------------------------------------
    % Gap distribution, 'Ia' tracks
    %----------------------------------------------------------------------------
    N = numel(idx_Ia);
    
    % Matrices of all 'Ia' tracks aligned to beginning & end
    Ms = NaN(N, data.movieLength);
    Me = NaN(N, data.movieLength);
    
    for k = 1:N
        i = idx_Ia(k);
        % gap frames: tracks(i).f(tracks(i).gapVect==1)
        Ms(k,1:trackLengths(i)) = 0;
        Ms(k,tracks(i).f(tracks(i).gapVect==1) - tracks(i).start + 1) = 1;
        % aligned at track end
        Me(k, end-trackLengths(i)+1:end) = 0;
        Me(k, end-trackLengths(i) + tracks(i).f(tracks(i).gapVect==1) - tracks(i).start + 1) = 1;
    end
    
    % global distributions not very informative (short tracks ambiguous) >>> broken down by cohort:
    colorV = jet(nc);
    cohortLabel = arrayfun(@(i) [' ' num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-data.framerate) ' s'], 1:nc, 'UniformOutput', false);
    
    figure('Color', 'w');
    ha1 = subplot(1,2,1); hold on;
    ha2 = subplot(1,2,2); hold on;
    hp = zeros(1,nc+1);
    for i = nc:-1:1 % plot shorter cohorts in front
        cIdx = (cohortBounds(i)<=trackLengths(idx_Ia) & trackLengths(idx_Ia)<cohortBounds(i+1));
        pGapStart = nansum(Ms(cIdx,:),1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
        pGapEnd = nansum(Me(cIdx,:),1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
        hp(i) = plot(ha1, pGapStart, 'Color', colorV(i,:), 'LineWidth', 1.5);
        plot(ha2, pGapEnd, 'Color', colorV(i,:), 'LineWidth', 1.5);
    end
    % all tracks
    pGapStart = nansum(Ms,1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
    pGapEnd = nansum(Me,1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
    hp(nc+1) = plot(ha1, pGapStart, '--', 'Color', 'k', 'LineWidth', 1.5);
    plot(ha2, pGapEnd, '--', 'Color', 'k', 'LineWidth', 1.5);
    
    
    
    set(ha1, 'XLim', [1 81], fset.axOpts{:}, 'XTick', 1+[1 20:20:80], 'XTickLabel', [1 20:20:80], 'YLim', [0 0.305]);
    set(ha2, 'XLim', [data.movieLength-80 data.movieLength], 'YAxisLocation', 'right',...
        fset.axOpts{:},...
        'XTick', data.movieLength - [80:-20:20 1], 'XTickLabel', [80:-20:20 1], 'YLim', [0 0.305]);
    axes(ha1);
    
    % legend
    hl = legend(hp, [cohortLabel {' all tracks'}], 'Location', 'NorthEast');
    set(hl, 'Box', 'off');
    pos = get(hl, 'Position');
    pos(1) = 0.5-pos(3)/2;
    set(hl, 'Position', pos);
    ht = get(hl, 'title');
    set(ht, 'string', 'Lifetime cohort', 'Visible', 'on', 'FontSize', 12);
    
    xlabel(ha1, '# frames after start', fset.tfont{:});
    ylabel(ha1, 'Gap frequency', fset.tfont{:});
    xlabel(ha2, '# frames before end', fset.tfont{:});
    
    % print('-depsc2', 'GapFrequency_Ia.eps');
    
    %%
    % % control w/o matrices: last 20 frames of tracks in 40-59 cohort
    % idx = idx_Ia(40<=trackLengths(idx_Ia) & trackLengths(idx_Ia)<60);
    % tmp = arrayfun(@(i) i.gapVect(1,end-39:end), tracks(idx), 'UniformOutput', false);
    % alt = sum(vertcat(tmp{:}),1);
    % figure; plot(alt);
    
    
    % examples of tracks with gaps at end-20, end-40, end-60 in cohorts 3-5, resp.
    % track in 60-80 cohort: 58 frames before end
    % cIdx = (60<=trackLengths(idx_Ia) & trackLengths(idx_Ia)<80);
    % cand = find(Me(:, end-58)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ia(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % % plotTrack(data, tracks(idx_Ia(cand(6))), 'DisplayMode', 'print', 'YTick', -20:20:180); % -> must be split in two!!
    % print('-depsc2', 'gapExample60-80_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ia(cand(1))));
    % plotTrackMontage(tracks(idx_Ia(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample60-80_montage.eps');
    %
    % cIdx = (40<=trackLengths(idx_Ia) & trackLengths(idx_Ia)<60);
    % cand = find(Me(:, end-39)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ia(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % print('-depsc2', 'gapExample40-60_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ia(cand(1))));
    % plotTrackMontage(tracks(idx_Ia(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample40-60_montage.eps');
    %
    %
    % cIdx = (20<=trackLengths(idx_Ia) & trackLengths(idx_Ia)<40);
    % cand = find(Me(:, end-23)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ia(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % print('-depsc2', 'gapExample20-40_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ia(cand(1))));
    % plotTrackMontage(tracks(idx_Ia(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample20-40_montage.eps');
    
    
    %%
    % Rinse, repeat:
    %----------------------------------------------------------------------------
    % Gap distribution, 'Ib' tracks, all gaps
    %----------------------------------------------------------------------------
    N = numel(idx_Ib);
    
    % Matrices of all 'Ib' tracks aligned to beginning & end
    Ms = NaN(N, data.movieLength);
    Me = NaN(N, data.movieLength);
    
    for k = 1:N
        i = idx_Ib(k);
        % gap frames: tracks(i).f(tracks(i).gapVect==1)
        Ms(k,1:trackLengths(i)) = 0;
        Ms(k,tracks(i).f(tracks(i).gapVect==1) - tracks(i).start + 1) = 1;
        % aligned at track end
        Me(k, end-trackLengths(i)+1:end) = 0;
        Me(k, end-trackLengths(i) + tracks(i).f(tracks(i).gapVect==1) - tracks(i).start + 1) = 1;
    end
    
    % global distributions not very informative (short tracks ambiguous) >>> broken down by cohort:
    
    figure('Color', 'w');
    ha1 = subplot(1,2,1); hold on;
    ha2 = subplot(1,2,2); hold on;
    hp = zeros(1,nc+1);
    for i = nc:-1:1 % plot shorter cohorts in front
        cIdx = (cohortBounds(i)<=trackLengths(idx_Ib) & trackLengths(idx_Ib)<cohortBounds(i+1));
        pGapStart = nansum(Ms(cIdx,:),1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
        pGapEnd = nansum(Me(cIdx,:),1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
        hp(i) = plot(ha1, pGapStart, 'Color', colorV(i,:), 'LineWidth', 1.5);
        plot(ha2, pGapEnd, 'Color', colorV(i,:), 'LineWidth', 1.5);
    end
    % all tracks
    pGapStart = nansum(Ms,1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
    pGapEnd = nansum(Me,1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
    hp(nc+1) = plot(ha1, pGapStart, '--', 'Color', 'k', 'LineWidth', 1.5);
    plot(ha2, pGapEnd, '--', 'Color', 'k', 'LineWidth', 1.5);
    
    set(ha1, 'XLim', [1 81], fset.axOpts{:}, 'LineWidth', 1.5, 'XTick', 1+[1 20:20:80], 'XTickLabel', [1 20:20:80], 'YLim', [0 0.305]);
    set(ha2, 'XLim', [data.movieLength-80 data.movieLength], 'YAxisLocation', 'right',...
        fset.axOpts{:},...
        'XTick', data.movieLength - [80:-20:20 1], 'XTickLabel', [80:-20:20 1], 'YLim', [0 0.305]);
    axes(ha1);
    
    % legend
    hl = legend(hp, [cohortLabel {' all tracks'}], 'Location', 'NorthEast');
    set(hl, 'Box', 'off');
    pos = get(hl, 'Position');
    pos(1) = 0.5-pos(3)/2;
    set(hl, 'Position', pos);
    ht = get(hl, 'title');
    set(ht, 'string', 'Lifetime cohort', 'Visible', 'on', 'FontSize', 12);
    
    xlabel(ha1, '# frames after start', fset.tfont{:});
    ylabel(ha1, 'Gap frequency', fset.tfont{:});
    xlabel(ha2, '# frames before end', fset.tfont{:});
    
    %print('-depsc2', 'GapFrequency_Ib_allGaps.eps');
    
    %%
    % cIdx = (60<=trackLengths(idx_Ib) & trackLengths(idx_Ib)<80);
    % cand = find(Me(:, end-58)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ib(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % print('-depsc2', 'gapExample60-80_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ib(cand(1))));
    % plotTrackMontage(tracks(idx_Ib(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample60-80_montage.eps');
    %
    % cIdx = (40<=trackLengths(idx_Ib) & trackLengths(idx_Ib)<60);
    % cand = find(Me(:, end-39)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ib(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % print('-depsc2', 'gapExample40-60_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ib(cand(1))));
    % plotTrackMontage(tracks(idx_Ib(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample40-60_montage.eps');
    %
    %
    % cIdx = (20<=trackLengths(idx_Ib) & trackLengths(idx_Ib)<40);
    % cand = find(Me(:, end-23)==1 & cIdx');
    % plotTrack(data, tracks(idx_Ib(cand(1))), 'DisplayMode', 'print', 'YTick', -20:20:180);
    % print('-depsc2', 'gapExample20-40_track.eps');
    % [stack, xa, ya] = getTrackStack(data, tracks(idx_Ib(cand(1))));
    % plotTrackMontage(tracks(idx_Ib(cand(1))), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
    % print('-depsc2', 'gapExample20-40_montage.eps');
    
    
    
    
    %%
    % %----------------------------------------------------------------------------
    % % Gap distribution, 'Ib' tracks, rejected gaps only
    % %----------------------------------------------------------------------------
    % N = numel(idx_Ib);
    %
    % % Matrices of all 'Ib' tracks aligned to beginning & end
    % Ms = NaN(N, data.movieLength);
    % Me = NaN(N, data.movieLength);
    %
    % for k = 1:N
    %     i = idx_Ib(k);
    %     % gap frames: tracks(i).gapIdx{tracks(i).gapStatus==5}
    %     Ms(k, 1:trackLengths(i)) = 0;
    %     Ms(k, [tracks(i).gapIdx{tracks(i).gapStatus==5}]) = 1;
    %
    %     % aligned at track end
    %     Me(k, end-trackLengths(i)+1:end) = 0;
    %     Me(k, end-trackLengths(i) + [tracks(i).gapIdx{tracks(i).gapStatus==5}]) = 1;
    % end
    %
    % % global distributions not very informative (short tracks ambiguous) >>> broken down by cohort:
    %
    % figure('Color', 'w');
    % ha1 = subplot(1,2,1); hold on;
    % ha2 = subplot(1,2,2); hold on;
    % hp = zeros(1,nc+1);
    % for i = nc:-1:1 % plot shorter cohorts in front
    %     cIdx = (cohortBounds(i)<=trackLengths(idx_Ib) & trackLengths(idx_Ib)<cohortBounds(i+1));
    %     pGapStart = nansum(Ms(cIdx,:),1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
    %     pGapEnd = nansum(Me(cIdx,:),1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
    %     hp(i) = plot(ha1, pGapStart, 'Color', colorV(i,:), 'LineWidth', 1.5);
    %     plot(ha2, pGapEnd, 'Color', colorV(i,:), 'LineWidth', 1.5);
    % end
    % % all tracks
    % pGapStart = nansum(Ms,1); pGapStart([1 end]) = NaN; pGapStart = pGapStart/nansum(pGapStart);
    % pGapEnd = nansum(Me,1); pGapEnd([1 end]) = NaN; pGapEnd = pGapEnd/nansum(pGapEnd);
    % hp(nc+1) = plot(ha1, pGapStart, '--', 'Color', 'k', 'LineWidth', 1.5);
    % plot(ha2, pGapEnd, '--', 'Color', 'k', 'LineWidth', 1.5);
    %
    % set(ha1, 'XLim', [1 81], plotOpts{:}, 'LineWidth', 1.5, 'XTick', 1+[1 20:20:80], 'XTickLabel', [1 20:20:80], 'YLim', [0 0.305]);
    % set(ha2, 'XLim', [data.movieLength-80 data.movieLength], 'YAxisLocation', 'right',...
    %     plotOpts{:},...
    %     'XTick', data.movieLength - [80:-20:20 1], 'XTickLabel', [80:-20:20 1], 'YLim', [0 0.305]);
    % axes(ha1);
    %
    % % legend
    % hl = legend(hp, [cohortLabel {' all tracks'}], 'Location', 'NorthEast');
    % set(hl, 'Box', 'off');
    % pos = get(hl, 'Position');
    % pos(1) = 0.5-pos(3)/2;
    % set(hl, 'Position', pos);
    % ht = get(hl, 'title');
    % set(ht, 'string', 'Lifetime cohort', 'Visible', 'on', 'FontSize', 12);
    %
    % xlabel(ha1, '# frames after start', fset.tfont{:});
    % ylabel(ha1, 'Gap frequency', fset.tfont{:});
    % xlabel(ha2, '# frames before end', fset.tfont{:});
    %
    % print('-depsc2', 'GapFrequency_Ib_rejectedGaps.eps');


    %----------------------------------------------------------------------------
    % Gap distance distribution, Ia vs. Ib tracks
    %----------------------------------------------------------------------------
    hasGaps = arrayfun(@(i) ~isempty(i.gapIdx) && numel(i.gapIdx)>1, tracks);
    gapDist_Ia = arrayfun(@(i) diff([i.gapIdx{:}])-1, tracks(mask_Ia & hasGaps), 'UniformOutput', false);
    gapDist_Ia = [gapDist_Ia{:}]; gapDist_Ia(gapDist_Ia<1) = [];
    
    gapDist_Ib = arrayfun(@(i) diff([i.gapIdx{:}])-1, tracks(mask_Ib & hasGaps), 'UniformOutput', false);
    gapDist_Ib = [gapDist_Ib{:}]; gapDist_Ib(gapDist_Ib<1) = [];
    
    ti = 1:50;
    hist_Ia = hist(gapDist_Ia, ti) / sum(mask_Ia & hasGaps);
    hist_Ib = hist(gapDist_Ib, ti) / sum(mask_Ib & hasGaps);
    figure;
    hold on;
    plot(ti, hist_Ia, 'k', 'LineWidth', 1.5);
    plot(ti, hist_Ib, 'r', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 15], fset.axOpts{:}, 'XTick', 1:20);
    xlabel('# frames between consecutive gaps', fset.tfont{:});
    ylabel('Occurence / track with > 1 gaps', fset.tfont{:});
    hl = legend('Single tracks (cat. Ia)', 'Single tracks, rej. gaps (cat. Ib)', 'Location', 'NorthEast');
    set(hl, 'Box', 'off');
    
    % print('-depsc2', 'GapDistStatistics.eps');
    
    
    %----------------------------------------------------------------------------
    % Gap occurrence, per cohort, Ia & Ib tracks
    %----------------------------------------------------------------------------
    gapsPerTrack_Ia = cell(1,nc);
    gapsPerTrack_Ib = cell(1,nc);
    
    for i = 1:nc
        cIdx = (cohortBounds(i)<=trackLengths & trackLengths<cohortBounds(i+1));
        gapsPerTrack_Ia{i} = arrayfun(@(i) sum(i.gapVect), tracks(mask_Ia & cIdx));
        gapsPerTrack_Ib{i} = arrayfun(@(i) sum(i.gapVect), tracks(mask_Ib & cIdx));
    end
    figure;
    barplot2([cellfun(@(i) mean(i), gapsPerTrack_Ia); cellfun(@(i) mean(i), gapsPerTrack_Ib)]',...
        [cellfun(@(i) std(i), gapsPerTrack_Ia); cellfun(@(i) std(i), gapsPerTrack_Ib)]',...
        'FaceColor', fset.cfTrackClasses(1:2,:), 'EdgeColor', fset.ceTrackClasses(1:2,:),...
        'XLabels', cohortLabel, 'XLabel', 'Lifetime cohort', 'YLabel', 'gap frames/track', 'YLim', [0 12]);
    
    hl = legend('Tracks w/ valid gaps', 'Tracks w/ rejected gaps', 'Location', 'NorthWest');
    set(hl, 'Box', 'off', fset.ifont{:});
    % print('-depsc2', 'GapFramesPerTrack.eps');

    
    %----------------------------------------------------------------------------
    % Gap occurrence, per cohort, Ia & Ib tracks (split into valid/rejected gaps)
    %----------------------------------------------------------------------------
    gapsPerTrack_Ia = cell(1,nc);
    gapsPerTrack_Ib4 = cell(1,nc);
    gapsPerTrack_Ib5 = cell(1,nc);
    
    for i = 1:nc
        cIdx = (cohortBounds(i)<=trackLengths & trackLengths<cohortBounds(i+1));
        gapsPerTrack_Ia{i} = arrayfun(@(i) numel(i.gapStatus), tracks(mask_Ia & cIdx));
        gapsPerTrack_Ib4{i} = arrayfun(@(i) sum(i.gapStatus==4), tracks(mask_Ib & cIdx));
        gapsPerTrack_Ib5{i} = arrayfun(@(i) sum(i.gapStatus==5), tracks(mask_Ib & cIdx));
    end
    figure;
    barplot2([cellfun(@(i) mean(i), gapsPerTrack_Ia);
        cellfun(@(i) mean(i), gapsPerTrack_Ib4);
        cellfun(@(i) mean(i), gapsPerTrack_Ib5)]',...
        [cellfun(@(i) std(i), gapsPerTrack_Ia);
        cellfun(@(i) std(i), gapsPerTrack_Ib4);
        cellfun(@(i) std(i), gapsPerTrack_Ib5)]',...
        'FaceColor', fset.cfTrackClasses([1 4 5],:), 'EdgeColor', fset.ceTrackClasses([1 4 5],:),...
        'XLabels', cohortLabel, 'XLabel', 'Lifetime cohort', 'YLabel', '#gaps / track', 'YLim', [0 12]);
    
    hl = legend('Tracks w/ valid gaps', 'Tracks w/ rejected gaps, valid gaps', 'Tracks w/ rejected gaps, rejected gaps', 'Location', 'NorthWest');
    set(hl, 'Box', 'off', fset.ifont{:});
    %print('-depsc2', 'GapsPerTrack.eps');

    
    %----------------------------------------------------------------------------
    % max int. dist. Ia vs Ib
    %----------------------------------------------------------------------------
    
    maxIntensities_Ia = cell(1,nc);
    maxIntensities_Ib = cell(1,nc);
    maxIntensities_Ia_start = cell(1,nc);
    maxIntensities_Ib_start = cell(1,nc);
    minIntensities_Ia_start = cell(1,nc);
    minIntensities_Ib_start = cell(1,nc);
    
    T = 4;
    
    for i = 1:nc
        cIdx = (cohortBounds(i)<=trackLengths & trackLengths<cohortBounds(i+1));
        maxIntensities_Ia{i} = arrayfun(@(i) max(i.A(mCh,:)), tracks(mask_Ia & cIdx));
        maxIntensities_Ib{i} = arrayfun(@(i) max(i.A(mCh,:)), tracks(mask_Ib & cIdx));
        
        maxIntensities_Ia_start{i} = arrayfun(@(i) max(i.A(mCh,1:T)), tracks(mask_Ia & cIdx));
        maxIntensities_Ib_start{i} = arrayfun(@(i) max(i.A(mCh,1:T)), tracks(mask_Ib & cIdx));
        minIntensities_Ia_start{i} = arrayfun(@(i) min(i.A(mCh,1:T)), tracks(mask_Ia & cIdx));
        minIntensities_Ib_start{i} = arrayfun(@(i) min(i.A(mCh,1:T)), tracks(mask_Ib & cIdx));
    end
    
    M_Ia = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', maxIntensities_Ia, 'UniformOutput', false);
    M_Ia = [M_Ia{:}];
    M_Ib = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', maxIntensities_Ib, 'UniformOutput', false);
    M_Ib = [M_Ib{:}];
    figure; h = boxplot2({M_Ia, M_Ib}, 'BarWidth', 0.8,...
        'FaceColor', fset.cfTrackClasses(1:2,:), 'EdgeColor', fset.ceTrackClasses(1:2,:),...
        'XLabels', cohortLabel, 'XLabel', 'Lifetime cohort', 'YLabel', 'Maximum intensity (A.U.)');
    hl = legend(h, 'Tracks w/ valid gaps', 'Tracks w/ rejected gaps', 'Location', 'NorthWest');
    set(hl, 'Box', 'off', fset.ifont{:});
    % print('-depsc2', 'MaxIntensityBoxPlot_Ia_vs_Ib.eps');
    
    
    M_Ia = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', maxIntensities_Ia_start, 'UniformOutput', false);
    M_Ia = [M_Ia{:}];
    M_Ib = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', maxIntensities_Ib_start, 'UniformOutput', false);
    M_Ib = [M_Ib{:}];
    figure; h = boxplot2({M_Ia, M_Ib}, 'BarWidth', 0.8,...
        'FaceColor', fset.cfTrackClasses(1:2,:), 'EdgeColor', fset.ceTrackClasses(1:2,:),...
        'XLabels', cohortLabel, 'XLabel', 'Lifetime cohort', 'YLabel', ['Maximum intensity, first ' num2str(T) ' frames (A.U.)']);
    hl = legend(h, 'Tracks w/ valid gaps', 'Tracks w/ rejected gaps', 'Location', 'NorthWest');
    set(hl, 'Box', 'off', fset.ifont{:});
    % print('-depsc2', 'MaxIntensityBoxPlot_Ia_vs_Ib_frames_1-4.eps');
    
    M_Ia = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', minIntensities_Ia_start, 'UniformOutput', false);
    M_Ia = [M_Ia{:}];
    M_Ib = cellfun(@(i) prctile(i, [50 25 75 2.5 97.5])', minIntensities_Ib_start, 'UniformOutput', false);
    M_Ib = [M_Ib{:}];
    figure; h = boxplot2({M_Ia, M_Ib}, 'BarWidth', 0.8,...
        'FaceColor', fset.cfTrackClasses(1:2,:), 'EdgeColor', fset.ceTrackClasses(1:2,:),...
        'XLabels', cohortLabel, 'XLabel', 'Lifetime cohort', 'YLabel', ['Minimum intensity, first ' num2str(T) ' frames (A.U.)']);
    hl = legend(h, 'Tracks w/ valid gaps', 'Tracks w/ rejected gaps', 'Location', 'NorthWest');
    set(hl, 'Box', 'off', fset.ifont{:});
    % print('-depsc2', 'MinIntensityBoxPlot_Ia_vs_Ib_frames_1-4.eps');
    
    
    plotTrackFirstFrameIntensities(data, tracks(idx_Ib), 'YLim', [0 80], 'Print', true);
    plotTrackFirstFrameIntensities(data, tracks(idx_Ib), 'YLim', [0 80], 'Print', true);
end

%%
%----------------------------------------------------------------------------
% Criteria for mapping
%----------------------------------------------------------------------------
% - max intensity must be within 2.5th percentile of max. intensity distribution for 'Ia' tracks
% - lifetime > 4 frames (at 4 frames: track = [x o o x])

% maxIntensities_Ia = cell(1,nc);
% for i = 1:nc
%     cIdx = (cohortBounds(i)<=trackLengths & trackLengths<cohortBounds(i+1));
%     maxIntensities_Ia{i} = arrayfun(@(i) max(i.A(mCh,:)), tracks(mask_Ia & cIdx));
% end

% maxIntDistr = maxIntensities_Ia -> clean up
mappingThresholdMaxInt = cellfun(@(i) prctile(i, 2.5), maxIntDistr);

%----------------------------------------------------------------------------
% Perform mapping
%----------------------------------------------------------------------------
% get lifetime histograms before change
lftHistsBefore = getLifetimeHistogram(data, tracks);
% tracksX = tracks;
% set all gap status values to '4' for tracks that match criteria
for k = 1:numel(idx_Ib);
    i = idx_Ib(k);

    % get cohort idx (logical)
    cIdx = cohortBounds(1:nc)<=tracks(i).lifetime_s & tracks(i).lifetime_s<cohortBounds(2:nc+1);
    
    if max(tracks(i).A(mCh,:)) >= mappingThresholdMaxInt(cIdx) && trackLengths(i)>4
        %tracksX(i).gapStatus(tracksX(i).gapStatus==5) = 4;
        tracks(i).catIdx = 1;
    end
end
lftHistsAfter = getLifetimeHistogram(data, tracks);

%----------------------------------------------------------------------------
% Display lifetime histograms before (Ia/Ib) and after
%----------------------------------------------------------------------------

if ip.Results.Display
    figure;
    hold on;
    plot(lftHistsBefore.t, lftHistsBefore.Ia, 'k.-');
    plot(lftHistsBefore.t, lftHistsBefore.Ib, 'g.-');
    plot(lftHistsAfter.t, lftHistsAfter.Ia, 'b.--');
    plot(lftHistsAfter.t, lftHistsAfter.Ib, 'r.--');
    set(gca, fset.axOpts{:}, fset.tfont{:});
    axis([0 100 0 0.15]);
    xlabel('Lifetime (s)', fset.tfont{:});
    ylabel('Frequency', fset.tfont{:});
    hl = legend('Valid tracks, before post-proc.', 'Rejected tracks, before post-proc.',...
           'Valid tracks, after post-proc.', 'Rejected tracks, after post-proc.');
    set(hl, 'Box', 'off', fset.ifont{:});
end

v = hist([tracks.catIdx], 1:8);
v = v/numel(tracks);
plotTrackClasses(v');

%%
%============================================================================
% III. Check buffer values
%============================================================================
% Condition: 2 frames must be below background noise level in both start and end buffer
Tbuffer = 2;


% loop through cat. Ia tracks
idx_Ia = find([tracks.catIdx]==1);
for k = 1:numel(idx_Ia)
    i = idx_Ia(k);
    
    % 'H' value: pval_Ar < 0.05; H0: A = background
    %tracks(i).startBuffer.pval_Ar
    %hval_Ar = tracks(i).startBuffer.pval_Ar < 0.05;
    %hval_Ar = pval_Ar < 0.05;
    
    % quick fix:
    if sum(tracks(i).startBuffer.A < tracks(i).startBuffer.sigma_r*kLevel) < Tbuffer ||...
            sum(tracks(i).endBuffer.A < tracks(i).endBuffer.sigma_r*kLevel) < Tbuffer
        tracks(i).catIdx = 2;
    end
    
end

v = hist([tracks.catIdx], 1:8);
v = v/numel(tracks);
plotTrackClasses(v');

% Examples:
% plotTrack(data, tracks(8452), 'DisplayMode', 'print'); print('-depsc2', 'track8452.eps');
% plotTrack(data, tracks(8451), 'DisplayMode', 'print'); print('-depsc2', 'track8451.eps');

% [stack, xa, ya] = getTrackStack(data, tracks(8451));
% plotTrackMontage(tracks(8451), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);
% 
% [stack, xa, ya] = getTrackStack(data, tracks(8452));
% plotTrackMontage(tracks(8452), stack, xa, ya, 'ShowMarkers', true, 'ShowDetection', false);

%%
%============================================================================
% IV. Switch ~isCCP tracks from Cat. Ia to Ib
%============================================================================
[tracks([tracks.catIdx]==1 & ~isCCP).catIdx] = deal(2);
v = hist([tracks.catIdx], 1:8);
v = v/numel(tracks);
plotTrackClasses(v');


%%
%============================================================================
% V. Apply threshold on max. intensity
%============================================================================




%%
%============================================================================
% VI. Cut tracks with sequential events (hotspots) into individual tracks
%============================================================================





%%
%============================================================================
% VII. Save output
%============================================================================
% ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
save([data.source 'Tracking' filesep 'tracksPost.mat'], 'tracks');




