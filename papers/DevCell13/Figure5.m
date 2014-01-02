%=========================================================================================
% Figure 5
%=========================================================================================
% This script generates the panels for Figure 5 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f5path = '/Users/aguet/Documents/Papers/DevCell13/Figure 5 - Dynamin/';

%%
%===========================================================
% 3rd data set: SK-MEL-2 w/ Dyn2-EGFP and tdTomato-LCa O/X
%===========================================================
%dpath = '/home/fa48/lccb-endocytosis/Marcel_fs003/121003_SKMEL2_dyn_dtTomLca/';
dpath = '~/Documents/MATLAB/Data/121003_SKMEL2_dyn_dtTomLca/';

dataDC3 = loadConditionData(dpath, {'568', '488'}, {'tdTomato', 'EGFP'}, 'Parameters', [1.49 100 6.7e-6]);
%runDetection(dataDC3, 'Overwrite', true);
%runTracking(dataDC3, loadTrackSettings(), 'Overwrite', true);
%runTrackProcessing(dataDC3, 'Overwrite', true);
%runSlaveChannelClassification(dataDC3, 'Overwrite', true, 'Mode', 'random');

%%
%=====================================================
% Panel A: Track examples
%=====================================================
tracks1 = loadTracks(dataDC3(1), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);
% tracks2 = loadTracks(dataDC3(2), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);
% tracks3 = loadTracks(dataDC3(3), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);
% tracks5 = loadTracks(dataDC3(5), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);
% tracks6 = loadTracks(dataDC3(6), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);
tracks7 = loadTracks(dataDC3(7), 'Category', 'Ia', 'Mask', true, 'Cutoff_f', 5);

ms = [6 2.5 0.5]; % marker sizes
sf = 4; % scaling factor relative to clathrin
ya = -25:25:125;

%------------------------------------
% Track A
%------------------------------------
track1 = tracks1(356);
track1.A(1,:) = sf*track1.A(1,:);
track1.A_pstd(1,:) = sf*track1.A_pstd(1,:);
track1.c(1,:) = sf*track1.c(1,:);
track1.c_pstd(1,:) = sf*track1.c_pstd(1,:);
track1.sigma_r(1,:) = sf*track1.sigma_r(1,:);
track1.startBuffer.A(1,:) = sf*track1.startBuffer.A(1,:);
track1.startBuffer.A_pstd(1,:) = sf*track1.startBuffer.A_pstd(1,:);
track1.startBuffer.c(1,:) = sf*track1.startBuffer.c(1,:);
track1.startBuffer.sigma_r(1,:) = sf*track1.startBuffer.sigma_r(1,:);
track1.endBuffer.A(1,:) = sf*track1.endBuffer.A(1,:);
track1.endBuffer.A_pstd(1,:) = sf*track1.endBuffer.A_pstd(1,:);
track1.endBuffer.c(1,:) = sf*track1.endBuffer.c(1,:);
track1.endBuffer.sigma_r(1,:) = sf*track1.endBuffer.sigma_r(1,:);

figure(fset.fOpts{:}, 'Position', [12 12 5.75 3.25]); % 8 5.5
ha = axes(fset.axOpts{:}, 'Position', [1.5 1 4 2]); % 6 3.5 
hold on;
plotTrack(dataDC3(1), track1, 1, 'Handle', ha, 'DisplayMode', 'print',...
    'MarkerSizes', ms,...
    'XLim', [-12 68+12], 'YTick', ya);
set(gca, 'TickLength', fset.TickLength*1.5);
ylabel('Fluo. int. (A.U.)', fset.lfont{:});
% print('-depsc', '-loose', 'track356_ch1.eps');
cla
plotTrack(dataDC3(1), track1, 2, 'Handle', ha, 'DisplayMode', 'print',...
    'MarkerSizes', ms, 'BackgroundConfidence', norminv(0.95,0,1)*10.5871, ...
    'XLim', [-12 68+12], 'YTick', ya);
set(gca, 'YTickLabel', [], 'TickLength', fset.TickLength*1.5);
ylabel('', fset.lfont{:});
% print('-depsc', '-loose', 'track356_ch2_bgConf.eps');

%------------------------------------
% Track B
%------------------------------------
track2 = tracks7(181);
track2.A(1,:) = sf*track2.A(1,:);
track2.A_pstd(1,:) = sf*track2.A_pstd(1,:);
track2.c(1,:) = sf*track2.c(1,:);
track2.sigma_r(1,:) = sf*track2.sigma_r(1,:);
track2.startBuffer.A(1,:) = sf*track2.startBuffer.A(1,:);
track2.startBuffer.A_pstd(1,:) = sf*track2.startBuffer.A_pstd(1,:);
track2.startBuffer.c(1,:) = sf*track2.startBuffer.c(1,:);
track2.startBuffer.sigma_r(1,:) = sf*track2.startBuffer.sigma_r(1,:);
track2.endBuffer.A(1,:) = sf*track2.endBuffer.A(1,:);
track2.endBuffer.A_pstd(1,:) = sf*track2.endBuffer.A_pstd(1,:);
track2.endBuffer.c(1,:) = sf*track2.endBuffer.c(1,:);
track2.endBuffer.sigma_r(1,:) = sf*track2.endBuffer.sigma_r(1,:);

figure(fset.fOpts{:}, 'Position', [12 12 5.75 3.25]); % 8 5.5
ha = axes(fset.axOpts{:}, 'Position', [1.5 1 4 2]); % 6 3.5 
hold on;
plotTrack(dataDC3(7), track2, 1, 'Handle', ha, 'DisplayMode', 'print',...
    'MarkerSizes', ms,...
    'XLim', [-12 68+12], 'YTick', ya);
set(gca, 'TickLength', fset.TickLength*1.5);
ylabel('Fluo. int. (A.U.)', fset.lfont{:});
xlabel('Time (s)', fset.lfont{:});
% print('-depsc', '-loose', 'track181_ch1.eps');
cla
plotTrack(dataDC3(7), track2, 2, 'Handle', ha, 'DisplayMode', 'print',...
    'MarkerSizes', ms, 'BackgroundConfidence', norminv(0.95,0,1)*10.5283,...
    'XLim', [-12 68+12], 'YTick', ya);
set(gca, 'YTickLabel', [], 'TickLength', fset.TickLength*1.5);
ylabel('', fset.lfont{:});
xlabel('Time (s)', fset.lfont{:});
% print('-depsc', '-loose', 'track181_ch2_bgConf.eps');

[stack1, xa1, ya1] = getTrackStack(dataDC3(1), track1);
[stack2, xa2, ya2] = getTrackStack(dataDC3(7), track2);
tmp1 = double(cat(3,stack1{1,:},stack2{1,:}));
tmp2 = double(cat(3,stack1{2,:},stack2{2,:}));
dRange = [min(tmp1(:)) max(tmp1(:)); min(tmp2(:)) max(tmp2(:))];
plotTrackMontage(track1, stack1, xa1, ya1, 'Width', 200, 'FramesPerRow', 20, 'DynamicRange', dRange);
% print('-depsc', '-loose', 'track356_montage.eps');
plotTrackMontage(track2, stack2, xa2, ya2, 'Width', 200, 'FramesPerRow', 20, 'DynamicRange', dRange);
% print('-depsc', '-loose', 'track181_montage.eps');


% % cumulative averaged montages (doesn't yield a visual difference)
% w = 1;
% nz = size(stack1,2);
% stack1int = cell(2,nz);
% for z = 1+w:nz-w
%     stack1int{1,z} = mean(cat(3, stack1{1,z-w:z+w}),3);
%     stack1int{2,z} = mean(cat(3, stack1{2,z-w:z+w}),3);
% end
% [stack1int{:,[1:w nz-w+1:nz]}] = deal(uint16(65535*ones(17,17)));
% 
% plotTrackMontage(track1, stack1, xa1, ya1, 'Width', 600, 'FramesPerRow', 20, 'DynamicRange', dRange);
% plotTrackMontage(track1, stack1int, xa1, ya1, 'Width', 600, 'FramesPerRow', 20, 'DynamicRange', dRange);


%%
%=====================================================
% Panel C: dyn2-EGFP at random positions
%=====================================================
track1_rand = getRandSlaveSignal(dataDC3(1), track1);

figure(fset.fOpts{:}, 'Position', [12 12 5.75 3.25]); % 8 5.5
ha = axes(fset.axOpts{:}, 'Position', [1.5 1 4 2]); % 6 3.5 
hold on;
plotTrack(dataDC3(1), track1_rand, 2, 'Handle', ha, 'DisplayMode', 'print',...
    'MarkerSizes', ms, 'BackgroundConfidence', norminv(0.95,0,1)*10.5871, ...
    'XLim', [-12 68+12], 'YTick', ya);
set(gca, 'TickLength', fset.TickLength*1.5);
ylabel('Fluo. int. (A.U.)', fset.lfont{:});
% set(gca, 'YTickLabel', [], 'TickLength', fset.TickLength*1.5);
% print('-depsc', '-loose', 'track356_ch2_randSlave.eps');

%%
%=====================================================
% Panel B: Background distribution
%=====================================================
bgDistr = getBackgroundFitDistr(dataDC3(1), 2);
dst = [bgDistr{:}];
lftDataDC3 = getLifetimeData(dataDC3, 'Overwrite', false, 'Cutoff_f', 5,...
    'Scale', true, 'RemoveOutliers', false, 'Mask', true);

% dynamin amplitudes, all tracks
A = arrayfun(@(i) i.A(:,:,2), lftDataDC3, 'unif', 0);
A = vertcat(A{:});

% noise s.d. level
S = arrayfun(@(i) i.sigma_r(:,:,2), lftDataDC3, 'unif', 0);
S = vertcat(S{:});
kLevel = norminv(1-0.05/2.0, 0, 1);

% gaps
G = vertcat(lftDataDC3.gapMat_Ia);

lftAll = vertcat(lftDataDC3.lifetime_s);
tlAll = vertcat(lftDataDC3.trackLengths);
signif = vertcat(lftDataDC3.significantMaster);

A2 = A(signif(:,2)==1,:);
A0 = A(signif(:,2)==0,:);
G2 = G(signif(:,2)==1,:);

lftV = lftAll(signif(:,2)==1);
tlength = tlAll(signif(:,2)==1);

da = 2;
av = -50:da:200;
fset = loadFigureSettings('print');


figure(fset.fOpts{:}, 'Position', [10 10 9.5 5.5]);
axes(fset.axOpts{:});
hold on;
ni = hist(dst, av)/numel(dst)/da;
h = zeros(3,1);
h(1) = stairsXT(av, ni, 'EdgeColor', 'k');
pct = prctile(dst, 95);
plot([pct pct], [0 0.1], 'k--', 'Color', 0*[1 1 1], 'LineWidth', 1);

td = kLevel*nanmean(S(:));
plot([td td], [0 0.1], 'g--', 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 1);

ni = hist(A0(~isnan(A0(:))), av);
ni = ni/sum(ni)/da;
h(3) = stairsXT(av, ni, 'EdgeColor', hsv2rgb([0 1 0.9]));

ni = hist(A2(~isnan(A2(:))), av);
ni = ni/sum(ni)/da;
h(2) = stairsXT(av, ni, 'EdgeColor', hsv2rgb([1/3 1 0.9]));
hl = legend(h, ' enDyn-EGFP background',...
    ' - enDyn2-EGFP ', ' +enDyn2-EGFP positions');
set(hl, 'Box', 'off', 'Position', [5.5 4 2 1.2]);  
set(gca, 'XLim', [-35 100], 'XTick', -50:25:150, 'YLim', [0 0.055]);
xlabel('Fluo. intensity (A.U.)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% print('-depsc2', '-loose', 'dynSignifDistributions.eps');


%%
%=====================================================
% Panel D: Intensity cohorts
%=====================================================
[cohorts, res] = plotIntensityCohorts(dataDC3, [1 2], 'ShowBackground', false,...
    'DisplayMode', 'print', 'ScaleSlaveChannel', false, 'ScalingFactor', [5 1],...
    'ShowLegend', false, 'ShowPct', false, 'YTick', 0:20:80,...
    'ExcludeVisitors', true, 'SlaveName', {'enDyn2-EGFP'});
% print('-depsc2', '-loose', 'dynCohortsDC3.eps');

%%
%=========================================================================================
% Panel E: Lifetime analysis
%=========================================================================================
lftResDC3 = runLifetimeAnalysis(dataDC3, 'MaxIntensityThreshold', 15.59,...
    'RemoveOutliers', false, 'LifetimeData', 'lifetimeData.mat', 'Cutoff_f', 5,...
    'DisplayMode', 'print', 'Display', 'off');

plotLifetimes(lftResDC3, 'ShowCargoDependent', true, 'SlaveNames', 'enDyn2-EGFP',...
    'DisplayMode', 'print', 'YTick', 0:0.005:0.025);

