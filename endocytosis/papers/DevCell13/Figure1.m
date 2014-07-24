%=========================================================================================
% Figure 1
%=========================================================================================
% This script generates the panels for Figure 1 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f1path = '/Users/aguet/Documents/Papers/DevCell13/Figure 1 - Detection method/';

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [180 180];

% path to data folder
% dpath = 'lccb-endocytosis/Marcel_fs003/110301_BSC1_rat_monkey/BSC1_LCa_ratox/';
dpath = 'Documents/MATLAB/Data/110301_BSC1_rat_monkey/BSC1_LCa_ratox/';

%%
% Data for top row of Panel A
dataOX = loadConditionData(dpath, {''}, {'EGFP'}, 'Parameters', [1.49 100 6.7e-6]);
k = 2; % Cell2!_1s

frameIdx = 50;
load([dataOX(k).source 'Detection/detection_v2.mat']);
frameInfoWT = load([dataOX(k).source 'Detection/detection_v1.mat']);
frameInfoWT = frameInfoWT.frameInfo;

%-------------------------------------
% I. Raw frame (full, with ROI)
%-------------------------------------
frameOX = double(imread(dataOX(k).framePaths{1}{frameIdx})); % EGFP-CLCa
[ny,nx] = size(frameOX);

% selected ROI
xi = 295+(0:160);
yi = 100+(0:160);

% dynamic range of ROI
limsOX = [xi(1)-0.5 xi(end)+0.5 yi(1)-0.5 yi(end)+0.5];
tmp = frameOX(yi,xi);
dRangeOX = [min(tmp(:)) max(tmp(:))];

w = 0.18; % rescale to fit screen
pos0 = [pos(1:2) w*nx w*ny];

figure('Position', pos0, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameOX);
colormap(gray(256));
hold on;
rectangle('Position', [xi(1) yi(1) xi(end)-xi(1) yi(end)-yi(1)], 'EdgeColor', [1 1 1], 'LineWidth', 1);
plotScaleBar(5e-6/(dataOX(k).pixelSize/dataOX(k).M), 'Location', 'SouthEast');
axis off;
caxis(dRangeOX);
if printFigures
    print('-depsc', '-loose', [f1path 'rawOX_full.eps']); %#ok<*UNRCH>
end

% Inverted contrast:
% imagesc(invertContrast(frameOX,dRangeOX));
% rectangle('Position', [xi(1) yi(1) xi(end)-xi(1) yi(end)-yi(1)], 'EdgeColor', [0 0 0], 'LineWidth', 1);
% plotScaleBar(5e-6/(dataOX(k).pixelSize/dataOX(k).M), 'Location', 'SouthEast', 'Color', [0 0 0]);


%-------------------------------------------------------------
% Raw data panels, individual & merged
%-------------------------------------------------------------
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameOX);
colormap(gray(256));
axis(limsOX);
axis off;
caxis(dRangeOX);
plotScaleBar(2e-6/(dataOX(k).pixelSize/dataOX(k).M), 'Location', 'SouthEast');
if printFigures
    print('-depsc', '-loose', [f1path 'rawOX.eps']);
end
% Inverted contrast:
% imagesc(invertContrast(frameOX, dRangeOX));
% plotScaleBar(2e-6/(dataOX(k).pixelSize/dataOX(k).M), 'Location', 'SouthEast', 'Color', [0 0 0]);
% print('-depsc', '-loose', [f1path 'rawOX_inv.eps']);

%%
%-------------------------------------------------------------
% Raw data panels, with detection
%-------------------------------------------------------------
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameOX);
colormap(gray(256));
hold on;
plot(frameInfo(frameIdx).x(1,:), frameInfo(frameIdx).y(1,:), 'go', 'MarkerSize', 8, 'LineWidth', 1);
axis(limsOX);
axis off;
caxis(dRangeOX);
if printFigures
    print('-depsc2', '-loose', [f1path 'detOX.eps']);
end
% Inverted contrast:
% imagesc(invertContrast(frameOX, dRangeOX));
% plot(frameInfo(frameIdx).x(1,:), frameInfo(frameIdx).y(1,:), 'o', 'Color', [0 0.8 0], 'MarkerSize', 8, 'LineWidth', 1);
% print('-depsc2', '-loose', [f1path 'detOX_inv.eps']);

%-------------------------------------------------------------
% Detection mask (significant pixels)
%-------------------------------------------------------------
% figure('Position', pos, 'PaperPositionMode', 'auto');
% axes('Position', [0 0 1 1]);
% mask = 255*(0~=double(imread(dataOX(k).maskPaths{frameIdx})));
% mask = rgbOverlay(frameOX, mask, [1 0 0], dRangeOX);
% imagesc(mask);
% colormap(gray(256));
% hold on;
% axis(limsOX);
% axis off;
% print('-depsc2', '-loose', [f1path 'detOX_mask.eps']);

%-------------------------------------------------------------
% Raw data panels, with tracks
%-------------------------------------------------------------
% tracks  = loadTracks(dataOX(k), 'Cutoff_f', 5, 'Category', 'all',...
%     'FileName', 'ProcessedTracks_gap=3_radius=3,6_a=2.5.mat');
% 
% figure('Position', pos, 'PaperPositionMode', 'auto');
% ha = axes('Position', [0 0 1 1]);
% plotFrame(dataOX(k), tracks, frameIdx, 1, 'Handle', ha, 'ShowGaps', false, 'DisplayType', 'random');
% axis(limsOX);
% axis off;
% caxis(dRangeOX);
% print('-depsc2', '-loose', [f1path 'detOX_tracks.eps']);

%%
%-------------------------------------------------------------
% Compare old detection (wavelets)
%-------------------------------------------------------------
% runWaveletDetection(data, varargin)
% [frameInfoWavelets, maskWavelets] = detectSpotsWT(frameOX, 4, 5, 1);

%-----------------------------------------------------------------------------------------
% Wavelet-based detection (with post-processing, as used in Loerke et al. 2009)
%-----------------------------------------------------------------------------------------
maskWT = double(imread([dataOX(k).source 'Detection/WaveletMasks/' 'wmask_' num2str(frameIdx, '%.3d') '.tif']));
maskWT = 255*(maskWT~=0);

% Using model-based approach as reference, identify
% 1) false positives (empty index from KD query in radius 2)
R = 2;
fpIdx = KDTreeBallQuery([frameInfo(frameIdx).x' frameInfo(frameIdx).y'], [frameInfoWT(frameIdx).xav(:,1) frameInfoWT(frameIdx).yav(:,1)], R);
fpIdx = cellfun(@(i) isempty(i), fpIdx);
% 2) false negatives
fnIdx = KDTreeBallQuery([frameInfoWT(frameIdx).xav(:,1) frameInfoWT(frameIdx).yav(:,1)], [frameInfo(frameIdx).x' frameInfo(frameIdx).y'], R);
fnIdx = cellfun(@(i) isempty(i), fnIdx);

frameRGB = rgbOverlay(frameOX, maskWT, [1 0 0], dRangeOX);

figure('Position', pos, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
axes('Position', [0 0 1 1]);
imagesc(frameRGB);
hold on;
plot(frameInfoWT(frameIdx).xav(fpIdx,1), frameInfoWT(frameIdx).yav(fpIdx,1), 'o', 'Color', hsv2rgb([0.1 1 1]), 'MarkerSize', 8, 'LineWidth', 1);
plot(frameInfo(frameIdx).x(fnIdx), frameInfo(frameIdx).y(fnIdx), 'o', 'Color', hsv2rgb([1/3 0 1]), 'MarkerSize', 8, 'LineWidth', 1);
axis(limsOX);
axis off;
if printFigures
    print('-depsc', '-loose', [f1path 'detOX_WT_L1.eps']);
end

%-------------------------------------------------------------
% Without post-processing (not used in figure)
%-------------------------------------------------------------
% [frameInfoWT, maskWavelets] = main283AUTO_standalone(frameOX, 0);
% 
% % Using model-based approach as reference, identify
% % 1) false positives (empty index from KD query in radius 2)
% R = 2;
% fpIdx = KDTreeBallQuery([frameInfo(frameIdx).x' frameInfo(frameIdx).y'], [frameInfoWT.xav(:,1) frameInfoWT.yav(:,1)], R);
% fpIdx = cellfun(@(i) isempty(i), fpIdx);
% % 2) false negatives
% fnIdx = KDTreeBallQuery([frameInfoWT.xav(:,1) frameInfoWT.yav(:,1)], [frameInfo(frameIdx).x' frameInfo(frameIdx).y'], R);
% fnIdx = cellfun(@(i) isempty(i), fnIdx);
% 
% maskWT = 255*(maskWavelets~=0);
% frameRGB = rgbOverlay(frameOX, maskWT, [1 0 0], dRangeOX);
% 
% figure('Position', pos, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
% axes('Position', [0 0 1 1]);
% imagesc(frameRGB);
% hold on;
% plot(frameInfoWT.xav(fpIdx,1), frameInfoWT.yav(fpIdx,1), 'o', 'Color', hsv2rgb([0.1 1 1]), 'MarkerSize', 8, 'LineWidth', 1);
% plot(frameInfo(frameIdx).x(fnIdx), frameInfo(frameIdx).y(fnIdx), 'o', 'Color', hsv2rgb([1/3 0 1]), 'MarkerSize', 8, 'LineWidth', 1);
% axis(limsOX);
% axis off;
% if printFigures
%     print('-depsc', '-loose', [f1path 'detOX_WT_L0.eps']);
% end

% histograms of the estimated intensities
% figure; hist(frameInfo(frameIdx).A, linspace(min(frameInfo(frameIdx).A), max(frameInfo(frameIdx).A), 30));
% Aavg = frameInfoWavelets.totalInt ./ frameInfoWavelets.area;
% figure; hist(Aavg, linspace(min(Aavg), max(Aavg), 30))

%%
%=========================================================================================
% Data for bottom row of Panel A
dpath = 'lccb-endocytosis/Marcel_fs003/110304_CLTa_ZNF/BSC1_CLTA_RFP/';
dataEN = loadConditionData(dpath, {''}, {'RFP'}, 'Parameters', [1.49 100 6.7e-6]);
k = 4; % Cell4_1s

% frameIdx = 50;
frameEN = double(imread(dataEN(k).framePaths{1}{frameIdx})); % enCLCa-RFP
[ny,nx] = size(frameEN);

tmp = load([dataEN(k).source 'Detection/detection_v2.mat']);
frameInfoEN = tmp.frameInfo;
frameInfoWT = load([dataEN(k).source 'Detection/detection_v1.mat']);
frameInfoWT = frameInfoWT.frameInfo;

%-------------------------------------
% I. Raw frame (full, with ROI)
%-------------------------------------

% selected ROI
xi = 410+(0:160);
yi = 430+(0:160);

% dynamic range of ROI
limsEN = [xi(1)-0.5 xi(end)+0.5 yi(1)-0.5 yi(end)+0.5];
tmp = frameEN(yi,xi);
dRangeEN = [min(tmp(:)) max(tmp(:))];

% adjust height relative to OX frame
w = 0.18;
pos0 = [pos(1:2) w*nx w*ny];

figure('Position', pos0, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameEN);
colormap(gray(256));
hold on;
rectangle('Position', [xi(1) yi(1) xi(end)-xi(1) yi(end)-yi(1)], 'EdgeColor', [1 1 1], 'LineWidth', 1);
plotScaleBar(5e-6/(dataEN(k).pixelSize/dataEN(k).M), 'Location', 'SouthEast');
axis off;
caxis(dRangeEN);
if printFigures
    print('-depsc', '-loose', [f1path 'rawEN_full.eps']);
end

%-------------------------------------------------------------
% Raw data panels, individual & merged
%-------------------------------------------------------------
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameEN);
colormap(gray(256));
axis(limsEN);
axis off;
caxis(dRangeEN);
plotScaleBar(2e-6/(dataEN(k).pixelSize/dataEN(k).M), 'Location', 'SouthEast');
if printFigures
    print('-depsc', '-loose', [f1path 'rawEN.eps']);
end

%-------------------------------------------------------------
% Raw data panels, with detection
%-------------------------------------------------------------
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameEN);
colormap(gray(256));
hold on;
plot(frameInfoEN(frameIdx).x(1,:), frameInfoEN(frameIdx).y(1,:), 'go', 'MarkerSize', 8, 'LineWidth', 1);
axis(limsEN);
axis off;
caxis(dRangeEN);
if printFigures
    print('-depsc2', '-loose', [f1path 'detEN.eps']);
end

%-----------------------------------------------------------------------------------------
% Wavelet-based detection (with post-processing, as used in Loerke et al. 2009
%-----------------------------------------------------------------------------------------
maskWT = double(imread([dataEN(k).source 'Detection/WaveletMasks/' 'wmask_' num2str(frameIdx, '%.3d') '.tif']));
maskWT = 255*(maskWT~=0);

% Using model-based approach as reference, identify
% 1) false positives (empty index from KD query in radius 2)
R = 2;
fpIdx = KDTreeBallQuery([frameInfoEN(frameIdx).x' frameInfoEN(frameIdx).y'], [frameInfoWT(frameIdx).xav(:,1) frameInfoWT(frameIdx).yav(:,1)], R);
fpIdx = cellfun(@(i) isempty(i), fpIdx);
% 2) false negatives
fnIdx = KDTreeBallQuery([frameInfoWT(frameIdx).xav(:,1) frameInfoWT(frameIdx).yav(:,1)], [frameInfoEN(frameIdx).x' frameInfoEN(frameIdx).y'], R);
fnIdx = cellfun(@(i) isempty(i), fnIdx);

frameRGB = rgbOverlay(frameEN, maskWT, [1 0 0], dRangeEN);

figure('Position', pos, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
axes('Position', [0 0 1 1]);
imagesc(frameRGB);
hold on;
plot(frameInfoWT(frameIdx).xav(fpIdx,1), frameInfoWT(frameIdx).yav(fpIdx,1), 'o', 'Color', hsv2rgb([0.1 1 1]), 'MarkerSize', 8, 'LineWidth', 1);
plot(frameInfoEN(frameIdx).x(fnIdx), frameInfoEN(frameIdx).y(fnIdx), 'o', 'Color', hsv2rgb([1/3 0 1]), 'MarkerSize', 8, 'LineWidth', 1);
axis(limsEN);
axis off;
if printFigures
    print('-depsc', '-loose', [f1path 'detEN_WT_L1.eps']);
end

%-------------------------------------------------------------
% Without post-processing (not used in figure)
%-------------------------------------------------------------
% [frameInfoWT, maskWavelets] = main283AUTO_standalone(frameEN,0);
% frameInfoWT.xcom = frameInfoWT.xav;
% frameInfoWT.ycom = frameInfoWT.yav;
% 
% % Using model-based approach as reference, identify
% % 1) false positives (empty index from KD query in radius 2)
% R = 2;
% fpIdx = KDTreeBallQuery([frameInfoEN(frameIdx).x' frameInfoEN(frameIdx).y'], [frameInfoWT.xcom(:,1) frameInfoWT.ycom(:,1)], R);
% fpIdx = cellfun(@(i) isempty(i), fpIdx);
% % 2) false negatives
% fnIdx = KDTreeBallQuery([frameInfoWT.xcom(:,1) frameInfoWT.ycom(:,1)], [frameInfoEN(frameIdx).x' frameInfoEN(frameIdx).y'], R);
% fnIdx = cellfun(@(i) isempty(i), fnIdx);
% 
% maskWT = 255*(maskWavelets~=0);
% frameRGB = rgbOverlay(frameEN, maskWT, [1 0 0], dRangeEN);
% 
% figure('Position', pos, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
% axes('Position', [0 0 1 1]);
% imagesc(frameRGB);
% hold on;
% plot(frameInfoWT.xcom(fpIdx,1), frameInfoWT.ycom(fpIdx,1), 'o', 'Color', hsv2rgb([0.1 1 1]), 'MarkerSize', 8, 'LineWidth', 1);
% plot(frameInfoEN(frameIdx).x(fnIdx), frameInfoEN(frameIdx).y(fnIdx), 'o', 'Color', hsv2rgb([1/3 0 1]), 'MarkerSize', 8, 'LineWidth', 1);
% axis(limsEN);
% axis off;
% if printFigures
%     print('-depsc', '-loose', [f1path 'detEN_WT_L0.eps']);
% end

%%
%=========================================================================================
% Panel B: Track examples
%=========================================================================================
% runTrackProcessing(dataOX, 'Buffer', [10 5], 'TrackerOutput', 'trackedFeatures_buffer=[10,5].mat',...
%     'FileName', 'ProcessedTracks_buffer=[10,5].mat', 'Overwrite', true);

k = 2; % Cell2!_1s
tracksOX = loadTracks(dataOX(k), 'Category', 'Ia', 'Cutoff_f', 5, 'FileName', 'ProcessedTracks_buffer=[10,5].mat');
tracksOX = tracksOX(40<=[tracksOX.lifetime_s] & [tracksOX.lifetime_s]<=50);
% cmeDataViewer(dataOX(k), tracksOX);
% 372 w/ [5 5] buffers
% itrackOX = tracksOX(387);
% itrackOX = tracksOX(370); % latest version, [5 5] buffers
itrackOX = tracksOX(366); % start 530, length 46


k = 4; % Cell4_1s
tracksEN = loadTracks(dataEN(k), 'Category', 'Ia', 'Cutoff_f', 5, 'FileName', 'ProcessedTracks_buffer=[10,5].mat');
% itrackEN = tracksEN(534); %too short
% itrackEN = tracksEN(514);
% 402
itrackEN = tracksEN(371);


fset = loadFigureSettings('print');
figure;
ha = axes(fset.axOpts{:}, 'Position', [2 2 7 3]);
hold on;
plotTrack(dataOX(k), itrackOX, 'Handle', ha, 'DisplayMode', 'print', 'MarkerSizes', [12 4 1.25]);
set(ha, 'TickLength', fset.TickLength*6/7, 'XTick', -10:10:50, 'YTick', -25:25:300);
axis([-12 itrackOX.t(end)-itrackOX.t(1)+7 -40 225]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
if printFigures
    print('-depsc2', [f1path 'trackExampleOX.eps']);
end

[stack, xa, ya] = getTrackStack(dataOX(2), itrackOX,... % k = 2
    'WindowWidth', 6, 'Reference', 'track');
stack = stack(6:end);
xtrack = itrackOX;
xtrack.startBuffer.t = xtrack.startBuffer.t(6:end);
plotTrackMontage(xtrack, stack, xa, ya,...
    'ShowMarkers', true,...
    'ShowDetection', false, 'FramesPerRow', 20);
if printFigures
    print('-depsc2', '-loose', [f1path 'trackMontageOX.eps']);
end


% ratio of track lengths
pr = (numel(itrackEN.t)+19) / (numel(itrackOX.t)+19);
pr2 = (numel(itrackEN.t)+4) / (numel(itrackOX.t)+4);

figure;
ha = axes(fset.axOpts{:}, 'Position', [2 2 7*pr 3]);
hold on;
plotTrack(dataEN(k), itrackEN, 'Handle', ha, 'DisplayMode', 'print', 'MarkerSizes', [12 4 1.25]);
set(ha, 'TickLength', fset.TickLength*6/(7*pr), 'XTick', -10:10:50, 'YTick', -25:25:300);%, 'Color', 'none');
axis([-12 itrackEN.t(end)-itrackEN.t(1)+7 -40 175]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
if printFigures
    print('-depsc2', [f1path 'trackExampleEN.eps']);
end

[stack, xa, ya] = getTrackStack(dataEN(k), itrackEN,...
    'WindowWidth', 6, 'Reference', 'track');
stack = stack(6:end);
xtrack = itrackEN;
xtrack.startBuffer.t = xtrack.startBuffer.t(6:end);
plotTrackMontage(xtrack, stack, xa, ya,...
    'ShowMarkers', true, 'Width', 600,...
    'ShowDetection', false, 'FramesPerRow', 20);
if printFigures
    print('-depsc2', '-loose', [f1path 'trackMontageEN.eps']);
end
