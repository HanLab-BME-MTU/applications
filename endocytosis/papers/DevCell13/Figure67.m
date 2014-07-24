%=========================================================================================
% Figures 6 & 7a,b
%=========================================================================================
% This script generates the panels for Figure 6 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f6path = '/Users/aguet/Documents/Papers/DevCell13/Figure 6/';

%%
% root = '/home/fa48/lccb-endocytosis/Costin_fs003/RPE_LCa_GFP/';
root = '/home/fa48/files/LCCB/endocytosis/Costin_fs003/RPE_LCa_GFP/';

%-----------------------------------------------------------------------------------------
% Load data
%-----------------------------------------------------------------------------------------
% Control: non-specific siRNA (3 experiments);
dpath = [root 'AP_2_BFP_stables/BFP_only_HIGH_siC7/'];
dataBFP_ctrl = loadConditionData(dpath, {'greenCLC'}, {'EGFP'});

% KD #1: alpha-adaptin KD w/o BFP
dpath = [root 'siRNA_only/sia-adaptin/061311_RPE_Cla_GFp_sialadaptin/'];
data_siRNA_a_adaptin = loadConditionData(dpath, {'greenCLC'}, {'EGFP'});

% KD #2: alpha-adaptin KD + BFP
dpath = [root 'AP_2_BFP_stables/BFP_only_HIGH_si-alpha-ad/092211_RPELCaGFP_BFP_only_si-alpha/'];
dataBFP_KD1 = loadConditionData(dpath, {'greenCLC'}, {'EGFP'});
dpath = [root 'AP_2_BFP_stables/BFP_only_HIGH_si-alpha-ad/100311_RPELCaGFP_BFP_only_si-alpha/'];
dataBFP_KD2 = loadConditionData(dpath, {''}, {'EGFP'});
dataBFP_KD = [dataBFP_KD1 dataBFP_KD2];

% Full-length rescue
dpath = [root 'AP_2_BFP_stables/FL-alpha-ad_HIGH_si-alpha_ad/'];
dataBFP_rescue = loadConditionData(dpath, {'greenCLC'}, {'EGFP'});

% Earless rescue
dpath = [root 'AP_2_BFP_stables/dEar-alpha-ad_HIGH_si-alpha_ad/'];
dataBFP_earless = loadConditionData(dpath, {'greenCLC'}, {'EGFP'});

% Discard obvious outliers
dataBFP_ctrl(9) = [];
dataBFP_earless([24 25]) = [];
dataBFP_rescue([26 27]) = [];

%%
% Run analysis
dataAll = [dataBFP_ctrl data_siRNA_a_adaptin dataBFP_KD dataBFP_rescue dataBFP_earless];
% runDetection(dataAll, 'Overwrite', true);
% runTracking(dataAll, loadTrackSettings(), 'Overwrite', true);
% runTrackProcessing(dataAll, 'Overwrite', true);

%%
%===========================================================
% Panel 6C
%===========================================================
M = xlsread([f6path 'Alpha adaptin ear cells -Tfn-uptake summary -08-15-12.xlsx'], 'D13:X16');
M(isnan(M)) = [];
M = reshape(M, [4 16])';
xa = [0 2.5 5 10];
mu = reshape(M(:,1), [4 4]);
sigma = reshape(M(:,4), [4 4]);
legendText = {' Control (NR siRNA)', ' \alpha-ad. siRNA (no rescue)', ' \alpha-ad. siRNA + FL \alpha-ad.', ' \alpha-ad. siRNA + \DeltaAD \alpha-ad.'};

figure(fset.fOpts{:}, 'Position', [4 4 6 6]);
axes(fset.axOpts{:}, 'Position', [1.5 1.5 4 4], 'TickLength', fset.TickLength*6/4);
hold on;

b = 1;
s = 1;
cv = hsv2rgb([0 0 0;
              0.99 1 b;
              0.3 1 b;
              0.55 1 b]);

for k = 1:4
    he = errorbar(xa, mu(:,k), sigma(:,k), 'Color', cv(k,:), 'LineWidth', 1, 'LineStyle', 'none', 'HandleVisibility', 'off');
    setErrorbarStyle(he, 0.15);
    plot(xa, mu(:,k), '.-', 'Color', cv(k,:), 'LineWidth', 1, 'MarkerSize', 10);

end
axis([0 10.2 0 1.2*1.02]);
set(gca, 'XTick', 0:2.5:10, 'YTickLabel', ['0' arrayfun(@(i) num2str(i, '%.1f'), 0.2:0.2:1.2, 'UniformOutput', false)]);
xlabel('Time of internalization (min)', fset.lfont{:});
ylabel('');
hl = legend(legendText, 'Location', 'NorthEast');
set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [1.6 4.4 1.8 1.2]);

if printFigures
    print('-depsc2', '-loose', [f6path 'TfnInternalization.eps']); %#ok<*UNRCH>
end
%%
%===========================================================
% Panel 6D: fixed cell data comparison
%===========================================================
% sb = 5e-6/(dataBFP_ctrl(1).pixelSize/dataBFP_ctrl(1).M);
% 
% fi = 10;
% frame = double(imread(dataBFP_earless(20).framePaths{1}{fi}));
% [ny nx] = size(frame);
% figure('Position', [50 50 nx ny],...
%     'InvertHardcopy', 'off', 'PaperUnits', 'Points', 'PaperSize', [nx ny],...
%     'PaperPosition', [0 0 nx ny], 'PaperPositionMode', 'auto',...
%     'DefaultLineLineSmoothing','on', 'DefaultPatchLineSmoothing','on');
% axes('Position', [0 0 1 1]);
% imagesc(frame);
% colormap(gray(256));
% plotScaleBar(sb, 'Location', 'SouthEast');
% axis off;
% print('-depsc2', '-loose', 'frame_earless.eps');


sb = 5e-6/(6.7e-6/100);

fpath = '/Users/aguet/Documents/MATLAB/endocytosis/NCBpaper/Figure 6/GFP-LCa_RPE_earCells/';
img_ctrl = double(imread([fpath 'control/TIRF 488 CA-1.tif'])); % 1 or 3
dRange = [min(img_ctrl(:)) max(img_ctrl(:))];
dRange(1) = 1.1*dRange(1);
% figure; imagesc(img_ctrl); axis image; colormap(gray(256)); colorbar;

xi = 1:320;
nx = numel(xi);
figure('Position', [200 200 nx nx], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(img_ctrl(430+xi, 680+xi));
colormap(gray(256));
axis image off;
caxis(dRange);
plotScaleBar(sb, 'Location', 'SouthWest');
if printFigures
    print('-depsc', '-loose', [f6path 'fixedCell_ctrl.eps']);
end

img_earless = double(imread([fpath 'dEar-a-adaptin/TIRF 488 CA-3.tif'])); % 3
% figure; imagesc(img_earless); axis image; colormap(gray(256)); colorbar;
figure('Position', [200 200 nx nx], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(img_earless(380+xi,500+xi));
colormap(gray(256));
axis image off;
caxis(dRange);
if printFigures
    print('-depsc2', '-loose', [f6path 'fixedCell_earless.eps']);
end

img_FL = double(imread([fpath 'FL-a-adaptin/TIRF 488 CA-1.tif'])); % 1
% figure; imagesc(img_FL); axis image; colormap(gray(256)); colorbar;
figure('Position', [200 200 nx nx], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(img_FL(375+xi,605+xi));
colormap(gray(256));
axis image off;
caxis(dRange);
if printFigures
    print('-depsc2', '-loose', [f6path 'fixedCell_FLrescue.eps']);
end

img_KD = double(imread([fpath 'siRNA alpha-adaptin only/TIRF 488 CA-1.tif'])); % 1
% figure; imagesc(img_KD); axis image; colormap(gray(256)); colorbar;
figure('Position', [200 200 nx nx], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
axes('Position', [0 0 1 1]);
frame = img_KD(120+xi,780+xi);
% imagesc(img_KD(370+xi,655+xi));
imagesc(frame);
colormap(gray(256));
axis image off;
tmp = img_KD(120+xi,780+xi);
dRangeX = [0.05 0.15].*diff(dRange)+dRange(1);

hold on;
dx = 15;
rectangle('Position', [dx dx nx-2*dx nx-2*dx], 'EdgeColor', 'w', 'LineStyle', '--', 'LineWidth', 1.5)

caxis(dRange);
print('-depsc', '-loose', [f6path 'fixedCell_KD.eps']);

hold off;
imagesc(frame(dx:nx-dx,dx:nx-dx));
plotScaleBar(sb, 'Location', 'SouthWest');
axis image off;
caxis(dRangeX);
if printFigures
    print('-depsc', '-loose', [f6path 'fixedCell_KD_highContrast.eps']);
end


figure(fset.fOpts{:}, 'Color', 'w');
axes(fset.axOpts{:}, 'Position', [1.5 1.5 0.01 2], 'XColor', 'w', 'YColor', 'w');
colormap(gray(256));
hc = colorbar;
pos = get(hc, 'Position');
pos(3) = pos(3)*1.5;
set(hc, 'Position', pos, 'YTick', []);
if printFigures
    print('-depsc2', [f6path 'fixedCell_colorbar.eps']);
end

%%
%===========================================================
% Panel 6E: initiation density
%===========================================================
% lifetime analysis on control to determine threshold
% runLifetimeAnalysis(data_siRNA_a_adaptin);
%runLifetimeAnalysis(dataBFP_ctrl([1:7 9:end]), 'Display', 'off');

opts = {'MaxIntensityThreshold', T, 'Display', 'off'};
lftRes_siRNA_a_adaptin = runLifetimeAnalysis(data_siRNA_a_adaptin, opts{:});
lftRes_ctrl = runLifetimeAnalysis(dataBFP_ctrl([1:7 9:end]), opts{:});
lftRes_KD = runLifetimeAnalysis(dataBFP_KD([1:3 5:end]), opts{:});
lftRes_rescue = runLifetimeAnalysis(dataBFP_rescue([1:12 14:end]), opts{:});
lftRes_earless = runLifetimeAnalysis(dataBFP_earless([1:20 23:end]), opts{:});


fset = loadFigureSettings('print');
cv = hsv2rgb([0 0 0;
              0.99 1 1;
              0.11 1 1;
              0.3 1 1;
              0.55 1 1]);
cf = rgb2hsv(cv);

cf(:,2) = 0.2;
cf(1,:) = [0 0 0.4];
cf = hsv2rgb(cf);
xlabels = {'Control ', '\alpha-ad. siRNA ', '\alpha-ad. siRNA + BFP ',...
    '\alpha-ad. siRNA + FL \alpha-ad.', '\alpha-ad. siRNA + \DeltaAD \alpha-ad.'};

ha = plotInitiationDensity({lftRes_ctrl, lftRes_siRNA_a_adaptin, lftRes_KD, lftRes_rescue, lftRes_earless},...
    xlabels, cv, cf, {[1 4; 1 5], [1 4; 4 5]});
set(ha(2), 'YTickLabel', {'0', [], '0.1', [], '0.2', [], '0.3', []});
if printFigures
    print('-depsc2', '-loose', 'InitiationDensityComparison_2panels.eps');
end
%%
% Calculate p-values

M_Ia = {lftRes_ctrl.initDensityIa(:,1),...
    lftRes_siRNA_a_adaptin.initDensityIa(:,1),...
    lftRes_KD.initDensityIa(:,1),...
    lftRes_rescue.initDensityIa(:,1),...
    lftRes_earless.initDensityIa(:,1)};

M_above = {lftRes_ctrl.initDensityCCP(:,1),...
     lftRes_siRNA_a_adaptin.initDensityCCP(:,1),...
     lftRes_KD.initDensityCCP(:,1),...
     lftRes_rescue.initDensityCCP(:,1),...
     lftRes_earless.initDensityCCP(:,1),};

%[hval pval] = kstest2(lftRes_rescue.initDensity_above(:,1), lftRes_earless.initDensity_above(:,1)) 
%[hval pval] = kstest2(lftRes_ctrl.initDensity_above(:,1), lftRes_rescue.initDensity_above(:,1)) 

[~,pvalP] = permTest(lftRes_ctrl.initDensityIa(:,1), lftRes_earless.initDensityIa(:,1), 'CmpFunction', @median);
[pvalR,~] = ranksum(lftRes_ctrl.initDensityIa(:,1), lftRes_earless.initDensityIa(:,1));
fprintf('Median density, all detections:  %.3f (control) vs. %.3f (-ear), p = %.3f (rank-sum) p = %.3f (permutation)\n',...
  prctile(M_Ia{1},50), prctile(M_Ia{5}, 50), pvalR, pvalP);

[~,pvalP] = permTest(lftRes_ctrl.initDensityIa(:,1), lftRes_rescue.initDensityIa(:,1), 'CmpFunction', @median);
[pvalR,~] = ranksum(lftRes_ctrl.initDensityIa(:,1), lftRes_rescue.initDensityIa(:,1));
fprintf('Median density, all detections:  %.3f (control) vs. %.3f (rescue), p = %.3f (rank-sum) p = %.3f (permutation)\n',...
  prctile(M_Ia{1},50), prctile(M_Ia{4}, 50), pvalR, pvalP);

% CCPs
[~,pvalP] = permTest(lftRes_ctrl.initDensityCCP(:,1), lftRes_rescue.initDensityCCP(:,1), 'CmpFunction', @median);
[pvalR,~] = ranksum(lftRes_ctrl.initDensityCCP(:,1), lftRes_rescue.initDensityCCP(:,1));
fprintf('Median density, >T detections: %.3f (control) vs. %.3f (rescue), p = %.3f (rank-sum) p = %.3f (permutation)\n',...
  prctile(M_above{1},50), prctile(M_above{4}, 50), pvalR, pvalP);

[~,pvalP] = permTest(lftRes_ctrl.initDensityCCP(:,1), lftRes_earless.initDensityCCP(:,1), 'CmpFunction', @median);
[pvalR,~] = ranksum(lftRes_ctrl.initDensityCCP(:,1), lftRes_earless.initDensityCCP(:,1));
fprintf('Median density, >T detections: %.3f (control) vs. %.3f (-ear), p = %.3e (rank-sum) p = %.3e (permutation)\n',...
    prctile(M_above{1},50), prctile(M_above{5},50), pvalR, pvalP);

%%
%===========================================================
% Panels 7 A,B
%===========================================================
opts = {'MaxIntensityThreshold', 100, 'Display', 'off'};
[lftRes_siRNA_a_adaptin, lftRes2_siRNA_a_adaptin] = runLifetimeAnalysis(data_siRNA_a_adaptin, opts{:});
[lftRes_ctrl, lftRes2_ctrl] = runLifetimeAnalysis(dataBFP_ctrl([1:7 9:end]), opts{:}); %remove 8?
[lftRes_KD, lftRes2_KD] = runLifetimeAnalysis(dataBFP_KD([1:3 5:end]), opts{:});
[lftRes_rescue, lftRes2_rescue] = runLifetimeAnalysis(dataBFP_rescue([1:12 14:end]), opts{:});
[lftRes_earless, lftRes2_earless] = runLifetimeAnalysis(dataBFP_earless([1:20 22:end]), opts{:});
%save lftResDeltaAD lftRes_ctrl lftRes_rescue lftRes_earless lftRes_siRNA_a_adaptin lftRes_KD lftRes2_ctrl lftRes2_rescue lftRes2_earless lftRes2_siRNA_a_adaptin lftRes2_KD

lftRes = {lftRes_ctrl, lftRes_rescue, lftRes_earless, lftRes_siRNA_a_adaptin, lftRes_KD};
leg = {' Control', ' \alpha-ad. siRNA + FL \alpha-ad.',...
    ' \alpha-ad. siRNA + \DeltaAD \alpha-ad.', ' \alpha-ad. siRNA',...
    ' \alpha-ad. siRNA + BFP'};
plotLifetimeComparison(lftRes, leg, 'Frequency', 'relative', 'PlotOrder', [4 5 1 2 3], 'PlotAll', true);

% print('-depsc2', '-loose', 'lftFull_alpha-ad_siRNA.eps');
% print('-depsc2', '-loose', 'lftAbove_alpha-ad_siRNA.eps');

% [H, pValue] = permKSTestMeanEDF(cumsum(lftRes_ctrl.lftHist_Ia,2), cumsum(lftRes_rescue.lftHist_Ia,2))
% [H, pValue] = permKSTestMeanEDF(cumsum(lftRes_ctrl.lftHist_Ia,2), cumsum(lftRes_earless.lftHist_Ia,2))

%%
%=============================================================
% Panel 7G: % persistent CCPs
%=============================================================
opts = {'ReturnValidOnly', false, 'Scale', true, 'Cutoff_f', 5, 'ExcludeVisitors', true};
lftData_ctrl = getLifetimeData(dataBFP_ctrl([1:7 9:end]), opts{:});
lftData_earless = getLifetimeData(dataBFP_earless([1:20 23:end]), opts{:});
catCtrl = {lftData_ctrl.catIdx_all};
catDelta = {lftData_earless.catIdx_all};

T = 100;
idxCtrl = arrayfun(@(i) nanmax(i.A,[],2)>=T, lftData_ctrl, 'unif', 0);
idxDelta = arrayfun(@(i) nanmax(i.A,[],2)>=T, lftData_earless, 'unif', 0);
catCtrlCCP = cellfun(@(x,y) x(y), catCtrl, idxCtrl, 'unif', 0);
catDeltaCCP = cellfun(@(x,y) x(y), catDelta, idxDelta, 'unif', 0);

v = catCtrlCCP;
v = arrayfun(@(i) hist(v{i}, 1:8)/numel(v{i}), 1:numel(v), 'unif', 0);
Vctrl = vertcat(v{:});
v = catDeltaCCP;
v = arrayfun(@(i) hist(v{i}, 1:8)/numel(v{i}), 1:numel(v), 'unif', 0);
Vdelta = vertcat(v{:});

[~,pvalClus] = permTest(Vctrl(:,5), Vdelta(:,5), 'CmpFunction', @mean);
[~,pvalPers] = permTest(Vctrl(:,4), Vdelta(:,4), 'CmpFunction', @mean);
fprintf('Persistent CCPs: %.3f, Clusters: %.3f\n', pvalPers, pvalClus);
% [hval,pval] = ttest2(Vctrl(:,5), Vdelta(:,5))


fset = loadFigureSettings('print');
figure(fset.fOpts{:}, 'Position', [15 5 6 5.5]);
axes(fset.axOpts{:}, 'Position', [1.8 1.5 2 3.5], 'TickLength', fset.TickLength/3.5*6);
hold on;

cv = hsv2rgb([0 0 0;
              0.99 1 1;              
              0.11 1 1;
              0.3 1 1;
              0.55 1 1]);
cf = rgb2hsv(cv);
cf(:,2) = 0.2;
cf(1,:) = [0 0 0.6];
cf = hsv2rgb(cf);

xlabels = {'Control ', '\DeltaAD \alpha-ad.'};

boxplot2({Vctrl(:,4)*100, Vdelta(:,4)*100}, 'BarWidth', 0.6, 'GroupDistance', 0.8,...
    'XTickLabel', xlabels, 'YLim', [0 4.5],...
    'FaceColor', cf([1 5],:), 'EdgeColor', cv([1 5],:),...
    'LineWidth', 1, 'Annotations', [1 2], 'DetectOutliers', false);


ylabel('% persistent CCPs', fset.lfont{:});
formatTickLabels(gca);
if printFigures
    print('-depsc2', '-loose', [f6path 'trackPct.eps']);
end

