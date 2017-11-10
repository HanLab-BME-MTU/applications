%=========================================================================================
% Figure 4
%=========================================================================================
% This script generates the panels for Figure 4 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f4path = '/Users/aguet/Documents/Papers/DevCell13/Figure 4/';


root = '/home/fa48/lccb-endocytosis/';
fset = loadFigureSettings('print');

%%
%=====================================================================
% Panel a): Lifetime dist. OX vs EN
%=====================================================================
dpath = [root 'Marcel_fs003/110301_BSC1_rat_monkey/BSC1_LCa_ratox/'];
dataOX = loadConditionData(dpath, {''}, {'EGFP'});

dpath = [root 'Marcel_fs003/110304_CLTa_ZNF/BSC1_CLTA_RFP/'];
dataEN = loadConditionData(dpath, {''}, {'RFP'});

opts = {'Display', 'off', 'ProcessedTracks', 'ProcessedTracks_gap=2_radius=5,10_a=2.5.mat'...
    'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat'};
lftResOX = runLifetimeAnalysis(dataOX, opts{:});
lftResEN = runLifetimeAnalysis(dataEN, 'Display', 'off');

%%
figure(fset.fOpts{:}, 'Name', 'Lifetime histograms (raw)', 'Position', [15 5 8 6]);
axes(fset.axOpts{:});
hold on;
plot(lftResEN.t, mean(vertcat(lftResEN.lftHist_Ia), 1), 'Color', hsv2rgb([0 1 0.9]), 'LineWidth', 1.5);
plot(lftResOX.t, mean(vertcat(lftResOX.lftHist_Ia), 1), 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 1.5);
ya = 0:0.01:0.05;
axis([0 min(120, lftResOX.t(end)) 0 ya(end)]);
set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% print('-depsc2', '-loose', 'LftRaw_ENvsOX.eps');


pOX = [mean(lftResOX.pctCS) mean(lftResOX.pctCCP)];
pOXstd = [std(lftResOX.pctCS) std(lftResOX.pctCCP)];
pOXstd = pOXstd/sum(pOX);
pOX = pOX/sum(pOX);

pEN = [mean(lftResEN.pctCS) mean(lftResEN.pctCCP)];
pENstd = [std(lftResEN.pctCS) std(lftResEN.pctCCP)];
pENstd = pENstd/sum(pEN);
pEN = pEN/sum(pEN);

cla;
hp = zeros(1,2);
hp(2) = plot(lftResEN.t, lftResEN.meanLftHistCCP, 'Color', hsv2rgb([0 1 0.9]), 'LineWidth', 1.5);
hp(1) = plot(lftResOX.t, lftResOX.meanLftHistCCP, 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 1.5);

hl = legend(hp, ' EGFP-CLCa O/X', ' enCLCa-RFP', 'Location', 'NorthEast');
set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [1.5 5.25 1.5 0.75]);

%-------------------------------------------------
% Inset with PSNR distributions
%-------------------------------------------------
ha = axes(fset.axOpts{:}, 'Position', [4.5 3.2 3 1.8], 'TickLength', fset.TickLength*6/3.5);
hold on;
set(gca, 'xscale', 'log');
plotPSNRDistribution(dataEN, ha, 'Color', hsv2rgb([0 1 0.9]), 'DisplayMode', 'print');
plotPSNRDistribution(dataOX, ha, 'DisplayMode', 'print');
set(ha, 'YTick', 0:0.01:0.04, 'YTickLabel', 0:0.01:0.04, 'YLim', [0 0.04]);
if printFigures
    print('-depsc2', '-loose', [f4path 'LftThresholded_ENvsOX_PSNRinset.eps']);
end

%%
%=====================================================================
% Panel b): Missed gaps simulation
%=====================================================================
% 1) Fit gamma distribution to experimental data, or use experimental distr.

% t = lftResOX.t(1:160);
% f = lftResOX.meanLftHist_A(1:160);
% f = f/sum(f);
% [prmVect, a, BIC] = fitGammaMixture(t, f, 'kna', 'FitMode', 'PDF', 'N', 1);
% ti = 0:120;
% g = gammaMixture(ti, prmVect, 'Mode', 'PDF');
% figure;
% hold on;
% plot(t,f);
% plot(ti, g, 'r');
prmVect = [0.0484 2.2885 1];

% Compute gap probabilities, smoothen
w = 20;
resGP = getGapProbability(dataOX, 'WindowWidth', w);
% resGP = getGapProbability(dataEN, 'WindowWidth', w);

fs = resGP.Pstart(2:end);
fe = resGP.Pend(1:end-1);
combP = [fs fe];
combP = 1.1*filterGauss1D(combP, 1, 'replicate');
fs = combP(1:w-1);
fe = combP(w:end);
% xi = 1:numel(combP);
% pps = csaps(xi,combP, 0.2);
% pps = spaps(xi,combP, 0.01);
% tmp = fnval(pps, xi);
% fs = [fs(1) tmp(2:w-1)];
% fe = [tmp(w:end-1) fe(end)];
meanP = mean([fs(end) fe(1)]);

% generate lifetime samples
t = 1:160;
nt = numel(t);
pdf = gammaMixture(t, prmVect, 'Mode', 'PDF');
cdf = gammaMixture(t, prmVect, 'Mode', 'CDF');

nExp = 1e6;
lft = round(interp1(cdf, t, rand(1,nExp)));
lft(isnan(lft)) = [];
nExp = numel(lft);

% for each sample, determine gap PDF
gapPM = zeros(nExp, nt);
for k = 1:nExp
    n = max(0,lft(k)-2*(w-1));
    gapPM(k,1:w-1) = fs;
    if n>=1
        gapPM(k,w+(0:n-1)) = meanP;
    end
    gapPM(k,w-1+n+(1:w-1)) = fe;
end
gapMask = gapPM > rand(size(gapPM));

% Gaps following experimental distribution (# given by distribution)
cutLft = cell(1,nExp);
nGaps = zeros(1,nExp);
parfor k = 1:nExp
    v = gapMask(k,1:lft(k));
    idx = find(v~=0);
    %if numel(idx)>1
        %idx = idx(1);
        %idx = unidrnd(numel(idx),1,1);
    %end
    nGaps(k) = numel(idx);
    cutLft{k} = diff([0 idx lft(k)+1])-1;
end
cutLft = [cutLft{:}];
cutLft(cutLft<1) = [];
% figure; hist(cutLft, 200);
% figure; hist(nGaps, 0:20);

% Single gap per frame, uniform distr.
cutLftS = zeros(nExp,2);
parfor k = 1:nExp
    gapIdx = unidrnd(lft(k),1,1);
    cutLftS(k,:) = [gapIdx-1 lft(k)-gapIdx];
end
cutLftS(cutLftS<1) = [];
% figure; hist(cutLftS(:),200);

%%
w = resGP.w;
dx = 10;
xs = 2:w;
xe = w+dx+ (1:w-1);

cutPDF = hist(cutLft, t);
cutPDF = cutPDF / sum(cutPDF);

cutPDFsingle = hist(cutLftS, t);
cutPDFsingle = cutPDFsingle / sum(cutPDFsingle);

figure(fset.fOpts{:});
axes(fset.axOpts{:});
hold on;
hp(3) = plot(t, cutPDFsingle, 'k', 'LineWidth', 1);
hp(2) = plot(t, cutPDF, '-', 'Color', fset.ceB, 'LineWidth', 1);
hp(1) = plot(t, pdf, 'r', 'LineWidth', 1);
axis([0 120 0 0.05]);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:})

hl = legend(hp, ' True distribution', ' Missed gaps: 1/track', ' Missed gaps: P(gap)', 'Location', 'NorthWest');
set(hl, 'Box', 'on', 'EdgeColor', [1 1 1], fset.ifont{:});
lpos = get(hl, 'Position');
% lpos(3) = 2;
lpos = [1.7 4.5 1 1];
set(hl, 'Position', lpos);

ha = axes(fset.axOpts{:}, 'Position', [5 3.5 2.5 1.5]);
hold on;
plot(xs, fs, '-', 'Color', fset.ceB, 'LineWidth', 1);
plot(xe, fe, '-', 'Color', fset.ceB, 'LineWidth', 1);
plot([xs(end)+1 xe(1)-1], meanP*[1 1], ':', 'Color', fset.ceB, 'LineWidth', 1);
set(gca, 'XTick', [2 10 20 [21 31 39]+dx], 'XTickLabel', {'2' '10' '20' 'n-20' 'n-10' 'n-2'},...
    'YTick', 0:0.1:0.3, fset.ifont{:}, 'TickLength', fset.TickLength*2);
axis([1 xe(end)+1 0 0.35]);
xlabel('Frame #', fset.ifont{:});
% ylabel('P(gap)', fset.ifont{:});
rotateXTickLabels(ha, 'AdjustFigure', false, 'YOffset', 0.05);

print('-depsc', '-loose', 'gapSimulation.eps');




%%
%=============================
% Panel C
%=============================
f4path = '/Users/aguet/Documents/MATLAB/endocytosis/CMEpaper/Figure 4 - Dual channel detection/';

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [180 180];

% ZFN
dpath = '/Users/aguet/Documents/MATLAB/endocytosis/testData/Marcel_fs003/111018_BSC1_dual_channel/BSC1_LCaZNF_EGFPm2/';
dataEN = loadConditionData(dpath, {'green_2s', 'red_2s'}, {'EGFP', 'RFP'});
k = 1; % Cell8_nice_2s

% Control, tdTomato-LCa (transient)
% dpath = '/Users/aguet/Documents/MATLAB/endocytosis/testData/Marcel_fs003/111221_BSC1_muEGFP2/control_dual_channel/';
% dataAP2Master = loadConditionData(dpath, {'green_2s', 'red_2s'}, {'EGFP', 'tdTomato'});
% data = dataAP2Master(3);
frameIdx = 50;
load([dataEN(k).source 'Detection/detection_v2.mat']);


frameMaster = double(imread(dataEN(k).framePaths{1}{frameIdx})); % µ2-EGFP
frameSlave = double(imread(dataEN(k).framePaths{2}{frameIdx})); % tdTomato-CLC
[ny,nx] = size(frameMaster);

% xi = 85+(0:160);
% yi = 305+(0:160);
xi = 85+(0:159);
yi = 305+(0:159);


lims = [xi(1)-0.5 xi(end)+0.5 yi(1)-0.5 yi(end)+0.5];

% dynamic ranges of the cutout regions
tmp = frameMaster(yi,xi);
dRangeM = [min(tmp(:)) max(tmp(:))];
tmp = frameSlave(yi,xi);
dRangeS = [min(tmp(:)) max(tmp(:))];

frameRGB = uint8(cat(3, scaleContrast(frameSlave, dRangeS, [0 255]),...
    scaleContrast(frameMaster, dRangeM, [0 255]), zeros(size(frameMaster))));



% full frame with inset
w = 0.7;
pos0 = [pos(1:2) nx*w ny*w];

figure('Position', pos0, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameMaster);
colormap(gray(256));
hold on;
caxis(dRangeM);
rectangle('Position', [xi(1) yi(1) 151 151], 'EdgeColor', [1 1 1], 'LineWidth', 1.5);
plotScaleBar(5e-6/(dataEN(k).pixelSize/dataEN(k).M), 'Location', 'SouthEast');
axis off;
% print('-depsc', '-loose', [f4path 'rawMaster_full.eps']);

figure('Position', pos0, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
% imagesc(frameSlave);
imagesc(frameRGB);
colormap(gray(256));
hold on;
caxis(dRangeS);
rectangle('Position', [xi(1) yi(1) 151 151], 'EdgeColor', [1 1 1], 'LineWidth', 1.5);
plotScaleBar(5e-6/(dataEN(k).pixelSize/dataEN(k).M), 'Location', 'SouthEast');
axis off;
% print('-depsc', '-loose', [f4path 'rawSlave_full.eps']);
print('-depsc', '-loose', ['test3.eps']);



%-------------------------------------------------------------
% Raw data panels, individual & merged
%-------------------------------------------------------------
% raw frame, ch 1
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameMaster);
colormap(gray(256));
axis(lims, 'off');
caxis(dRangeM);
h = plotScaleBar(2e-6/(dataEN(k).pixelSize/dataEN(k).M), 'Location', 'SouthEast');
% print('-depsc', '-loose', [f4path 'raw_ch1NEW.eps']);
delete(h);
hold on;
plot(frameInfo(frameIdx).x(1,:), frameInfo(frameIdx).y(1,:), 'go', 'MarkerSize', 8, 'LineWidth', 1);
% print('-depsc', '-loose', [f4path 'det_ch1NEW.eps']);


% raw frame, ch 2
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameSlave);
colormap(gray(256));
axis(lims, 'off');
caxis(dRangeS);
% print('-depsc2', '-loose', [f4path 'raw_ch2.eps']);
hold on;
plot(frameInfo(frameIdx).x(2,:), frameInfo(frameIdx).y(2,:), 'ro', 'MarkerSize', 8, 'LineWidth', 1);
% print('-depsc2', '-loose', [f4path 'det_ch2.eps']);


% raw frame, RGB
figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(frameRGB);
axis(lims, 'off');
% print('-depsc2', '-loose', [f4path 'raw_RGB.eps']);
hold on;
plot(frameInfo(frameIdx).x(2,:), frameInfo(frameIdx).y(2,:), 'ro', 'MarkerSize', 8, 'LineWidth', 1);
plot(frameInfo(frameIdx).x(1,:), frameInfo(frameIdx).y(1,:), 'go', 'MarkerSize', 8, 'LineWidth', 1);
% print('-depsc2', '-loose', [f4path 'det_RGB.eps']);


    
%%
%=====================================================================
% Panel D: Lifetime dist. for AP2/LCa dual channel
%=====================================================================
dpath = '/home/fa48/lccb-endocytosis/Marcel_fs003/111018_BSC1_dual_channel/BSC1_LCaZNF_EGFPm2/';
dataAP2Master = loadConditionData(dpath, {'green_2s', 'red_2s'}, {'EGFP', 'RFP'});
dataLCaMaster = loadConditionData(dpath, {'red_2s', 'green_2s'}, {'RFP', 'EGFP'});
% exclude static cell (7) & outlier (10)
dataAP2Master = dataAP2Master([1:6 8:9]);
dataLCaMaster = dataLCaMaster([1:6 8:9]);

opts = {'Cutoff_f', 0, 'RemoveOutliers', false, 'Display', 'off'};
lftResAP2Master0 = runLifetimeAnalysis(dataAP2Master, opts{:}, 'MaxIntensityThreshold', 83);
lftResLCaMaster0 = runLifetimeAnalysis(dataLCaMaster, opts{:}, 'MaxIntensityThreshold', 14);

[map, masterIdx, slaveIdx, unassignedMasterIdx, unassignedSlaveIdx] = assignSlaveTracksToMaster(dataLCaMaster, dataAP2Master,...
    'InputMode', 'data', 'MinOverlap', 1);

opts = {'ExcludeVisitors', false, 'Cutoff_f', 0, 'RemoveOutliers', false, 'Display', 'off'};
lftResAP2Master = runLifetimeAnalysis(dataAP2Master, opts{:}, 'MaxIntensityThreshold', 83, 'SelectIndex', masterIdx);
lftResLCaMaster = runLifetimeAnalysis(dataLCaMaster, opts{:}, 'MaxIntensityThreshold', 14, 'SelectIndex', slaveIdx);

%%
figure(fset.fOpts{:}, 'Name', 'Lifetime histograms (raw)');
axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 3.5/5*3.1]);
hold on;

xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});

ya = 0:0.01:0.03;
axis([0 min(120, lftResAP2Master.t(end)) 0 0.031]);
set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);

pLCa = [mean(lftResLCaMaster.pctCS) mean(lftResLCaMaster.pctCCP)];
pLCastd = [std(lftResLCaMaster.pctCS) std(lftResLCaMaster.pctCCP)];
pLCastd = pLCastd/sum(pLCa);
pLCa = pLCa/sum(pLCa);

pAP2 = [mean(lftResAP2Master.pctCS) mean(lftResAP2Master.pctCCP)];
pAP2std = [std(lftResAP2Master.pctCS) std(lftResAP2Master.pctCCP)];
pAP2std = pAP2std/sum(pAP2);
pAP2 = pAP2/sum(pAP2);
hp = zeros(2,1);
hp(2) = plot(lftResLCaMaster.t, lftResLCaMaster.meanLftHistCCP, 'Color', hsv2rgb([0.07 1 1]), 'LineWidth', 1.5);
hp(1) = plot(lftResAP2Master.t, lftResAP2Master.meanLftHistCCP, 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 1.5);
hl = legend(hp, [' ' char(181) '2-EGFP'], ' enCLCa-RFP', 'Location', 'NorthEast');
set(hl, 'Box', 'off', 'Position', [5 3 1.5 0.8]);
% print('-depsc2', '-loose', 'LftThresholded_ENvsAP2.eps');


%%
%=====================================================================
% Panel E: Mean values (bootstrapped)
%=====================================================================

figure(fset.fOpts{:}, 'Position', [15 5 3.5 6]);
axes(fset.axOpts{:}, 'Position', [1.5 2 1.75 3.5], 'TickLength', fset.TickLength*6/3.5);
hold on;
% 1s data set:
% EGFP-CLCa O/X

clear mu sigma;
[mu(1), sigma(1)] = getLifetimeDistributionMean(lftResOX);
[mu(2), sigma(2)] = getLifetimeDistributionMean(lftResAP2Master);
[mu(3), sigma(3)] = getLifetimeDistributionMean(lftResLCaMaster);

hues = [0 0.06 0.33]'; % OX, EN, AP2
cv = [0 0 0;
      0.33 1 0.8;
      0.07 1 1];
  cv = hsv2rgb(cv);

for k = 1:3
    he = errorbar(k, mu(k), sigma(k), 'Color', cv(k,:,:), 'LineStyle', 'none', 'LineWidth', 1.2);
    setErrorbarStyle(he, 0.125);
    plot(k, mu(k), '.', 'Color', cv(k,:), 'MarkerSize', 10);
end

axis([0.5 3.5 0 60]);
txt = {'EGFP-CLCa O/X', [char(181) '2-EGFP'], 'enCLCa-RFP'};
set(gca, 'YTick', 0:10:60, 'XTick', 1:3, 'XTickLabel', txt);
rotateXTickLabels(gca, 'AdjustFigure', false, 'Angle', 50);
ylabel('Mean lifetime (s)', fset.lfont{:});
% print('-depsc2', '-loose', 'LifetimeMeans.eps');

%%
sum(lftResOX.meanLftHist_A.*lftResOX.t)
sum(lftResEN.meanLftHist_A.*lftResEN.t)

% sum(mean(vertcat(lftResOX.lftHist_Ia), 1).*lftResOX.t)
% sum(mean(vertcat(lftResEN.lftHist_Ia), 1).*lftResEN.t)



%%
%=====================================================================
% Panel f: matched track examples
%=====================================================================
%dpath = '/Users/aguet/Documents/MATLAB/endocytosis/testData/Marcel_fs003/111018_BSC1_dual_channel/BSC1_LCaZNF_EGFPm2/';
dpath = '/home/fa48/lccb-endocytosis/Marcel_fs003/111018_BSC1_dual_channel/BSC1_LCaZNF_EGFPm2/';
dataAP2Master = loadConditionData(dpath, {'green_2s', 'red_2s'}, {'EGFP', 'RFP'});
dataCLaMaster = loadConditionData(dpath, {'red_2s', 'green_2s'}, {'RFP', 'EGFP'});

% select cell10
tracksAP2Master = loadTracks(dataAP2Master(10), 'Category', 'Ia', 'Cutoff_f', 2);
tracksCLaMaster = loadTracks(dataCLaMaster(10), 'Category', 'Ia', 'Cutoff_f', 2);

% idx are slave fragments; if empty: master has no detected slaves
idx = assignTracksToMaster(tracksCLaMaster, tracksAP2Master);

% matchedIdx: master tracks that have slave signal
matchedIdx = find(cellfun(@(i) ~isempty(i), idx)); % defined with cutoff = 4 !!
unMatchedIdx = setdiff(1:numel(tracksAP2Master), matchedIdx);

%%
% trackDisplayGUI(dataAP2Master(1), tracksAP2Master(matchedIdx));

% determine max. intensities
% maxCh1 = max(arrayfun(@(i) max(i.A(1,:)+i.c(1,:)-mean(i.c(1,:))), tracksAP2Master(candIdx)))
% maxCh2 = max(arrayfun(@(i) max(i.A(2,:)+i.c(2,:)-mean(i.c(2,:))), tracksAP2Master(candIdx)))
% minCh1 = max(arrayfun(@(i) min(i.A(1,:)+i.c(1,:)-mean(i.c(1,:))), tracksAP2Master(candIdx)))
% minCh2 = max(arrayfun(@(i) min(i.A(2,:)+i.c(2,:)-mean(i.c(2,:))), tracksAP2Master(candIdx)))
% candIdx = matchedIdx([102 142 139 137 138 309 127 421 605 624 850]);

% candIdx = [127 102 142 139 137 138 309  421 605 624 850];
candIdx = [127 138 309 605 421 850];

YLim = {[-80 450],[-10 60]};
YTick = {0:150:450, 0:15:60};
XTick = 0:20:120;

dt = 2; % frame rate
XLim = [-7*dt 100];

for i = 1:numel(candIdx)
    k = candIdx(i);
    ha = plotMatchedTracks(dataAP2Master(1), tracksAP2Master(matchedIdx(k)), tracksCLaMaster(idx{matchedIdx(k)}),...
        'YLim', YLim, 'YTick', YTick, 'XTick', XTick, 'XLim', XLim, 'Hues', [0.33 0 0.07]);
    if i~=1
        set(ha, 'YTickLabel', []);
    end
    %print('-depsc2', [f4path 'matchedTracks_' num2str(k, '%.4d') '.eps']);
end
% close all

