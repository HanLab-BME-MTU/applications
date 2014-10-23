%=========================================================================================
% Figure 7
%=========================================================================================
% This script generates the epi:TIR panels for Figure 7 from Aguet et al., Dev. Cell, 2013.


%============================================================
% Panels D: 2nd epi/TIRF experiment
%============================================================
% 1) Load data
% dpath = '/home/fa48/files/LCCB/endocytosis/Marcel_fs003/121025_RPEegfpCLC_Epi_TIRF/';
dpath = '/groups/lccb/fa48/121025_RPEegfpCLC_Epi_TIRF/';

epiTIRFdata_ctrl = loadConditionData([dpath 'Control'], {'TIRF 488', 'Epi 488'}, {'EGFP', 'EGFP'});
epiTIRFdata_earless = loadConditionData([dpath 'deltaearrescue'], {'TIRF 488', 'Epi 488'}, {'EGFP', 'EGFP'});
epiTIRFdata_KD = loadConditionData([dpath 'control_norescue'], {'TIRF 488', 'Epi 488'}, {'EGFP', 'EGFP'});
% data = [epiTIRFdata_ctrl epiTIRFdata_earless epiTIRFdata_KD];
%getGaussianPSFsigmaFromData(epiTIRFdata_ctrl(4).framePaths{1}(1:10:51))
%getGaussianPSFsigmaFromData(epiTIRFdata_ctrl(4).framePaths{2}(1:10:51))
% runDetection(data, 'Sigma', [1.9 1.9]);
% runTracking(data, loadTrackSettings(), 'Overwrite', true);
% runTrackProcessing(data, 'Overwrite', true, 'FileName', 'ProcessedTracks_!isPSF.mat');

T = 74.5;
opts = {'MaxIntensityThreshold', T, 'Display', 'off', 'Rescale', true, 'RemoveOutliers', false, 'Overwrite', false};
lftResET_ctrl = runLifetimeAnalysis(epiTIRFdata_ctrl, opts{:});
lftResET_earlessH = runLifetimeAnalysis(epiTIRFdata_earless([1 3 5]), opts{:});
lftResET_earlessL = runLifetimeAnalysis(epiTIRFdata_earless([2 4 6]), opts{:});
lftResET_KD = runLifetimeAnalysis(epiTIRFdata_KD, opts{:});

%---------------------------------------------------------------
% 2) Calibration: intensity distr. of 1st frame (CCPs, >T only!)
%---------------------------------------------------------------
opts = {'ReturnValidOnly', true, 'Cutoff_f', 5, 'ExcludeVisitors', true};
lftDataEpiTIRctrl = getLifetimeData(epiTIRFdata_ctrl, opts{:});
lftDataEpiTIRdelta = getLifetimeData(epiTIRFdata_earless, opts{:});
% test scaling factors
% a = rescaleEDFs(maxActrlTIR, 'Display', true); % OK, std: 0.056
% a = rescaleEDFs(maxAdeltaTIR, 'Display', true) % std: 0.26
% a = rescaleEDFs(maxActrlEpi, 'Display', true); % OK, std: 0.047
% a = rescaleEDFs(maxAdeltaEpi, 'Display', true) % two groups
% two groups: high: 1,3,5 low: 2,4,6

maxA = arrayfun(@(i) max(i.A(:,:,1),[],2), lftDataEpiTIRctrl, 'unif', 0);
a = scaleEDFs(maxA, 'Display', false); %scaling looks good

maxA = arrayfun(@(i) max(i.A(:,:,1),[],2), lftDataEpiTIRdelta([1 3 5]), 'unif', 0);
aH = scaleEDFs(maxA, 'Display', false); %scaling looks good
maxA = arrayfun(@(i) max(i.A(:,:,1),[],2), lftDataEpiTIRdelta([2 4 6]), 'unif', 0);
aL = scaleEDFs(maxA, 'Display', false); %scaling looks good
aDelta([1 3 5]) = aH; aDelta([2 4 6]) = aL;

% apply to both channels
for i = 1:numel(lftDataEpiTIRctrl)
    lftDataEpiTIRctrl(i).A =   a(i) * lftDataEpiTIRctrl(i).A;
    lftDataEpiTIRctrl(i).sbA = a(i) * lftDataEpiTIRctrl(i).sbA;
    lftDataEpiTIRctrl(i).ebA = a(i) * lftDataEpiTIRctrl(i).ebA;
end
for i = 1:numel(lftDataEpiTIRdelta)
    lftDataEpiTIRdelta(i).A =   aDelta(i) * lftDataEpiTIRdelta(i).A;
    lftDataEpiTIRdelta(i).sbA = aDelta(i) * lftDataEpiTIRdelta(i).sbA;
    lftDataEpiTIRdelta(i).ebA = aDelta(i) * lftDataEpiTIRdelta(i).ebA;
end

Actrl = vertcat(lftDataEpiTIRctrl.A);
maxTIRCtrl = nanmax(Actrl(:,:,1),[],2);
[maxEpiCtrl, idx] = nanmax(Actrl(:,:,2),[],2); 
[nt,nf,nc] = size(Actrl);
TIRatMaxEpiCtrl = Actrl(sub2ind([nt,nf,nc], (1:nt)', idx, ones(nt,1)));

Adelta = vertcat(lftDataEpiTIRdelta([1 3 5]).A);
maxTIRDelta = nanmax(Adelta(:,:,1),[],2);
[maxEpiDelta, idx] = nanmax(Adelta(:,:,2),[],2); 
[nt,nf,nc] = size(Adelta);
TIRatMaxEpiDelta = Adelta(sub2ind([nt,nf,nc], (1:nt)', idx, ones(nt,1)));

% per-cell arrays for lifetime analysis
iMaxTIRCtrl = arrayfun(@(i) nanmax(i.A(:,:,1),[],2), lftDataEpiTIRctrl, 'unif', 0);
iMaxEpiCtrl = arrayfun(@(i) nanmax(i.A(:,:,2),[],2), lftDataEpiTIRctrl, 'unif', 0);
iMaxTIRDelta = arrayfun(@(i) nanmax(i.A(:,:,1),[],2), lftDataEpiTIRdelta, 'unif', 0);
iMaxEpiDelta = arrayfun(@(i) nanmax(i.A(:,:,2),[],2), lftDataEpiTIRdelta, 'unif', 0);

%------------------------------------------
% 2) Smoothing spline interpolation
%------------------------------------------
tlengths = vertcat(lftDataEpiTIRctrl.trackLengths);
AInterpCtrl = NaN(size(Actrl));
for ti = 1:numel(tlengths);
    idx = 1:tlengths(ti);
    AInterpCtrl(ti,idx,1) = interpB3SmoothingSpline1D(idx, Actrl(ti,idx,1), 1.5);
    AInterpCtrl(ti,idx,2) = interpB3SmoothingSpline1D(idx, Actrl(ti,idx,2), 1.5);
end

tlengths = vertcat(lftDataEpiTIRdelta([1 3 5]).trackLengths);
AInterpDelta = NaN(size(Adelta));
for ti = 1:numel(tlengths);
    idx = 1:tlengths(ti);
    AInterpDelta(ti,idx,1) = interpB3SmoothingSpline1D(idx, Adelta(ti,idx,1), 1.5);
    AInterpDelta(ti,idx,2) = interpB3SmoothingSpline1D(idx, Adelta(ti,idx,2), 1.5);
end

rLS = sum(AInterpCtrl(:,1:5,1).*AInterpCtrl(:,1:5,2),2)./sum(AInterpCtrl(:,1:5,1).^2,2);
rLSdAD = sum(AInterpDelta(:,1:5,1).*AInterpDelta(:,1:5,2),2)./sum(AInterpDelta(:,1:5,1).^2,2);
% noisier:
% rLS = sum(Actrl(:,1:5,1).*Actrl(:,1:5,2),2)./sum(Actrl(:,1:5,1).^2,2);
% rLSdAD = sum(Adelta(:,1:5,1).*Adelta(:,1:5,2),2)./sum(Adelta(:,1:5,1).^2,2);

maxTIRCtrlSmooth = nanmax(AInterpCtrl(:,:,1),[],2);
[maxEpiCtrlSmooth, idx] = nanmax(AInterpCtrl(:,:,2),[],2); 
maxTIRDeltaSmooth = nanmax(AInterpDelta(:,:,1),[],2);
[maxEpiDeltaSmooth, idx] = nanmax(AInterpDelta(:,:,2),[],2); 

%%
%==================================================
% Panel D: Epi/TIRF ratio vs. lifetime
%==================================================
fset = loadFigureSettings('print');
rv = 0:0.05:3;
lv = 0:2:120;
cutoff_s = 10;
T = 74.5;

mode = 'global';
switch mode
    case 'global'
        ratioVctrl = maxEpiCtrl ./ maxTIRCtrl * 5.5;
        ratioVdelta = maxEpiDelta ./ maxTIRDelta * 5.5;
        pct = prctile(maxTIRCtrl, 0);
    case 'smooth'
        ratioVctrl = maxEpiCtrlSmooth ./ maxTIRCtrlSmooth ./ rLS;
        ratioVdelta = maxEpiDeltaSmooth ./ maxTIRDeltaSmooth ./ rLSdAD;
        pct = prctile(maxTIRCtrlSmooth, 0);
end

% Control
lftVctrl = vertcat(lftDataEpiTIRctrl.lifetime_s);

% remove CCPs<T and CCPs with negative ratio
switch mode
    case 'global'
        idx = maxTIRCtrl>=T & ratioVctrl>0;
    case 'smooth'
        idx = maxTIRCtrlSmooth>=T & ratioVctrl>0;
end
lftVctrl = lftVctrl(idx);
ratioVctrl = ratioVctrl(idx);

% 99% confidence interval for ratio=1
tmp = ratioVctrl;
tmp = tmp(tmp<1 & lftVctrl>20)-1;
ciC = norminv(1-0.01/2,0,1)*std([tmp; -tmp]);

% CCPs > R=1+ci and in [1-ci,1+ci]
nd = numel(lftDataEpiTIRctrl);
cellIdxCtrl = arrayfun(@(i) i*ones(numel(lftDataEpiTIRctrl(i).lifetime_s),1), 1:nd, 'unif', 0);
cellIdxCtrl = vertcat(cellIdxCtrl{:});
na = arrayfun(@(i) sum(ratioVctrl(cellIdxCtrl(idx)==i)>1+ciC), 1:nd);
nb = arrayfun(@(i) sum(ratioVctrl(cellIdxCtrl(idx)==i)>1-ciC & ratioVctrl(cellIdxCtrl(idx)==i)<1+ciC), 1:nd);
fprintf('Curved CCPs WT:  %.1f±%.1f, flat CCPs: %.1f±%.1f (%.1f±%.1f%%; %d cells)\n',...
    mean(na), std(na), mean(nb), std(nb), mean(na./nb)*100, std(na./nb)*100, nd);

figure(fset.fOpts{:}, 'Position', [15 5 6.5 6.5]);
axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
hold on;
[~, rangeCtrl] = densityplot(lftVctrl, ratioVctrl, lv, rv,...
    'DisplayFunction', @sqrt, 'Div', numel(lftDataEpiTIRctrl));

plot([0 120], 1-ciC*[1 1], 'w--', 'LineWidth', 1);
plot([0 120], 1+ciC*[1 1], 'w--', 'LineWidth', 1);

set(gca, 'XTick', 0:20:120);
formatTickLabels(gca);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Epi:TIR ratio', fset.lfont{:});
%print('-depsc', '-loose', ['EpiTIRFRatioVsLft_ctrl_' mode '.eps']);


% Delta AD cells
smode = 'H';
switch smode
    case 'All'
        ai = 1:6;
    case 'H'
        ai = [1 3 5];
    case 'L'
        ai = [2 4 6];
end

lftVdelta = vertcat(lftDataEpiTIRdelta(ai).lifetime_s);

% remove CCPs<T and CCPs with negative ratio
switch mode
    case 'global'
        idx = maxTIRDelta>=T & ratioVdelta>0;
    case 'smooth'
        idx = maxTIRDeltaSmooth>=T & ratioVdelta>0;
end
lftVdelta = lftVdelta(idx);
ratioVdelta = ratioVdelta(idx);

% 99% confidence interval for ratio=1
tmp = ratioVdelta;
tmp = tmp(tmp<1 & lftVdelta>20)-1;
ciD = norminv(1-0.01/2,0,1)*std([tmp; -tmp]);

% CCPs > R=1+ci and in [1-ci,1+ci]
nd = numel(ai);
cellIdxDelta = arrayfun(@(i) i*ones(numel(lftDataEpiTIRdelta(ai(i)).lifetime_s),1), 1:nd, 'unif', 0);
cellIdxDelta = vertcat(cellIdxDelta{:});
na = arrayfun(@(i) sum(ratioVdelta(cellIdxDelta(idx)==i)>1+ciD), 1:nd);
nb = arrayfun(@(i) sum(ratioVdelta(cellIdxDelta(idx)==i)>1-ciD & ratioVdelta(cellIdxDelta(idx)==i)<1+ciD), 1:nd);
fprintf('Curved CCPs dAD: %.1f±%.1f, flat CCPs: %.1f±%.1f (%.1f±%.1f%%; %d cells)\n',...
    mean(na), std(na), mean(nb), std(nb), mean(na./nb)*100, std(na./nb)*100, nd);

figure(fset.fOpts{:}, 'Position', [15 5 6.5 6.5]);
axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
hold on;
[~, rangeDelta] = densityplot(lftVdelta, ratioVdelta, lv, rv,...
    'DisplayFunction', @sqrt, 'Div', numel(ai));

plot([0 120], 1-ciD*[1 1], 'w--', 'LineWidth', 1);
plot([0 120], 1+ciD*[1 1], 'w--', 'LineWidth', 1);
set(gca, 'XTick', 0:20:120);
formatTickLabels(gca);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Epi:TIR ratio', fset.lfont{:});
% title('\DeltaAD cells', fset.lfont{:});
%print('-depsc', '-loose', ['EpiTIRFRatioVsLft_delta_' mode '.eps']);


%%
%=============================================================
% Panel E: Lifetime comparison: > or < than Epi:TIR = 1
%=============================================================
cv = hsv2rgb([0 0 0;
              0.99 1 1;              
              0.55 1 1]);
T = 74.5;

mode = 'global';
switch mode
    case 'global'
        cf = 1/5.5;
        ratioVctrl = maxEpiCtrl ./ maxTIRCtrl ./ cf;
        ratioVdelta = maxEpiDelta ./ maxTIRDelta ./ cf;
    case 'smooth'
        ratioVctrl = maxEpiCtrlSmooth ./ maxTIRCtrlSmooth ./ rLS;
        ratioVdelta = maxEpiDeltaSmooth ./ maxTIRDeltaSmooth ./ rLSdAD;
end
iRatioCtrl = arrayfun(@(i) ratioVctrl(cellIdxCtrl==i), 1:max(cellIdxCtrl), 'unif', 0);
iRatioDelta = arrayfun(@(i) ratioVdelta(cellIdxDelta==i), 1:max(cellIdxDelta), 'unif', 0);

% rRange = [0.8:0.2:1.4 1.5];
rRange = 1.5;
% rRange = 1+max(ciC,ciD);
for ri = 1:numel(rRange)
    rf = rRange(ri);
    
    % create cell arrays for track indexes
    ratioCtrlGT = cellfun(@(x,y) x>=rf & y>=T, iRatioCtrl, iMaxTIRCtrl, 'unif', 0);
    ratioCtrlLT = cellfun(@(x,y) x<rf & y>=T, iRatioCtrl, iMaxTIRCtrl, 'unif', 0);
    ratioDeltaGT = cellfun(@(x,y) x>=rf & y>=T, iRatioDelta, iMaxTIRDelta(ai), 'unif', 0);
    ratioDeltaLT = cellfun(@(x,y) x<rf & y>=T, iRatioDelta, iMaxTIRDelta(ai), 'unif', 0);
    opts = {'MaxIntensityThreshold', T, 'Display', 'off', 'Rescale', false,...
        'ExcludeVisitors', true, 'RemoveOutliers', false, 'Overwrite', false};
    lftRes_ctrlGT = runLifetimeAnalysis(epiTIRFdata_ctrl, opts{:}, 'SelectIndex', ratioCtrlGT);
    lftRes_ctrlLT = runLifetimeAnalysis(epiTIRFdata_ctrl, opts{:}, 'SelectIndex', ratioCtrlLT);
    lftRes_deltaGT = runLifetimeAnalysis(epiTIRFdata_earless(ai), opts{:}, 'SelectIndex', ratioDeltaGT);
    lftRes_deltaLT = runLifetimeAnalysis(epiTIRFdata_earless(ai), opts{:}, 'SelectIndex', ratioDeltaLT);
    
    figure(fset.fOpts{:}, 'Position', [5 5 8 7]);
    axes(fset.axOpts{:});
    hold on;
    h(1) = plot(lftRes_ctrlGT.t, lftRes_ctrlGT.meanLftHistCCP, 'k', 'LineWidth', 1);
    h(2) = plot(lftRes_ctrlLT.t, lftRes_ctrlLT.meanLftHistCCP, '-', 'Color', 0.6*[1 1 1], 'LineWidth', 1);
    
    h(3) = plot(lftRes_deltaGT.t, lftRes_deltaGT.meanLftHistCCP, 'Color', cv(3,:), 'LineWidth', 1);
    h(4) = plot(lftRes_deltaLT.t, lftRes_deltaLT.meanLftHistCCP, 'Color', cv(2,:), 'LineWidth', 1);

    %median lifetimes
    [tmp, idx] = unique(2*cumsum(mean(lftRes_deltaGT.meanLftHistCCP,1)));
    interp1(tmp, lftRes_deltaGT.t(idx), 0.50)
    [tmp, idx] = unique(2*cumsum(mean(lftRes_deltaLT.meanLftHistCCP,1)));
    interp1(tmp, lftRes_deltaLT.t(idx), 0.50)

    axis([0 120 0 0.06]);
    set(gca, 'YTick', 0:0.01:0.06);
    
    na1 = cellfun(@(x) sum(x), ratioCtrlGT);
    na2 = cellfun(@(x) sum(x), ratioCtrlLT);
    na3 = cellfun(@(x) sum(x), ratioDeltaGT);
    na4 = cellfun(@(x) sum(x), ratioDeltaLT);
    rstd = @(x) 1/norminv(0.75, 0, 1) * mad(x, 1);
    rmean = @(x) round(mean(x)/100)*100;
    
    str = num2str(rf, '%.1f');
    fmt = '%.1f';
    fmt2 = '%.2f';
    hl = legend(h,...
        [' Control, Epi:TIR > ' str ' (~' num2str(mean(na1./(na1+na2))*100,fmt) '% CCPs)'],...
        [' Control, Epi:TIR < ' str ],...
        [' \DeltaAD \alpha-ad., Epi:TIR > ' str ' (~' num2str(mean(na3./(na3+na4))*100,fmt) '% CCPs)'],...
        [' \DeltaAD \alpha-ad., Epi:TIR < ' str]);
    set(hl, 'Box', 'off', 'Position', [1.5 5.25 2 1.5]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    %print('-depsc2', '-loose', ['EpiTIRFRatioLifetimes_avg_rT=' num2str(rf, '%.1f') '_' mode 'NEW.eps']);
end


%%
%=============================================================
% Panel F: Tfn uptake with LatA
%=============================================================
f7path = '/Users/aguet/Documents/MATLAB/endocytosis/CMEpaper/Figure 7/';
[num,txt,raw] = xlsread([f7path 'TfnRUptake_RPE_LatA_summary.xlsx'], 'I2:J23');

fset = loadFigureSettings('print');
xa = [0 2.5 5 10];
% labels = {};

mu = num(:,1); mu(isnan(mu)) = []; mu = reshape(mu, [4 4]);
sigma = num(:,2); sigma(isnan(sigma)) = []; sigma = reshape(sigma, [4 4]);

%normalize by FL
r = mu(4,1);
mu = mu/r*100;
sigma = sigma/r*100;

cv = hsv2rgb([0 0 0;
              0.99 1 1;              
              0.11 1 1;
              0.3 1 1;
              0.55 1 1]);

figure(fset.fOpts{:}, 'Position', [5 5 6.5 6.5]);
axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
hold on;

he(1) = errorbar(xa, mu(:,1), sigma(:,1), 'Color', [0 0 0], 'LineWidth', 1.5); % FL
he(2) = errorbar(xa, mu(:,2), sigma(:,2), '--', 'Color', 0.4*[1 1 1], 'LineWidth', 1.5); % FL LatA
he(3) = errorbar(xa, mu(:,3), sigma(:,3), 'Color', hsv2rgb([0.55 1 0.9]), 'LineWidth', 1.5); % dAD
he(4) = errorbar(xa, mu(:,4), sigma(:,4), '--', 'Color', hsv2rgb([0.6 1 1]), 'LineWidth', 1.5); % dAD LatA
arrayfun(@(i) setErrorbarStyle(i, 0.2), he);

set(gca, 'XTick', 0:2:10);
axis([0 10+10/50 0 130]);
xlabel('Uptake (min)', fset.lfont{:});
lg = {' FL \alpha-ad.', ' FL \alpha-ad. + LatA', ' \DeltaAD \alpha-ad.', ' \DeltaAD \alpha-ad. + LatA'};
hl = legend(lg, 'Location', 'SouthEast');
set(hl, 'Box', 'off', 'Position', [3.2 1.8 1.8 1.2]);
% print('-depsc2', '-loose', [f7path 'LatA.eps']);
