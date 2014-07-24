%=========================================================================================
% Figure 2
%=========================================================================================
% This script generates the panels for Figure 2 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f2path = '/Users/aguet/Documents/Papers/DevCell13/Figure 2 - Lifetime analysis/';

% root = '/Users/aguet/Documents/MATLAB/endocytosis/testData/';
root = '/home/fa48/lccb-endocytosis/';
% root = '/home/fa48/files/LCCB/endocytosis/';
dpath = [root 'Marcel_fs003/110301_BSC1_rat_monkey/BSC1_LCa_ratox/'];
dataOX = loadConditionData(dpath, {''}, {'EGFP'});


%%
%=====================================================
% Panels a & b: raw lifetime distributions
%=====================================================
opts = {'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat', 'Cutoff_f', 5};
lftResOX_0 = runLifetimeAnalysis(dataOX, opts{:}, 'ExcludeVisitors', false, 'Display', 'all');

% print('-depsc2', '-loose', 'lftMean+CDF_dataOX.eps');

%%
%=====================================================
% Panel c: EDF scaling
%=====================================================
lftDataOX = getLifetimeData(dataOX, 'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat',...
    'Cutoff_f', 0, 'ExcludeVisitors', false, 'Overwrite', false,...
    'Scale', false, 'DisplayScaling', true);
maxA_all = arrayfun(@(i) nanmax(i.A(:,:,1),[],2), lftDataOX, 'unif', 0);
[a, c] = scaleEDFs(maxA_all, 'Display', true, 'Reference', 'med',...
    'DisplayMode', 'print', 'XTick', 0:40:320);
if printFigures
    print('-depsc2', '-loose', [f2path 'ScaledEDF_dataOX.eps']);
end

%%
%=================================================================================
% Panel d: intensity vs. time heatmaps, per cohort. All cells.
%=================================================================================
% plotIntensityDistributions(dataOX);
% print('-depsc2', '-loose', 'trackIntensityMaps10s_all.eps');
plotIntensityDistributions(dataOX, 'ShowPct', true, 'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat');
if printFigures
    print('-depsc2', '-loose', [f2path 'trackIntensityMaps10s_all_med.eps']);
end
%%
%=====================================================
% Panel e: Max. intensity distributions
%=====================================================
plotMaxIntensityDistribution(dataOX, 'DisplayMode', 'print', 'XTick', 0:40:320,...
    'ShowSignificance', false, 'ShowGaussians', true, 'ShowFirstFrame', false,...
    'CohortLB', [1  11 16 21 41 61], 'CohortUB', [10 15 20 40 60 120], 'FirstNFrames', 6,...
    'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat');
if printFigures
    print('-depsc2', '-loose', [f2path 'MaxIntDistr_dataOX_6c_first6_gaussians.eps']);
end

%%
%==================================================================================
% Panel f: lifetime distribution after thresholding, before visitor exclusion
%==================================================================================
plotLifetimes(lftResOX_0, 'DisplayMode', 'print', 'PlotAll', true);
if printFigures
    print('-depsc2', '-loose', [f2path 'lftDist_MaxIntensityThresholdOnly.eps']);
end

% Exponential fit to CSs
% w = mean(lftResOX_0.pctCS);
% [mu,~,Aexp] = fitExpToHist(lftResOX_0.t, lftResOX_0.meanLftHistCS);
% plot(lftResOX_0.t, w*Aexp/mu*exp(-1/mu*lftResOX_0.t), 'r-', 'LineWidth', 1);
