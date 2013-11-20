%plotLifetimeComparison(lftRes, varargin) plots CCP lifetime distributions from different conditions
%
% Input:
%         lftRes : cell array of structures returned by runLifetimeAnalysis()
%         legend : string array of the same size as 'lftRes'
%
% Options ('specifier', value):
%      'Frequency' : 'relative' normalizes the distributions by proportion of CCPs
%        'PlotAll' : true|{false} also displays the lifetime distributions for all CCSs
%          'Color' : color matrix, Nx3
%
% Note: the control condition should be first in the 'lftRes' array

% Francois Aguet, 06/03/2013

function [ha] = plotLifetimeComparison(lftRes, varargin)

N = numel(lftRes);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftRes', @iscell);
ip.addOptional('legend', arrayfun(@(i) [' Condition ' num2str(i)], 1:N, 'unif', 0));
ip.addParamValue('PlotOrder', 1:N);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('Frequency', '', @(x) strcmpi(x, 'relative'));
ip.addParamValue('Color', []);
ip.addParamValue('SlaveName', [], @iscell);
ip.addParamValue('Control', []);
ip.parse(lftRes, varargin{:});

fset = loadFigureSettings('print');

% normalization
if strcmpi(ip.Results.Frequency, 'relative')
    w = cellfun(@(i) mean(i.pctCCP), lftRes);
    w = w/w(1);
else
    w = ones(N,1);
end
hp = zeros(N,1);

cv = ip.Results.Color;
if isempty(cv)
    cv = hsv2rgb([0 0 0;
        0.3 1 1;
        0.55 1 1;
        0.99 1 1;
        0.11 1 1]);
    if N>5
        cv = [0 0 0; jet(N-1)];
    end
end

% 1) Plot CCP distributions
if ~all(cellfun(@(i) isfield(i, 'lftHistSlaveCCP'), lftRes))
    ha = setupFigure;
    for i = ip.Results.PlotOrder
        hp(i) = plot(lftRes{i}.t, w(i)*lftRes{i}.meanLftHistCCP, '-', 'LineWidth', 1, 'Color', cv(i,:));
    end
    
else
    na = 2;
    % adapted from plotIntensityCohorts
    
%     na = 2;
%     ah = 1;
%     sigCombIdx = [1 0]';
%     
%     pct = zeros(nd,2);
%     for i = 1:nd
%         s =  lftData(i).significantMaster;
%         idx = lftData(i).maxA(:,1)>ip.Results.MaxIntensityThreshold;
%         pct(i,:) = sum([s(idx,2) ~s(idx,2)],1)/sum(idx);
%     end
%     meanPct = mean(pct,1);
%     stdPct = std(pct,[],1);
    
    pctS = cellfun(@(i) i.pctSlaveCCP/sum(i.pctSlaveCCP), lftRes, 'unif', 0);
    pctS = vertcat(pctS{:});
    meanPct = mean(pctS,1);
    stdPct = std(pctS,[],1);
    
    SlaveName = ip.Results.SlaveName;
    if isempty(SlaveName)
        SlaveName = cell(1,na);
    end
    sigCombIdx = [1 0]';
    tmp = sigCombIdx;
    tmp(tmp==1) = '+';
    tmp(tmp==0) = '-';
    atext = cell(1,na);
    switch na
        case 2
            for a = 1:na
                atext{a} = [tmp(a,1) SlaveName{1} ': ' num2str(meanPct(a)*100, '%.1f') '±' num2str(stdPct(a)*100, '%.1f') '%'];
            end
        case 3
            for a = 1:na
                atext{a} = [SlaveName{1} tmp(a,1) ' / ' SlaveName{2} tmp(a,2) ': '...
                    num2str(meanPct(a)*100, '%.1f') '±' num2str(stdPct(a)*100, '%.1f') '%'];
            end
    end
    
    
    ha = setupFigure(1,2, 'SameAxes', true, 'YSpace', [1.5 1 0.75]);
    
    legendText = cell(1,N);
    
    if ~isempty(ip.Results.Control)
        plot(ha(1), ip.Results.Control.t, ip.Results.Control.meanLftHistCCP, 'k', 'LineWidth', 1);
    end
    for i = ip.Results.PlotOrder
        hp(i) = plot(ha(1), lftRes{i}.t, w(i)*lftRes{i}.lftHistSlaveCCP{1}, '-', 'LineWidth', 1, 'Color', cv(i,:));
        legendText{i} = [' +' SlaveName{1} ', ' ip.Results.legend{i} ', ' num2str(100*pctS(i,1), '%.1f') '%'];
    end
    hl = legend(ha(1), hp, legendText);
    set(hl, 'Box', 'off', fset.sfont{:}, 'Units', 'centimeters',...
        'Position', [2.5 5-(N-1)*0.4 2 N*0.4]);
    
    if ~isempty(ip.Results.Control)
        plot(ha(2), ip.Results.Control.t, ip.Results.Control.meanLftHistCCP, 'k', 'LineWidth', 1);
    end
    for i = ip.Results.PlotOrder
        hp(i) = plot(ha(2), lftRes{i}.t, w(i)*lftRes{i}.lftHistSlaveCCP{2}, '-', 'LineWidth', 1, 'Color', cv(i,:));
        legendText{i} = [' -' SlaveName{1} ', ' ip.Results.legend{i} ', ' num2str(100*pctS(i,2), '%.1f') '%'];
    end
    hl = legend(ha(2), hp, legendText);
    get(hl)
    set(hl, 'Box', 'on', fset.sfont{:}, 'Units', 'centimeters',...
        'Position', [9.5 5-(N-1)*0.4 2 N*0.4], 'EdgeColor', 'w');
end

% pctS
% hl = legend(ha(end), hp, cellfun(@(i) [' ' i], legendText, 'unif', 0));
% set(hl, 'Box', 'off', fset.sfont{:}, 'Units', 'centimeters',...
%     'Position', [3.5 5-(N)*0.4 2 N*0.4]);
XLim = [0 160];
YLim = [0 0.03];
set(ha, 'XLim', XLim, 'YLim', YLim, 'XTick', 0:20:200, 'YTick', 0:0.01:0.1);
arrayfun(@(i) xlabel(i, 'Lifetime (s)', fset.lfont{:}), ha);
if strcmpi(ip.Results.Frequency, 'relative')
    ylabel(ha(1), 'Relative frequency', fset.lfont{:});
else
    ylabel(ha(1), 'Frequency', fset.lfont{:});
end

% for a = 1:na
%     text(XLim(1), 1.05*YLim(2), atext{a}, fset.sfont{:}, 'Parent', ha(a),...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
% end



% 2) Optional: plot CCS distributions (all objects)
if ip.Results.PlotAll
    setupFigure;
    hold on;
    for i = ip.Results.PlotOrder
        hp(i) = plot(lftRes{i}.t, nanmean(lftRes{i}.lftHist_Ia,1), '-', 'LineWidth', 1, 'Color', cv(i,:));
    end
    hl = legend(hp, ip.Results.legend);
    set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [4 5-N*0.3 2 N*0.3]);
    axis([0 120 0 0.05]);
    set(gca, 'YTick', 0:0.01:0.05);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
end
