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

function plotLifetimeComparison(lftRes, varargin)

N = numel(lftRes);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftRes', @iscell);
ip.addOptional('legend', arrayfun(@(i) [' Condition ' num2str(i)], 1:N, 'unif', 0));
ip.addParamValue('PlotOrder', 1:N);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('Frequency', '', @(x) strcmpi(x, 'relative'));
ip.addParamValue('Color', []);
ip.parse(lftRes, varargin{:});

fset = loadFigureSettings('print');

% 1) Plot CCP distributions
setupFigure;
hold on;

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

for i = ip.Results.PlotOrder
    hp(i) = plot(lftRes{i}.t, w(i)*lftRes{i}.meanLftHistCCP, '-', 'LineWidth', 1, 'Color', cv(i,:));
end

hl = legend(hp, ip.Results.legend);
set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [4 5-N*0.3 2 N*0.3]);
axis([0 120 0 0.05]);
set(gca, 'YTick', 0:0.01:0.05);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});


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
