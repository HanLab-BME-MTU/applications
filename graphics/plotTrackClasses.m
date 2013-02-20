% Francois Aguet, 02/01/2012

function [mu, sigma, hf] = plotTrackClasses(v, varargin)

fset = loadFigureSettings('print');

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('v');
ip.addOptional('c', []);
ip.addParamValue('Handle', []);
ip.addParamValue('YLim', [0 0.8]);
ip.addParamValue('YTick', 0:0.1:1);
ip.addParamValue('FaceColor', fset.cfTrackClasses);
ip.addParamValue('EdgeColor', fset.ceTrackClasses);
ip.parse(v, varargin{:});
c = ip.Results.c;

xlabels = {'single tracks', 'single tracks, rejected', 'single tracks, cut', 'single tracks, persistent',...
    'comp. tracks', 'comp. tracks, rejected', 'comp. tracks, cut', 'comp. tracks, persistent'};


% Setup figure window
if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
    hf = get(ha, 'Parent');
else
    hf = figure(fset.fOpts{:});
    ha = axes(fset.axOpts{:});
end

if iscell(v)
    vpos = cellfun(@(i,j) i(j==1), v, c, 'unif', 0);
    vneg = cellfun(@(i,j) i(j==0), v, c, 'unif', 0);
    v = arrayfun(@(i) hist(vpos{i}, 1:8)/numel(v{i}), 1:numel(v), 'unif', 0);
    v = vertcat(v{:});
end

mu = mean(v,1);
sigma = std(v,[],1);

barplot2(mu, sigma, 'Handle', ha, 'BarWidth', 0.6, 'GroupDistance', 0.8,...
    'FaceColor', ip.Results.FaceColor, 'EdgeColor', ip.Results.EdgeColor,...
    'XLabels', xlabels, 'YTick', 0:0.2:1, 'YLim', ip.Results.YLim);
ylabel('% tracks', fset.lfont{:});

% inset
% pos = get(gca, 'Position');
% pos = [pos(1)+pos(3)-160 pos(2)+pos(4)-110 150 100];
% ha = axes('Units', 'pixels', 'Position', pos);
% barplot2(v, v_std, 'Handle', ha, 'BarWidth', 1.5, 'GroupDistance', 1,...
%     'FaceColor', ip.Results.FaceColor, 'EdgeColor', ip.Results.EdgeColor,...
%     'AdjustFigure', false, 'XLabels', [], 'XLabel', 'Plaques', 'LabelFontSize', 14);


