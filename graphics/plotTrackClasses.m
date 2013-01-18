% Francois Aguet, 02/01/2012

function hf = plotTrackClasses(v, varargin)

fset = loadFigureSettings();

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('v');
ip.addOptional('v_std', []);
ip.addParamValue('Handle', []);
ip.addParamValue('YLim', [0 0.8]);
ip.addParamValue('YTick', 0:0.1:1);
ip.addParamValue('FaceColor', fset.cfTrackClasses);
ip.addParamValue('EdgeColor', fset.ceTrackClasses);
ip.parse(v, varargin{:});

xlabels = {'single tracks', 'single tracks, rejected', 'single tracks, cut', 'single tracks, persistent',...
    'comp. tracks', 'comp. tracks, rejected', 'comp. tracks, cut', 'comp. tracks, persistent'};


% Setup figure window
if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
    hf = get(ha, 'Parent');
else
    hf = figure('Position', [440 378 400 300], 'PaperPositionMode', 'auto');
    ha = gca;
end


v_std = ip.Results.v_std;
if all(v_std==0)
    v_std = [];
end

barplot2(v, v_std, 'Handle', ha, 'BarWidth', 0.8,...
    'FaceColor', ip.Results.FaceColor, 'EdgeColor', ip.Results.EdgeColor,...
    'XLabels', xlabels, 'YLabel', '% tracks');
set(ha, 'LineWidth', 2);


% inset
% pos = get(gca, 'Position');
% pos = [pos(1)+pos(3)-160 pos(2)+pos(4)-110 150 100];
% ha = axes('Units', 'pixels', 'Position', pos);
% barplot2(v, v_std, 'Handle', ha, 'BarWidth', 1.5, 'GroupDistance', 1,...
%     'FaceColor', ip.Results.FaceColor, 'EdgeColor', ip.Results.EdgeColor,...
%     'AdjustFigure', false, 'XLabels', [], 'XLabel', 'Plaques', 'LabelFontSize', 14);


