% Francois Aguet, 02/01/2012

function plotTrackClasses(v, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('v');
ip.addOptional('v_std', []);
ip.addParamValue('Handle', []);
ip.addParamValue('YLim', [0 0.8]);
ip.parse(v, varargin{:});

xlabels = {'single tracks', 'single tracks, rej. gaps', 'single tracks, cut', 'single tracks, persistent',...
    'comp. tracks', 'comp. tracks, rej. gaps', 'comp. tracks, cut', 'comp. tracks, persistent'};


% Setup figure window
if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
else
    figure('Position', [440 378 400 300], 'PaperPositionMode', 'auto');
    ha = gca;
end

fset = loadFigureSettings();

barplot2(v, ip.Results.v_std, 'Handle', ha, 'BarWidth', 1.5, 'GroupDistance', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB,...
    'XLabels', xlabels, 'YTick', 0:0.1:0.8, 'YLim', ip.Results.YLim, 'YLabel', '% tracks');
