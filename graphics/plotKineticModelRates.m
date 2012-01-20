%
% Inputs:
%         k : rate vector 
%     k_std : propagated standard deviation
%   corrMat : correlation matrix 
%
%

function XTickLabel = plotKineticModelRates(k, k_std, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('k');
ip.addRequired('k_std');
ip.addOptional('CorrMat', []);
ip.addParamValue('Handle', [], @ishandle);
ip.parse(k, k_std, varargin{:});

nk = numel(k);

fset = loadFigureSettings();


% sort rate vector depending on model (alpha: productive path; beta, other
switch nk
    case 1
        alpha = k;
        beta = [];
    case 3
        alpha = k(2:3);
        beta = k(1);
    case 5
        alpha = k([2 4 5]);
        beta = k([1 3]);
    case 7
        alpha = k([2 4 6 7]);
        beta = k([1 3 5]);
end

av = arrayfun(@(i) ['\alpha_' num2str(i)], 1:numel(alpha), 'UniformOutput', false);
bv = arrayfun(@(i) ['\beta_' num2str(i)], 1:numel(beta), 'UniformOutput', false);
XTickLabel = cell(1,nk);

ai = 1:2:nk;
bi = 2:2:nk;
XTickLabel(ai) = av;
XTickLabel(bi) = bv;


if isempty(ip.Results.Handle)
    figure('Position', [440 378 700 320], 'PaperPositionMode', 'auto');
    ha(1) = axes('Units', 'Pixels', 'Position', [70 60 360 240]);
else
    ha(1) = ip.Results.Handle;
end   
hold on;

xa = 1:nk;
plot(xa(ai), k(ai), '.', 'Color', 'k', 'MarkerSize', 20);
plot(xa(bi), k(bi), '.', 'Color', 0.4*[1 1 1], 'MarkerSize', 20);
errorbar(ha(1), xa(ai), k(ai), k_std(ai), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2);
errorbar(ha(1), xa(bi), k(bi), k_std(bi), 'Color', 0.4*[1 1 1], 'LineStyle', 'none', 'LineWidth', 2);

set(ha(1), 'XLim', [0.5 nk+0.5], 'XTickLabel', [], 'LineWidth', 2, fset.sfont{:}, 'XTick', 1:nk);
YLim = get(gca, 'YLim');
arrayfun(@(k) text(xa(k), YLim(1)-0.02*diff(YLim), XTickLabel{k},...
    fset.sfont{:},...
    'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
    'Interpreter', 'TeX'),...
    1:length(xa), 'UniformOutput', false);
ylabel('Rate (s^{-1})');

if ~isempty(ip.Results.CorrMat) && isempty(ip.Results.Handle)
    ha(2) = axes('Units', 'Pixels', 'Position', [480 60 180 240]);
    plotCorrelationMatrix(ip.Results.CorrMat, 'Handle', ha(2), 'TickLabels', XTickLabel);
end