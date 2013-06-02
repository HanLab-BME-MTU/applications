%[psnr] = getPSNRDistribution(data, varargin) plots the PSNR distribution for the detections in 'data'.

% Francois Aguet, 12/18/12

function ha = plotPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ha', [], @isnumeric);
ip.addParamValue('Color', hsv2rgb([0.33 1 0.9]));
ip.addParamValue('XLim', [1 300]);
ip.parse(data, varargin{:});
ha = ip.Results.ha;

psnr = arrayfun(@(i) getPSNRDistribution(i, 'Mode', 'all'), data, 'unif', 0);
psnr = [psnr{:}];

% PSNR values, in [dB]
dx = 0.2;
xi = 10.^((0:dx:30)/10); 
fi = hist(psnr, xi)/numel(psnr);

fset = loadFigureSettings('print');
if isempty(ha)
    figure(fset.fOpts{:});
    ha = axes(fset.axOpts{:}, 'xscale', 'log');
    hold on;
end
plot(xi, fi, '-', 'Color', ip.Results.Color, 'LineWidth', 1);
set(ha, 'XLim', ip.Results.XLim, 'XTick', 10.^(0:2), 'XTickLabel', {'1', '10', '100'},...
    fset.sfont{:});
xlabel('PSNR', fset.sfont{:});
ylabel('Frequency', fset.sfont{:});
formatTickLabels();
