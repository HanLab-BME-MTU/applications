function plotMaxIntVsLifetime(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addOptional('lb', [1  11 16 21 41 61]);
ip.addOptional('ub', [10 15 20 40 60 120]);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off', 'all'})));
ip.addParamValue('ExcludeVisitors', false, @islogical);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.parse(data, varargin{:});

opts = {'Scale', true, 'ReturnValidOnly', true, 'Cutoff_f', ip.Results.Cutoff_f,...
    'ExcludeVisitors', ip.Results.ExcludeVisitors};
lftData = getLifetimeData(data, opts{:});

lvec = 0:1:120;
avec = 0:2:300;

fset = loadFigureSettings('print');
maxACtrl = vertcat(lftData.maxA);
lftCtrl = vertcat(lftData.lifetime_s);

figure(fset.fOpts{:}, 'Position', [10 10 6.5 6.5]);
axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
hold on;
densityplot(lftCtrl, maxACtrl, lvec, avec, 'DisplayFunction', @log);
ylabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
xlabel('Lifetime (s)', fset.lfont{:});
set(gca, 'XTick', 0:20:120);