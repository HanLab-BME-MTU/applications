%plotMaxIntVsLifetime(data, varargin) plots maximum CCS intensity vs. lifetime
%
% Inputs:
%          data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%    
%      'ExcludeVisitors' : true|{false}
%             'Cutoff_f' : Minimum track length, in frames
%         'FirstNFrames' : Calculates maximum intensity over the first N frames only.
%                          If this is a vector, plots for the individual values are generated.

% Francois Aguet, 2011 (last modified 06/20/2013)

function plotMaxIntVsLifetime(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addParamValue('ExcludeVisitors', false, @islogical);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.addParamValue('FirstNFrames', [], @isvector);
ip.parse(data, varargin{:});

opts = {'Scale', true, 'ReturnValidOnly', true, 'Cutoff_f', ip.Results.Cutoff_f,...
    'ExcludeVisitors', ip.Results.ExcludeVisitors};
lftData = getLifetimeData(data, opts{:});

lvec = 0:1:120;
avec = 0:2:300;

fset = loadFigureSettings('print');
A = vertcat(lftData.A);
lft = vertcat(lftData.lifetime_s);

if isempty(ip.Results.FirstNFrames)
    maxA = max(A,[],2);
    
    figure(fset.fOpts{:}, 'Position', [10 10 6.5 6.5]);
    axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
    hold on;
    densityplot(lft, maxA, lvec, avec, 'DisplayFunction', @log);
    ylabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    set(gca, 'XTick', 0:20:120);
else
    for i = 1:numel(ip.Results.FirstNFrames)
        maxA = max(A(:,1:ip.Results.FirstNFrames(i)),[],2);
        figure(fset.fOpts{:}, 'Position', [10 10 6.5 6.5]);
        axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6); %#ok<LAXES>
        hold on;
        densityplot(lft, maxA, lvec, avec, 'DisplayFunction', @log);
        ylabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
        xlabel('Lifetime (s)', fset.lfont{:});
        set(gca, 'XTick', 0:20:120);
    end
end
