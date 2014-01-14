%plotMaxIntVsLifetime(data, varargin) plots maximum CCS intensity vs. lifetime
%
% Inputs:
%          data : structure returned by loadConditionData()
%
% Options:
%            xl : vector of lifetime bin centers
%            xa : vector of intensity bin centers
%
% Options ('specifier', value):
%    
%      'ExcludeVisitors' : true|{false}
%             'Cutoff_f' : Minimum track length, in frames
%         'FirstNFrames' : Calculates maximum intensity over the first N frames only.
%                          If this is a vector, plots for the individual values are generated.
%      'DisplayFunction' : Intensity scaling for display. Default: @sqrt
%                          For linear intensities, use @(x) x
%              'Channel' : Channel number. Default: 1
%               'Legend' : Cell array of legend strings matching input data
%                'NormX' : {true}|false. Normalize density as a function of lifetime.
%                'Width' : Width of the figure grid. Default: 4

% Francois Aguet, 2011 (last modified 06/20/2013)

function ha = plotMaxIntVsLifetime(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addOptional('xl', []);
ip.addOptional('xa', []);
ip.addParamValue('ExcludeVisitors', false, @islogical);
ip.addParamValue('Cutoff_f', 1, @isscalar);
ip.addParamValue('FirstNFrames', [], @isvector);
ip.addParamValue('DisplayFunction', @sqrt);
ip.addParamValue('Channel', 1, @isposint);
ip.addParamValue('Legend', []);
ip.addParamValue('Parent', []);
ip.addParamValue('LifetimeData', 'LifetimeData.mat');
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat');
ip.addParamValue('PlotIndividual', false, @islogical);
ip.addParamValue('NormX', true, @islogical);
ip.addParamValue('FontSize', 10);
ip.addParamValue('Width', 4, @isposint);
ip.addParamValue('AmplitudeCorrection', []);
ip.parse(data, varargin{:});
ch = ip.Results.Channel;

if ip.Results.PlotIndividual
    lftData = getLifetimeData(data, 'Overwrite', false, 'Mask', true,...
        'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
        'ReturnValidOnly', false, 'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
    data = arrayfun(@(i) i, data, 'unif', 0);
    lftData = arrayfun(@(i) i, lftData, 'unif', 0);
    nd = numel(data);
else
    if ~iscell(data)
        data = {data};
    end
    nd = numel(data);
    lftData = cell(1,nd);
    for i = 1:nd
        lftData{i} = getLifetimeData(data{i}, 'Overwrite', false, 'Mask', true,...
            'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
            'ReturnValidOnly', false, 'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
    end
end

legendText = ip.Results.Legend;
if isempty(legendText)
    legendText =  cellfun(@(i) getDirFromPath(getExpDir(i)), data, 'unif', 0);
end

maxA = cell(nd,1);
maxALft = cell(nd,1);
lft = cell(nd,1);
for k = 1:nd
    na = numel(data{k});
    
    nCh = numel(data{k}(1).channels);
    
    lft{k} = cell(1,na);
    for i = 1:na
        lft{k}{i} = lftData{k}(i).lifetime_s;
        % lifetime at intensity maximum
        for c = 1:nCh
            if ~isempty(ip.Results.FirstNFrames)
                [tmp, maxIdx] = max(lftData{k}(i).A(:,ip.Results.Cutoff_f:ip.Results.FirstNFrames,c),[],2);
            else
                [tmp, maxIdx] = max(lftData{k}(i).A(:,ip.Results.Cutoff_f:end,c),[],2);
            end
            maxA{k}{i}(c,:) = tmp;
            maxALft{k}{i}(c,:) = (maxIdx-1)*data{k}(i).framerate;
        end
    end
end

xa = ip.Results.xa;
if isempty(xa)
    tmp = cellfun(@(i) i(ch,:), [maxA{:}], 'unif', 0);
    xa = linspace(0,prctile([tmp{:}],99.9),40);
end
xl = ip.Results.xl;
if isempty(xl)
    tmp = [lft{:}]; tmp = vertcat(tmp{:});
    xl = linspace(0,max(tmp),40);
end

ha = ip.Results.Parent;
if isempty(ha)
    ha = setupFigure(ceil(nd/ip.Results.Width), ip.Results.Width, nd, 'SameAxes', true,...
        'AxesWidth', 3, 'AxesHeight', 3, 'InsetPosition', [],...
        'XSpace', [1.5 0.5 0.5], 'YSpace', [1.5 1 0.5]);
end

colormap(jet(256));
for k = 1:nd
    av = [maxA{k}{:}];
    if size(av,1)>=ch
        av = av(ch,:);
        lv = lft{k};
        lv = vertcat(lv{:});
        %lv = [maxALft{k}{:}];
        %lv = lv(ch,:);
        
        densityplot(lv, av, xl, xa, 'Parent', ha(k), 'DisplayFunction', ip.Results.DisplayFunction, 'NormX', ip.Results.NormX);
        text(xl(end)/2, xa(end), legendText{k}, 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center', 'Parent', ha(k), 'FontSize', ip.Results.FontSize);
    end
end
formatTickLabels(ha);
