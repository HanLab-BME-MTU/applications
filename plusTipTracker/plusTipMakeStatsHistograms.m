function plusTipMakeStatsHistograms(rawData,saveDir,varargin)
% plusTipMakeStatsHistograms compares statistics between groups
%
% Input:
% 
%   rawData - A structure containing the stats field or a cell array or a
%   cell array of cell array (if comparing nGroups)
%
%   saveDir - a string containing the directory where to save the result
%
%   plotStd - 0 or 1 if the user wants to plot the standard deviation as an
%   errorbar
%
%   plotSde - 0 or 1 if the user wants to plot the standard error as an
%   errorbar
%
% Called by plusTipPoolGroupData. 
% Sebastien Besson, July 2011

% Check input
ip = inputParser;
ip.addRequired('rawData',@(x) iscell(x) || isnumeric(x));
ip.addRequired('saveDir',@ischar);
ip.addParamValue('plotStd',0,@isscalar);
ip.addParamValue('plotSte',1,@isscalar);
ip.addParamValue('labels',{},@iscell);
ip.parse(rawData,saveDir,varargin{:});
plotStd=ip.Results.plotStd;
plotSte=ip.Results.plotSte;
labels=ip.Results.labels;

% Convert input to generic cell array of size Nx1 where N is the number of
% groups. Each cell contains M(i)x1 cells where M(i) is the number of
% samples for the i-th group. Each sample cell should be a Lx9 matrix.
if isnumeric(rawData), 
    rawData={{rawData}}; 
elseif isnumeric(rawData{1})
    rawData={rawData};
end
nGroups = numel(rawData);

if ~isdir(saveDir), mkdir(saveDir); end

% Create anonymous function to parse raw data for each group/cell
parseSampleRawData = @(groupData) cellfun(@(x) x.stats.growth_speed_median,...
    groupData,'UniformOutput',false);
plotData.rawData = cellfun(@(x) cell2mat(parseSampleRawData(x)),rawData,...
    'UniformOutput',false);

% Read raw data for each couple quantity event.;
plotData.avgData= cellfun(@(x) mean(x),plotData.rawData);
plotData.stdData= cellfun(@(x) std(x),plotData.rawData);
plotData.steData= cellfun(@(x) std(x)/sqrt(size(x,1)),plotData.rawData);

figure;
hold on;

x=1:nGroups;
hold on;
bar(x,plotData.avgData);

% Overlay standard error
validSte = plotData.steData~=0;
if ~isempty(find(validSte,1)) && plotSte
    errorbar(x,plotData.avgData,plotData.steData,'.k');
end
        
% Overlay standard error
validStd = plotData.stdData~=0;
if ~isempty(find(validStd,1)) && plotStd
    errorbar(x,plotData.avgData,plotData.stdData,'.k');
end

title('Growth speed median comparison')
set(gca,'XTick',1:nGroups,'XTickLabel',labels)
ylabel('Growth speed median (microns/min)');
saveas(gcf,[saveDir filesep 'histogram_growth_speed_media.tif'])
close(gcf)
