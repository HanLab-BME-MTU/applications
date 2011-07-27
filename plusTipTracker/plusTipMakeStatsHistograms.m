function plusTipMakeStatsHistograms(rawData,saveDir,varargin)
% plusTipMakeStatsHistograms compares statistics between groups
%
% Synopsis:  
%   plusTipMakeStatsHistograms(rawData,saveDir)
%   plusTipMakeStatsHistograms(rawData,saveDir,'plotStd',1,'plotSte',1)
%   plusTipMakeStatsHistograms(rawData,saveDir,'fields','growth_speed_median')
%   plusTipMakeStatsHistograms(rawData,saveDir,'fields','growth_speed_median','names','Growth speed median','units)
%   plusTipMakeStatsHistograms(rawData,saveDir,'fields',fieldnames(rawData{1}{1}.stats))
%
% Input:
% 
%   rawData - A structure containing the stats field or a cell array or a
%   cell array of cell array (if comparing nGroups)
%
%   saveDir - a string containing the directory where to save the result
%
%   Possible additional arguments (as parameter/value pairs)
%
%       plotStd - 0 or 1 if the user wants to plot the standard deviation
%       as an errorbar
%
%       plotSde - 0 or 1 if the user wants to plot the standard error as an
%       errorbar
%  
%       labels - a string or a cell array of string for the label of the
%       groups
%
%       fields - a string of a cell array of string containing the
%       statistics fields to plot as histograms
%
%       names - a string of a cell array of string containing the names of
%       the fields to plot (for axis labels and title)
%
%       units - a string of a cell array of string containing the units of
%       the fields to plot (for axis labels)
%
% Sebastien Besson, July 2011

% Check input
ip = inputParser;
ip.addRequired('rawData',@(x) iscell(x) || isnumeric(x));
ip.addRequired('saveDir',@ischar);
ip.addParamValue('plotStd',0,@isscalar);
ip.addParamValue('plotSte',1,@isscalar);
ip.addParamValue('labels',{},@iscell);
ip.addParamValue('fields','growth_speed_median',@(x) iscell(x) || ischar(x));
ip.addParamValue('names',{},@(x) iscell(x) || ischar(x));
ip.addParamValue('units',{},@(x) iscell(x) || ischar(x));
ip.parse(rawData,saveDir,varargin{:});
plotStd=ip.Results.plotStd;
plotSte=ip.Results.plotSte;
labels=ip.Results.labels;
fields=ip.Results.fields;
names=ip.Results.names;
units=ip.Results.units;

% Convert input to generic cell array of size Nx1 where N is the number of
% groups. Each cell contains M(i)x1 cells where M(i) is the number of
% samples for the i-th group. Each sample cell should be a Lx9 matrix.
if isnumeric(rawData), 
    rawData={{rawData}}; 
elseif isnumeric(rawData{1})
    rawData={rawData};
end

% Create save directory if absen
if ~isdir(saveDir), mkdir(saveDir); end

% Convert fields names  & units into cell arrays
if ischar(fields), fields={fields};  end
if ischar(names), names={names};  end
if ischar(units), units={units};  end

% Set default values to names and units
if isempty(names),names=fields; end
if isempty(units),units=cell(1,numel(fields)); end

% Call the main plotting function
arrayfun(@(i)plotHistogram(rawData,fields{i},saveDir,plotStd,plotSte,...
    labels,names{i},units{i}),1:numel(fields))

end

function plotHistogram(rawData,field,saveDir,plotStd,plotSte,labels,name,unit)

% Create anonymous function to parse raw data for each group/cell
parseSampleRawData = @(groupData) cellfun(@(x) x.stats.(field),...
    groupData,'UniformOutput',false);
rawData = cellfun(@(x) cell2mat(parseSampleRawData(x)),rawData,...
    'UniformOutput',false);

% Read raw data for each couple quantity event.;
plotData= cellfun(@(x) mean(x),rawData);
stdData= cellfun(@(x) std(x),rawData);
if nnz(stdData)~=0 && plotStd, plotData = [plotData;stdData]; end
steData= cellfun(@(x) std(x)/sqrt(size(x,1)),rawData);
if nnz(steData)~=0 && plotSte, plotData = [plotData;steData]; end


% Create initial bar plot
% figure;
% nGroups = numel(rawData);
% x=1:nGroups;
% bar(x,plotData.avgData);
% hold on;
% 
% % Overlay standard error
% validSte = plotData.steData~=0;
% if ~isempty(find(validSte,1)) && plotSte
%     errorbar(x(validSte),plotData.avgData(validSte),plotData.steData(validSte),'.k');
% end
%         
% % Overlay standard deviation
% validStd = plotData.stdData~=0;
% if ~isempty(find(validStd,1)) && plotStd
%     errorbar(x(validStd),plotData.avgData(validStd),plotData.stdData(validStd),'.k');
% end
% 
% % Additional graphic options (does not use the interpreter b
% title([name ' comparison'],'Interpreter','none');
% set(gca,'XTick',1:nGroups,'XTickLabel',labels);
% ylabel([name ' (' unit ')'],'Interpreter','none');
% saveas(gcf,[saveDir filesep 'histogram_' field '.tif'])
figure;barplot2(plotData,'xLabels',labels,...
    'ylabel',[name ' (' unit ')']);
print('-dtiff', '-r300',[saveDir filesep 'histogram_' field '.tif']);

close(gcf)
end