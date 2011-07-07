function plusTipMakeHistograms(rawData,saveDir,varargin)
% plusTipMakeHistograms saves speed, lifetime, and displacement histograms
% for growth, fgap, and bgap populations
%
% Input:
% 
%   rawData - A matrix containing the .
%
%   saveDir - a string containing the 
%
%   plotStd - 0 or 1 if the user wants to plot the standard deviation as an
%   errorbar
%
%   plotSde - 0 or 1 if the user wants to plot the standard error as an
%   errorbar
%
% Called by plusTipPoolGroupData.  See function for input format.
% Kathryn Applegate, Jan 2010
% Sebastien Besson, June 2011

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

% Define structures for each quantity/events
eventVals = [{'growth','fgap','bgap'};{[1 0 0],[0 0 1],[0 1 0]}];
eventType = cell2struct(eventVals,{'name','color'});

dataVals = [{'speed','lifetime','displacement'};...
    {'microns/min','sec','microns'};{1:3,4:6,7:9}];
dataType=cell2struct(dataVals,{'name','units','cols'});

plotData(3,3) = struct('rawData',[],'avgBins',[],'stdBins',[],'steBins',[]);

% Create anonymous function to parse raw data for each group/cell
parseSampleRawData = @(groupData,i,j) cellfun(@(x) x(:,dataType(i).cols(j)),...
    groupData,'UniformOutput',false);
parseRawData = @(i,j) cellfun(@(x) parseSampleRawData(x,i,j),rawData,...
    'UniformOutput',false);

% Read raw data for each couple quantity event.
for k=1:9
    [i,j]=ind2sub([3,3],k);
    plotData(i,j).rawData=parseRawData(i,j);
end

n=cell(1,3);
for i=1:3
    allCellsData=vertcat(plotData(i,:).rawData);
    allData = vertcat(allCellsData{:});
    n{i} = linspace(nanmin(vertcat(allData{:})), nanmax(vertcat(allData{:})),25);
end


for k=1:9
    [i,j]=ind2sub([3,3],k);
    % Create bin for each group/cell in a group
    binData = cellfun(@(x) cellfun(@(y) histc(y,n{i}),...
        x,'UniformOutput',false),plotData(i,j).rawData,'UniformOutput',false);
    plotData(i,j).binData = cellfun(@(x) horzcat(x{:}),binData,'UniformOutput',false);
    % Average the groups
    plotData(i,j).avgBins = cell2mat(cellfun(@(x) mean(x,2),plotData(i,j).binData,'Unif',false));
    plotData(i,j).stdBins = cell2mat(cellfun(@(x) std(x,0,2),plotData(i,j).binData,'Unif',false));
    plotData(i,j).steBins = cell2mat(cellfun(@(x) std(x,0,2)/sqrt(size(x,2)),plotData(i,j).binData,'Unif',false));
end

for k=1:9
    [i,j]=ind2sub([3,3],k);
    label=[regexprep(dataType(i).name,'(\<[a-z])','${upper($1)}') ' ('  ...
        dataType(i).units ')'];
    
    % Make individual plots for non-empty event data
    %     validBinData = find(~cellfun(@isempty,data(i).avgBins));
    %     for j=validBinData
    figure;hold on;
    if nGroups>1
        colors=varycolor(nGroups);
    else
        colors= eventType(j).color;
    end
    
    h=zeros(1,nGroups);
    for iGroup=1:nGroups
        dn=mean(diff(n{i})); % get bar spacing
        x = n{i}-dn/2+dn/(2*nGroups)+(iGroup-1)*dn/nGroups;
        hold on;
        h(iGroup)=bar(x,plotData(i,j).avgBins(:,iGroup),...
            'BarWidth',1/nGroups,'FaceColor',colors(iGroup,:));
        
        % Overlay standard error
        validSte = plotData(i,j).steBins(:,iGroup)~=0;
        if ~isempty(validSte) && plotSte
            errorbar(x(validSte),plotData(i,j).avgBins(validSte,iGroup),...
                plotData(i,j).steBins(validSte),'.k');
        end
        
        % Overlay standard error
        validStd = plotData(i,j).stdBins(:,iGroup)~=0;
        if ~isempty(validStd) && plotStd
            errorbar(x(validStd),plotData(i,j).avgBins(validStd,iGroup),...
                plotData(i,j).stdBins(validStd,iGroup),'.k');
        end
    end
    %         colormap(eventType(j).color);
    title([eventType(j).name ' ' dataType(i).name ' distribution'])
    xlabel(label);
    ylabel('Frequency of tracks');
    yLim = get(gca,'YLim');
    set(gca,'YLim',[0 yLim(2)]);
    if nGroups>1, legend(h,labels); end
    saveas(gcf,[saveDir filesep 'histogram_' dataType(i).name '_'...
        eventType(j).name '.tif'])
    close(gcf)

end

%     figure
%     bar(n,M,'stack')
%     colormap(vertcat(event.color))
%     legend({event.name},'Location','best')
%     title(['Stacked ' data(i).name ' distributions'])
%     xlabel(label);
%     ylabel('Frequency of tracks');
%     saveas(gcf,[saveDir filesep 'histogram_' data(i).name '_stacked.tif'])
%     close(gcf)