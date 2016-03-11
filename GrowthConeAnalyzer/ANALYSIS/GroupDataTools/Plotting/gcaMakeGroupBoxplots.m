function [ output_args ] = gcaMakeGroupBoxplots( toPlot,saveDir )
%INPUT: toPlot  a structure with fields of parameters to plot such that
% toPlot.params{iGroup} = dataMat where each row is an observation and each
% column is a cell/neurite in the group. the key to this format is these
% can be easily horizontally cat to make a large arrary of all groups that
% can be easily fed to boxplot to make individual boxplots and group
% boxplots
% toPlot.info gives information about the project Lists, the groupData
% names, the colors to use for the plots etc.

params = fieldnames(toPlot);
params = params(~strcmpi('info',params)) ;


for iParam = 1:length(params)
    % get the dataMat (all groups horizontally catenated)
    
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam})));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}),'uniformoutput',0);
    
    
    
    
    
    dataMatLargeC = horzcat(forDataMat{:});
    
    
    fsFigure(.75)
    setAxis
    projListAll = vertcat(toPlot.info.projList{:}) ;
    names = projListAll(:,2);
    names = strrep(names,'_',' ');
    colors = vertcat(toPlot.info.colors{:});
    %individual cells color coded by group.
    % if ~isempty(regexp(params{iParam},'Vel','Once'));
    %   dataMatLargeC =   dataMatLargeC.*216./5;
    % end
    
    
    h1 = boxplot(dataMatLargeC,'colorGroup',toPlot.info.grouping,'notch','on','outlierSize',1,'colors',colors,'Labels',names,'labelorientation','inline');
    set(h1(:),'Linewidth',2);
    
    saveas(gcf,[saveDir filesep 'perCellBoxplot' params{iParam} '.fig']);
    
    close gcf
    
    setAxis
    
    [~,numProjs] = size(dataMatLargeC);
    %values = arrayfun(@(x) prctile(dataMatLarge(:,x),95),1:numProjs);
    
    
    
    
    % get the per cell stats
    values = nanmean(dataMatLargeC,1);
    for iGroup = 1:numel(toPlot.info.names)
        y = values(toPlot.info.grouping == iGroup);
        x = repmat(iGroup,length(y),1);
        scatter(x,y,300,colors(iGroup,:),'filled','MarkerEdge','w');
        
        dataStructCell(iGroup).(params{iParam}) = y;
    end
    axis([0.5 numel(toPlot.info.names)+0.5 min(values(:))-1 max(values(:))+1]);
    
    % do some simple stats
    compMat = discriminationMatrix(dataStructCell);
    
    % get values above the diagonal
    pValues = compMat.(params{iParam})(1,2:end);
    pValuesTitle = mat2cell(pValues,1,ones(1,length(pValues)));
    pValuesTitle = cellfun(@(x) [ 'P-value PermTest = ' num2str(x,5)], pValuesTitle,'uniformoutput',0);
    %
    
    title(pValuesTitle,'FontSize',30);
    ylabel(params{iParam},'FontSize',30)
    set(gca,'XTick',[1:numel(toPlot.info.names)]);
    set(gca,'XTickLabels',[ toPlot.info.names  ]);
    saveas(gcf,[saveDir filesep 'perCellMeans' params{iParam}  '.fig']);
    close gcf
    
    %% pooled boxplots
    setAxis
    
    h = boxplot(dataMatLargeC,toPlot.info.grouping,'notch','on','outliersize',1,'colors',colors,'Labels',toPlot.info.names);
    set(h(:),'Linewidth',2);
    set(gca,'XTick',1:8);
    set(gca,'XTickLabel',toPlot.info.names,'FontSize',14);
    ylabel(params{iParam},'FontSize',30)
    % get the median of group 1
    values = toPlot.(params{iParam}){1}; % first group all projects
    forLine = nanmedian(values(:));
    %   if ~isempty(regexp(params{iParam},'Vel','Once'));
    %   forLine =   forLine.*216./5;
    % end
    
    line([0.5 numel(toPlot.info.names)+0.5] ,[forLine, forLine],'Linewidth',2,'color','k');
    
    forN = cellfun(@(x) x(:),toPlot.(params{iParam}),'uniformoutput',0);
    N = cellfun(@(x) length(x(~isnan(x))),forN);
    Nstring = num2cell(N);
    Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
    title(Nstring);
    
    saveas(gcf,[saveDir filesep 'pooled' params{iParam} '.fig']);
    
    
    
    
end


save('compMat.mat','compMat');






end

