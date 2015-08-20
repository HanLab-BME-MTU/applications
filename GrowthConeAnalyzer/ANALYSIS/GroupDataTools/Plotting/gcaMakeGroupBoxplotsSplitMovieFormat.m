function [ output_args ] = gcaMakeGroupBoxplots( toPlot,varargin)
%INPUT: toPlot  a structure with fields of parameters to plot such that
% toPlot.params{iGroup} = dataMat where each row is an observation and each
% column is data corresponding to the 1st and 2nd half of the cell/neurite
% movie in the group. the key to this format is these
% can be easily horizontally cat to make a large array of all groups.
% this can be easily fed to boxplot to make individual boxplots and group
% boxplots
% toPlot.info gives information about the project Lists, the groupData
% names, the colors to use for the plots etc.

%NOTE: Currently uses the split movie treatment format
% (2 columns of data - first 5 and last 5- per movie)
% - might change but for now this format works well for the ARP23
%  treatment data and can be generalized to the rest with the correct grouping
%  framework .


%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot'); 
ip.addParameter('OutputDirectory',pwd,@(x) ischar(x));
ip.addParameter('yLimOff',false,@(x) islogical(x)); 

ip.parse(toPlot,varargin{:}); 
%%

params = fieldnames(toPlot);
params = params(~strcmpi('info',params)) ;
saveDir = ip.Results.OutputDirectory; 
for iParam = 1:length(params)
    
    % get the dataMat (all groups horizontally catenated)
    
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam}).dataMat));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}).dataMat,'uniformoutput',0);
    
    dataMatLargeC = horzcat(forDataMat{:});
    
    turnOff = 0;
    if turnOff ~= 1
    %% Make the boxplot of the individual cells
    % Start boxplots : if user desires plot begin and end of movie separately for each neurite
    setAxis('off');
    projListAll = vertcat(toPlot.info.projList{:}) ;
    names = projListAll(:,2);
    names = strrep(names,'_',' ');
    
    pre=  arrayfun(@(x) [names{x} ' PreTreat'],1:numel(names),'uniformoutput',0);
    post = arrayfun(@(x) [names{x} ' PostTreat'],1:numel(names),'uniformoutput',0);
    
    labels =   arrayfun(@(x,y) [x; y],pre,post,'uniformoutput',0);
    labels = vertcat(labels{:});
    labels = cell(labels) ;
    
    % create color shades
    for iGroup= 1: numel(toPlot.info.names)
        colorShades{iGroup}=  [toPlot.info.colorShades{iGroup}(end,:) ; toPlot.info.colorShades{iGroup}(7,:)];
    end
    colorShadesFinal = vertcat(colorShades{:});
    
    %individual cells color coded by group.
    % if ~isempty(regexp(params{iParam},'Vel','Once'));
    %   dataMatLargeC =   dataMatLargeC.*216./5;
    % end
    
    h1 =  boxplot(dataMatLargeC,'colorGroup',toPlot.info.groupingPoolBeginEndMovie,...
        'colors',colorShadesFinal,'notch','on','outlierSize',1,'symbol','+','Labels',labels,'labelorientation','inline');
    if (~isfield(toPlot.(params{iParam}),'ylim') || ip.Results.yLimOff == true) ; 
    else  
    axis([0.5 size(dataMatLargeC,2)+ 0.5 0 toPlot.(params{iParam}).ylim]);
    end 
    if isfield(toPlot.(params{iParam}),'yLabel'); 
    ylabel(toPlot.(params{iParam}).yLabel);
    else 
        ylabel(strrep(params{iParam},'_',' ')); 
    end 
    
    
    % perform the individual stats 
%     
%     for i = 2:2:size(dataMatLargeC)
%         permTest([dataMatLargeC(i,:), dataMatLargeC(i+1,:)],'Cmp 
%         
%     end 
    
    
    
    % 
    
%     
    % h1 = boxplot(dataMatLargeC,'colorGroup',toPlot.info.grouping,'notch','on','outlierSize',1,'colors',colors,'Labels',names,'labelorientation','inline');
    set(h1(:),'Linewidth',2);
    
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.png']);
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.eps'],'psc2');
    
    
    %% pooled boxplots
    setAxis('off')
    
    
    prePool=  cellfun(@(x) [x ' PreTreat'],toPlot.info.names,'uniformoutput',0);
    postPool = cellfun(@(x) [x ' PostTreat'],toPlot.info.names,'uniformoutput',0);
    labelsPool =   arrayfun(@(x,y) [x; y],prePool,postPool,'uniformoutput',0);
    labelsPool = vertcat(labelsPool{:});
    labelsPool = cell(labelsPool) ;
    
    
    h = boxplot(dataMatLargeC,toPlot.info.groupingPoolBeginEndMovie,'notch','on',...
        'symbol', '+','outliersize',1,'colors',colorShadesFinal, ...
        'Labels',labelsPool,'labelorientation','inline');
    
    
    set(h(:),'Linewidth',2);
    nPooledGrps = length(unique(toPlot.info.groupingPoolBeginEndMovie));
    
    
    if (~isfield(toPlot.(params{iParam}),'ylim') || ip.Results.yLimOff == true) ;
    else
        axis([0.5 nPooledGrps + 0.5 0 toPlot.(params{iParam}).ylim]);
    end
    
    
    if isfield(toPlot.(params{iParam}),'yLabel');
        ylabel(toPlot.(params{iParam}).yLabel);
    else
        ylabel(strrep(params{iParam},'_',' '));
    end
    
    
    
    
    
    set(gca,'XTick',1:nPooledGrps);
    set(gca,'XTickLabel',labelsPool,'FontSize',10);
    %     ylabel(params{iParam},'FontSize',30)
    % get the median of group 1
    %     values = toPlot.(params{iParam}){1}; % first group all projects
    %     forLine = nanmedian(values(:));
    %   if ~isempty(regexp(params{iParam},'Vel','Once'));
    %   forLine =   forLine.*216./5;
    % end
    
    %     line([0.5 numel(toPlot.info.names)+0.5] ,[forLine, forLine],'Linewidth',2,'color','k');
    
    % get the values for each group before and after treat
    % forN will be a cell 2*the number of conditions long (for before and after the condition)
    % each forN will contain a double matrix : max observation number x n neurite matrix  
    forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPoolBeginEndMovie==x),1:nPooledGrps,'uniformoutput',0);
    % combine the data for all neurites in the condition. 
    forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
    %forN = cellfun(@(x) x(:),toPlot.(params{iParam}).dataMat,'uniformoutput',0);
    
    medianValuesAll = cellfun(@(x) nanmedian(x),forN) ; 
    
    
    
    N = cellfun(@(x) length(x(~isnan(x))),forN);
    Nstring = num2cell(N);
    Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
    title(Nstring);
    
    pairs = unique(toPlot.info.groupingPoolBeginEndMovie);
    % row is before after treat, column is the condition 
    pairs = reshape(pairs,2,length(pairs)/2);
    
   
    
    
    for ipair = 1:size(pairs,2)
        [hit(ipair),pValues(ipair)] = permTest(forN{pairs(1,ipair)},forN{pairs(2,ipair)},'CmpFunction',@nanmedian,'nrep',5000);
        
        % plot line
        forLine = nanmedian(forN{pairs(1,ipair)});
        
        line([pairs(1,ipair)-0.25, pairs(2,ipair)+0.25],[forLine,forLine],...
            'color',colorShadesFinal(pairs(1,ipair),:),'linewidth',2);
        
        if hit(ipair)==0
            text(pairs(1,ipair),double(nanmedian(forN{pairs(1,ipair)}+2)),'NS','FontSize',10);
            % text(pairs(2,ipair),nanmedian(forN{pairs(2,ipair)}),'NS','FontSize',10);
        else
            text(pairs(1,ipair),double(nanmedian(forN{pairs(1,ipair)}+2)),num2str(pValues(ipair),4),'FontSize',10);
            
        end
        % test the difference between groups
        percentChange(ipair) = ( nanmedian(forN{pairs(2,ipair)}) -nanmedian(forN{pairs(1,ipair)})) / nanmedian(forN{pairs(1,ipair)});
     
    end
    
    toPlot.(params{iParam}).percentChange= percentChange;
    toPlot.(params{iParam}).pValues = pValues; % add a percent change field 
    clear percentChange
    
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.png']);
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.eps'],'psc2');
    end
    
    turnOff2 = 1; 
    if turnOff2 ~=1
    %% start pooled.
    setAxis('off')
    colors  = vertcat(toPlot.info.color{:});
    grpNames = vertcat(toPlot.info.names(:));
    h1 =  boxplot(dataMatLargeC,toPlot.info.groupingPoolWholeMovie,...
        'colors',colors,'notch','on','outlierSize',1,'symbol','+','Labels',grpNames,'labelorientation','inline');
    
    
%     axis([0.5 numel(grpNames)+ 0.5 0 toPlot.(params{iParam}).ylim]);
%     ylabel(toPlot.(params{iParam}).yLabel);
    set(h1(:),'Linewidth',2);
    
    forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPoolWholeMovie==x),1:numel(grpNames),'uniformoutput',0);
    forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
    %     forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
    N = cellfun(@(x) length(x(~isnan(x))),forN);
    Nstring = num2cell(N);
    Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
    title(Nstring);
    line([0.5,numel(grpNames)+0.5],[nanmedian(forN{1}),nanmedian(forN{1})],'color','k','linewidth',2);
    
    % Perform some simple stats against control
    for i = 2:numel(grpNames)
        [hit(i),pValue(i)] =   permTest(forN{1},forN{i},'CmpFunction',@median);
        if hit(i) == 0
            text(i,nanmedian(forN{i}),'NS','color','k','FontSize',10);
        else
            if pValue(i) == 0
                text( i,double(nanmedian(forN{i})),num2str(pValue(i),'%01d'),'color','k','FontSize',10);
            else
                text(i,double(nanmedian(forN{i})),num2str(pValue(i),'%04d'),'color','k','FontSize',10);
            end
        end % if hit
    end % for iGroup
    
%     axis([0.5 length(toPlot.info.names) + 0.5 0 toPlot.(params{iParam}).ylim]);
%     ylabel(toPlot.(params{iParam}).yLabel);
    
    saveas(gcf,[saveDir filesep 'pooled' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'pooled' params{iParam} '.png']);
    saveas(gcf,[saveDir filesep 'pooled' params{iParam} '.eps'],'psc2');
   
    %% plot per cell
    
    setAxis('off'); 
    nCells = length(unique(toPlot.info.groupingPerCell));
    forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPerCell==x),1:nCells,'uniformoutput',0);
    forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
    
    % get the means of each neurite for the entire dataMat (condigion info
    % lost unfortunately)
    means = cellfun(@(x) nanmean(x), forN);
    
    % add back the condition grouping info
    meansPerGp =  arrayfun(@(x) means(toPlot.info.grouping==x)',1:numel(toPlot.info.names),'uniformoutput',0);
    
    % reformat cell to double for notBoxplot
    % each column should now be the means/medians for each  cells in a condition
    % group
    perCellDataMat = reformatDataCell(meansPerGp);
    
    h = notBoxPlot(perCellDataMat);
    
    %set the colors
    arrayfun(@(i) set(h(i).data,'markerFaceColor', toPlot.info.color{i}),1:numel(toPlot.info.names));
    arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(toPlot.info.names));
    arrayfun(@(i) set(h(i).mu,'color',toPlot.info.color{i}),1:numel(toPlot.info.names));
    arrayfun(@(i) set(h(i).semPtch,'faceColor',toPlot.info.colorShades{i}(4,:)),1:numel(toPlot.info.names));
    arrayfun(@(i) set(h(i).sdPtch,'faceColor',toPlot.info.colorShades{i}(1,:)),1:numel(toPlot.info.names));
    
    
    % perform some quick stat tests
    for i = 2:numel(toPlot.info.names)
        [hit(i),pValues(i)] = permTest(meansPerGp{1},meansPerGp{i});
        if hit(i) == 0
            text(i,double(mean(meansPerGp{i})),'NS');
        else
            if pValues(i) == 0 
                text(i,double(mean(meansPerGp{i})),num2str(pValues(i),1)); 
            else 
            text(i,double(mean(meansPerGp{i})),num2str(pValues(i),4));
            end 
        end
    end
    
   % ylabel(toPlot.(params{iParam}).yLabel);
    set(gca,'XTick',1:numel(toPlot.info.names));
    set(gca,'XTickLabel',toPlot.info.names,'FontSize',10);
    
%         axis([0.5 length(toPlot.info.names)+ 0.5 0 toPlot.(params{iParam}).ylim]);
%     ylabel(toPlot.(params{iParam}).yLabel);
    
    saveas(gcf,[saveDir filesep 'perCell' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'perCell' params{iParam} '.png']);
    saveas(gcf,[saveDir filesep 'perCell' params{iParam} '.eps'],'psc2');
    
end

end

%save('compMat.mat','compMat');

% collect data 
 forHeatMap = arrayfun(@(x) toPlot.(params{x}).percentChange,1:numel(params),'uniformoutput',0); 
 forHeatMap = vertcat(forHeatMap{:}); 
close all  
 colorMap = HeatMap(forHeatMap,'RowLabels',params,'ColumnLabels',toPlot.info.names,... 
    'colormap','redbluecmap','ColumnLabelsRotate',45); 
save('colorMap','colorMap'); 
  saveas(gcf,'SummaryFigure.png');
  saveas(gcf,'SummaryFigure.fig'); 
  saveas(gcf,'SummaryFigure.eps','psc2'); 

save('toPlotWithStats','toPlot'); 



end

