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
    
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam}).dataMat));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}).dataMat,'uniformoutput',0);
    
    dataMatLargeC = horzcat(forDataMat{:});
   
    
%% Make the boxplot of the individual cells     
    %% Start boxplots : if user desires plot begin and end of movie separately for each neurite  
     setAxis('on'); 
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

   axis([0.5 size(dataMatLargeC,2)+ 0.5 0 toPlot.(params{iParam}).ylim]); 
   ylabel(toPlot.(params{iParam}).yLabel); 
    
   % h1 = boxplot(dataMatLargeC,'colorGroup',toPlot.info.grouping,'notch','on','outlierSize',1,'colors',colors,'Labels',names,'labelorientation','inline');
    set(h1(:),'Linewidth',2);
    
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.png']); 
    saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.eps'],'psc2'); 

    
    %% pooled boxplots
    setAxis('on')
    
    
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
    axis([0.5 nPooledGrps + 0.5 0 toPlot.(params{iParam}).ylim]); 
    ylabel(toPlot.(params{iParam}).yLabel); 
   
     
    
     
    
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
    
forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPoolBeginEndMovie==x),1:nPooledGrps,'uniformoutput',0); 
forN = cellfun(@(x) x(:),forN,'uniformoutput',0); 
    %forN = cellfun(@(x) x(:),toPlot.(params{iParam}).dataMat,'uniformoutput',0);

    N = cellfun(@(x) length(x(~isnan(x))),forN);
    Nstring = num2cell(N);
    Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
    title(Nstring);
    
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.fig']);
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.png']); 
    saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.eps'],'psc2'); 
    
    
    
end



%save('compMat.mat','compMat');






end

