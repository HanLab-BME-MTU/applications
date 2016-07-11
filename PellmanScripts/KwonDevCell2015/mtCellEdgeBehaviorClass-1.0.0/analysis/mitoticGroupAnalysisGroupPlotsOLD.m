function [ output_args ] = mitoticGroupAnalysisGroupPlots(toPlot,varargin)
%% GCAGroupAnalysisGroupBoxplots: 
% Function Designed for Data Exploration: Plot one or more 
% 
% INPUT: toPlot  a structure with fields of parameters to plot such that
% toPlot.params{iGroup} = dataMat where each row is an observation and each
% column is data corresponding to a cell
% movie in the group. the key to this format is these
% can be easily horizontally cat to make a large array of all groups.
% this can be easily fed to boxplot to make individual boxplots and group
% boxplots
% toPlot.info gives information about the project Lists, the groupData
% names, the colors to use for the plots etc.

ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');

%ip.addParameter('Interactive',true,@(x) islogical(x)); 
ip.addParameter('Measurements',[]); % cell of measurements to be plotted  
ip.addParameter('outputDirectory',pwd,@(x) ischar(x));
% ip.addParameter('yLimOff',false,@(x) islogical(x)); 
% ip.addParameter('splitMovie',false);
% ip.addParameter('order',[]); % order to plot data from the cluster object 
ip.addParameter('plotType','pooled'); %'pooled' 
ip.addParameter('perMovieStat','nanmean' ); 
ip.parse(toPlot,varargin{:}); 



%%


% if ~isempty(ip.Results.Measurements);
%     params = ip.Results.Measurements;
% else
    params = fieldnames(toPlot);
    params = params(~strcmpi('info',params)) ;
    
%     if ip.Results.Interactive
%         paramSelect  = listSelectGUI(params,[],'move');
%         params  = params(paramSelect);
%     end
% end

perMovieStat = str2func(ip.Results.perMovieStat); 

%%
for iParam = 1:length(params)
    
    % get the dataMat (all groups horizontally catenated)
    
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam})));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}),'uniformoutput',0);
    
    dataMatLargeC = horzcat(forDataMat{:});
    
   
            %% 
       
            switch ip.Results.plotType
                case 'pooled'
                    %             %% start pooled.
                    setAxis('off')
                    colors  = vertcat(toPlot.info.colors{:});
                    grpNames = vertcat(toPlot.info.names(:));
                    
                    h1 =  boxplot(dataMatLargeC,toPlot.info.grouping,...
                        'colors',colors,'notch','on','outlierSize',1,'symbol','+','Labels',grpNames,'labelorientation','inline');
                    
                    set(h1(:),'Linewidth',2);
                    
                    forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.grouping==x),1:numel(grpNames),'uniformoutput',0);
                    forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
                    
                    cmpFunc = @nanmedian;
                    % for N median
                    percentChange = arrayfun(@(x) (cmpFunc(forN{x})-cmpFunc(forN{1}))./cmpFunc(forN{1}),2:numel(forN));
                    values = arrayfun(@(x) cmpFunc(forN{x}),1:numel(forN));
                    %toPlot.(params{iParam}).percentChange= percentChange;
                    
                    %     forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
                    N = cellfun(@(x) length(x(~isnan(x))),forN);
                    Nstring = num2cell(N);
                    Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
                    title(Nstring);
                    line([0.5,numel(grpNames)+0.5],[nanmedian(forN{1}),nanmedian(forN{1})],'color','k','linewidth',2);
                    
                    % Perform some simple stats against control
                    for i = 2:numel(grpNames)
                        [hit(i),pValue(i)] =   permTest(forN{1},forN{i},'CmpFunction',@nanmedian);
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
                    %
                    %     axis([0.5 length(toPlot.info.names) + 0.5 0 toPlot.(params{iParam}).ylim]);
                    % ylabel(toPlot.(params{iParam}).yLabel);
                    ylabel(strrep(params{iParam},'_',' '));
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.fig']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.png']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.eps'],'psc2');
                    
                    %% plot per cell
                case 'perCell'
                    % create the group folder (Veil/Branch/Filo) 
%                     if isfield(toPlot.(params{iParam}),'yGroup');
%                         cDir = [saveDir filesep toPlot.(params{iParam}).yGroup];
%                         if ~isdir(cDir)
%                             mkdir(cDir);
%                         end
%                         
%                     else
%                         cDir = saveDir;
%                     end
                    setAxis('on',0.75,20);
                    % get the cell array of per Movie values 
                    perCellValues = cellfun(@(x) perMovieStat(x,1)',toPlot.(params{iParam}),'uniformoutput',0); 
                    perCellDataMat = horzcat(perCellValues{:}); 
                    
                    
                    %nCells = length(unique(toPlot.info.groupingPerCell)); 
                    % collect the data forN (note 2016030 not sure why I
                    % formated this way-  it it seems odd 
                    % I compiled the groups only to again separate them...
                    % might just be due to history... check the split movie
                    % format) 
                    
                    % put into a 1xc cell where c is the total number of 
                    % neurite movies sampled over all groups 
                    % each cell contains N measurements compiled per movie
%                  
%                     forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPerCell==x),1:nCells,'uniformoutput',0);
%                     forN = cellfun(@(x) x(:),forN,'uniformoutput',0); 
                    
                    
                    
                    
                    % get the perNeurite stat of each neurite for the entire dataMat (condition info
                    % lost unfortunately) 1xc double where c the total number
                    % of neurite movies : the single row is distribution
                    % stat per neurite movie 
                    
%                     neuriteStatEachMovieCompile = cellfun(@(x) perMovieStat(x), forN);
                    
                    % add back the condition grouping info ( again I know a
                    % bit circuitous) 1xc cell where c is the number of
                    % perturbation groups, each cell is an rx1 double 
                    % where r is the number of movies sampled per
                    % perturbation- and each row is a distribution stat as 
                    % calculated for the entire neurite movie 
%                     neuriteStatEachMovieGrouped =  arrayfun(@(x) neuriteStatEachMovieCompile(toPlot.info.grouping==x)',1:numel(toPlot.info.names),'uniformoutput',0);
%                     
%                     % reformat cell to a padded double for notBoxplot
%                     % rxc matrix where r is the single measurement per movie (defined by the perMovieStat)
%                     % and c is the perturbation group 
%                     perCellDataMat = reformatDataCell(neuriteStatEachMovieGrouped);
                    
                    % sort the order of the perturbation groups if required
                    % (note this is helpful so that the order of groups reflects the 
                    % clustergram output order for easy comparison) 
%                     if ~isempty(ip.Results.order)
%                         namesC = ip.Results.order;
%                         IDSort = cellfun(@(x) find(strcmpi(x,toPlot.info.names)),namesC);
%                         IDSort = [1, IDSort];
%                         perCellDataMat = perCellDataMat(:,IDSort);
%                         names = ['KD Rac1' ; namesC']; 
%                         %names = ['Control' ; namesC'];
%                         colors = toPlot.info.color(IDSort);
%                         colorShades = toPlot.info.colorShades(IDSort);
%                         neuriteStatEachMovieGrouped = neuriteStatEachMovieGrouped(IDSort);
%                     else
                        colors = toPlot.info.colors;
                        names = toPlot.info.names;
%                         colorShades = toPlot.info.colorShades;
%                     end
                    
                    
                    h = notBoxPlot(perCellDataMat);
                  
                    %set the colors
                    arrayfun(@(i) set(h(i).data,'markerFaceColor', colors{i}),1:numel(names));
                    %arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(names));
                    arrayfun(@(i) set(h(i).mu,'color',colors{i}),1:numel(names));
                    %                     arrayfun(@(i) set(h(i).semPtch,'faceColor',colorShades{i}(4,:)),1:numel(names));
                    %                     arrayfun(@(i) set(h(i).sdPtch,'faceColor',colorShades{i}(1,:)),1:numel(names));
                    arrayfun(@(i) set(h(i).semPtch,'faceColor','none'),1:numel(names));
                    arrayfun(@(i) set(h(i).semPtch,'edgeColor','none'),1:numel(names));
                    arrayfun(@(i) set(h(i).sdPtch,'faceColor','none'),1:numel(names));
                    
                    % perform some quick stat tests
                    for i = 2:numel(toPlot.info.names)
                        [hit(i),pValues(i)] = permTest(perCellValues{1},perCellValues{i});
                        if hit(i) == 0
                            text(i,double(mean(perCellValues{i})),'NS','FontSize',14);
                        else
                            if pValues(i) == 0
                                text(i,double(mean(perCellValues{i})),num2str(pValues(i),1),'FontSize',14);
                            else
                                text(i,double(mean(perCellValues{i})),num2str(pValues(i),4),'FontSize',14);
                            end
                        end
                    end
                    forLabel = strrep(params{iParam},'_',' ');
                    
                    
                    ylabel(forLabel);
                    %ylabel(toPlot.(params{iParam}).yLabel);
                    set(gca,'XTick',1:numel(names));
                    set(gca,'XTickLabel',names,'FontSize',20);
                    
                    %axis([0.5 length(toPlot.info.names)+ 0.5 0 toPlot.(params{iParam}).ylim]);
                    
                    if isfield(toPlot.(params{iParam}),'yLabel');
                        ylabel(toPlot.(params{iParam}).yLabel);
                    end
                    set(gca,'XTickLabelRotation',45);
                    
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perCell' params{iParam} '.fig']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perCell' params{iParam} '.png']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perCell' params{iParam} '.eps'],'psc2');
                    close gcf 
            end  % switch ip.Results.plotType
        
  
end
end 

