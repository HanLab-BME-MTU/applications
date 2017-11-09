function [ pValues] = mitoticGroupAnalysisGroupPlots(toPlot,varargin)
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
ip.addParameter('plotType','perProject'); %'pooled' or 'perProject'
ip.addParameter('perMovieStat','nanmean' ); 
ip.addParameter('writeCSV',true,@(x) islogical(x)); 
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
    
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam}).dataMat));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}).dataMat,'uniformoutput',0);
    
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
                    for i = 1:numel(grpNames)
                        [hit(i),pValues(i)] =   permTest(forN{1},forN{i},'CmpFunction',@nanmedian);
                        if i>1
                            if hit(i) == 0
                                text(i,nanmedian(forN{i}),'NS','color','k','FontSize',10);
                            else
                                if pValue(i) == 0
                                    text( i,double(nanmedian(forN{i})),num2str(pValues(i),'%01d'),'color','k','FontSize',10);
                                else
                                    text(i,double(nanmedian(forN{i})),num2str(pValues(i),'%04d'),'color','k','FontSize',10);
                                end
                            end % if hit
                        end
                    end % for iGroup
                    %
                    %     axis([0.5 length(toPlot.info.names) + 0.5 0 toPlot.(params{iParam}).ylim]);
                    % ylabel(toPlot.(params{iParam}).yLabel);
                    ylabel(strrep(params{iParam},'_',' '));
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.fig']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.png']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'pooled' params{iParam} '.eps'],'psc2');
                    
                    if ip.Results.writeCSV
                        % dataMatLargeC 
                        % test operating system 
                        % if windows write xlssheets 
                        % if linux write each data set as an individual csv file. 
                        
                        csvDir = [ip.Results.outputDirectory filesep 'CSVFiles' filesep 'pooled' filesep params{iParam}];
                        if ~isdir(csvDir)
                            mkdir(csvDir)
                        end
                            
                        for iGroup = 1:numel(toPlot.info.names)                      
                            % transfer the dataMat per project to a dataset
                            % array for easy export to a .csv file. 
                           data = toPlot.(params{iParam}).dataMat{iGroup}'; 
                              toExport = mat2dataset(data ,'ObsNames',toPlot.info.projList{iGroup}(:,1));
                              
                              export(toExport,'file',[csvDir filesep toPlot.info.names{iGroup} '_' params{iParam}  '.csv']);  
                        end 
                    
                    end % ip.Results.writeCSV
                    
                    %% plot per cell
                case 'perProject'
                  
                    setAxis('on',0.75,20);
                    % get the cell array of per Movie values 
                    perCellValues = cellfun(@(x) perMovieStat(x,1)',toPlot.(params{iParam}).dataMat,'uniformoutput',0); 
                    perCellDataMat = horzcat(perCellValues{:}); 
                    
                    colors = toPlot.info.colors;
                    names = toPlot.info.names;
                   
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
                    for i = 1:numel(toPlot.info.names)
                        [hit(i),pValues(i)] = permTest(perCellValues{1},perCellValues{i});
                        if i > 1
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
                    end
                    
                   
                    
                    forLabel = strrep(params{iParam},'_',' ');
                    
                    
                    ylabel(forLabel);
                  
                    set(gca,'XTick',1:numel(names));
                    set(gca,'XTickLabel',names,'FontSize',20);
                    
                    %axis([0.5 length(toPlot.info.names)+ 0.5 0 toPlot.(params{iParam}).ylim]);
                    
                    if isfield(toPlot.(params{iParam}),'yLabel');
                        ylabel(toPlot.(params{iParam}).yLabel);
                    end
                    set(gca,'XTickLabelRotation',45);
                    
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perProject' params{iParam} '.fig']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perProject' params{iParam} '.png']);
                    saveas(gcf,[ip.Results.outputDirectory filesep 'perProject' params{iParam} '.eps'],'psc2');
                    close gcf 
                    
                    if ip.Results.writeCSV
                        % dataMatLargeC 
                        % test operating system 
                        % if windows write xlssheets 
                        % if linux write each data set as an individual csv file. 
                        
                        csvDir = [ip.Results.outputDirectory filesep 'CSVFiles' filesep 'perProject'];
                        if ~isdir(csvDir)
                            mkdir(csvDir)
                        end
                            
                        for iGroup = 1:numel(toPlot.info.names)                      
                            % transfer the dataMat per project to a dataset
                            % array for easy export to a .csv file. 
                           
                              toExport = mat2dataset(perCellValues{iGroup},'ObsNames',toPlot.info.projList{iGroup}(:,1));
                              
                              export(toExport,'file',[csvDir filesep toPlot.info.names{iGroup} '_' params{iParam}  '.csv']);  
                        end 
                    
                    end % ip.Results.writeCSV
                    
                    
                    
            end  % switch ip.Results.plotType
        
  
end
save('pValues.mat','pValues') 
end 

