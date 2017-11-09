function [ zValuesInd ] = gcaHeatMap(toPlot,varargin)
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

%% Note as of 20160329 please revamp split movie format 
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');
ip.addParameter('Interactive',true,@(x) islogical(x));
ip.addParameter('OutputDirectory',pwd,@(x) ischar(x));
ip.addParameter('clearUnusedFields',true);  % Will rewrite the toPlot file just with 

% the user defined measurements 

% ip.addParameter('yLimOff',false,@(x) islogical(x));
ip.addParameter('diffMetric','perCellZ',@(x) ischar(x));           
            % 'perCellZ : the perCell metric is performed for each neurite movie 
                        %  to obtain a distribution N = number of cells per group  
                        %  the difference between the control distribution and the perturbation 
                        %  distribution is calculated. 
                        %  This difference is normalized by the standard deviation of
                        %  the control distribution.
            % 'pooledZ'   : 
            
%ip.addParameter('splitMovie',false);
ip.addParameter('perNeuriteStat','nanmean');
ip.addParameter('calcZ',true); % adds a z-score to the structure

% Plotting Parameters 
ip.addParameter('makePlots',true); % make group plots of those that are significant 
ip.addParameter('plotCutOff',1); % value for making plots  (Currently assumes one std) 

ip.addParameter('cutOffPValue',[]); % filter z-score by permTest p-value before clustering. 
% default don't filter 

ip.addParameter('filterNonSigFeatures',false); 

ip.addParameter('plotImagesc',true);% 

ip.parse(toPlot,varargin{:});
%%
toPlot = helperAddYLabels(toPlot); 
perNeuriteStat = str2func(ip.Results.perNeuriteStat); 

params = fieldnames(toPlot);
params = params(~strcmpi('info',params)) ;
paramSelectLogic = false(size(params,1),1); 

if ip.Results.Interactive
    paramSelect  = listSelectGUI(params,[],'move');
    
    params  = params(paramSelect);
    paramSelectLogic(paramSelect) = true; 
else 
    paramSelectLogic(1:length(paramSelectLogic)) = true; 
end

switch ip.Results.diffMetric
    case 'perCellZ'
        outDir = [ip.Results.OutputDirectory filesep 'HierarchicalCluster' ...
            filesep ip.Results.diffMetric  filesep ip.Results.perNeuriteStat ];
    case 'pooledZ'
        outDir = [ip.Results.OutputDirectory filesep 'HierarchicalCluster' ...
            filesep ip.Results.diffMetric  ];
end
    

if ~isdir(outDir) 
    mkdir(outDir)
end 

if ip.Results.calcZ
% %% SPLIT MOVIE     
%     if ip.Results.splitMovie
%         for iParam = 1:length(params)
%   %% KS TEST : SPLIT            
%             if strcmpi(ip.Results.diffMetric,'ks');
%                 for iGroup = 1:numel(toPlot.info.names);
%                     
%                     % for now 20151015 lets just make the before and after z score
%                     dataMat = toPlot.(params{iParam}).dataMat{iGroup};
%                     
%                     %ks = arrayfun(@(x) distribTest(dataMat(:,x),dataMat(:,x-1)),2:2:size(dataMat,2));
%                     [~,~,ks] = arrayfun(@(x) kstest2(dataMat(:,x),dataMat(:,x-1)),2:2:size(dataMat,2));
%                     % get a directionality : for now just just use the
%                     % delta in the mean ? 
%                     
%                     signs =  arrayfun(@(x) sign(nanmedian(dataMat(:,x))-nanmedian(dataMat(:,x-1))),2:2:size(dataMat,2));   %
%                     
%                     ks = signs.*ks; 
%                     
%                     toPlot.(params{iParam}).ks{iGroup} = ks;
%                     clear ks
%                 end
%             end % strcmpi
%             
%  %% Permutation Test of the Medians   SPLIT           
%             if strcmpi(ip.Results.diffMetric,'medPermTest');
%                 for iGroup = 1:numel(toPlot.info.names);
%                     
%                     % for now 20151015 lets just make the before and after z score
%                     dataMat = toPlot.(params{iParam}).dataMat{iGroup};
%                     
%                     [~,diffMed] = arrayfun(@(x) permTest(dataMat(:,x),dataMat(:,x-1),'CmpFunction',@median),2:2:size(dataMat,2));
%                    
%                     toPlot.(params{iParam}).permTestMed{iGroup} = diffMed;
%                     clear ks
%                 end
%                 
%                 forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).permTestMed{:}),1:numel(params),'uniformoutput',0); 
%                 
%                 
%             end % strcmpi
%  %% Per Cell "Z-Scores"           : SPLIT  NOTE CHECK            
%             
%             if strcmpi(ip.Results.diffMetric,'perCellZ');
%                 for iGroup = 1:numel(toPlot.info.names);
%                     
%                     % for now 20151015 lets just make the before and after z score
%                     dataMat = toPlot.(params{iParam}).dataMat{iGroup};
%                     
%             % CHECK THIS 20160329        
%                     zValuesPerNeurite = arrayfun(@(x) (nanmean(dataMat(:,x)) - nanmean(dataMat(:,x-1)))./nanstd(dataMat(:,x-1)),2:2:size(dataMat,2));
%                     zValuesInd.(params{iParam}){iGroup} = zValuesPerNeurite; 
%                     toPlot.(params{iParam}).zPerNeurite{iGroup} = zValuesPerNeurite;
%                     clear zValuesPerNeurite
%                 end
%                 
%                      forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).zPerNeurite{:}),1:numel(params),'uniformoutput',0); 
%               
%             end % strcmpi
%             
%         end
%         
%         forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).(ip.Results.diffMetric){:}),1:numel(params),'uniformoutput',0); 
%         forHeatMap =vertcat(forHeatMap{:});
%         projList =  vertcat(toPlot.info.projList{:});
%         names = projList(:,2);
%         
%         % clustergram(forHeatMap,'ColorMap','redbluecmap','RowLabels',params,'ColumnLabels',names);
%         
%         
%         %          colorMap = HeatMapMine(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
%         %       'colormap','redbluecmap','ColumnLabelsRotate',45,'Symmetric',1);
%         colorMap = HeatMap(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
%             'colormap','redbluecmap','ColumnLabelsRotate',45,'DisplayRange',4,'Symmetric',1);
%         %
% %         save('zValuesInd.mat','zValuesInd'); 
%     else % not splitMovie
%% No Partitioning         
        
        for iParam = 1:length(params)
            
            % calculate z-score
            
            control = toPlot.(params{iParam}).dataMat{1}(:);
            control = control(~isnan(control));
            
            % [x,um,std] = zscore(
            
            [x,um,std] = zscore(control);
            % quick sanity check notes: the .dataMat is in the format of
            % r x c such that each column is a cell/neurite (ie observation) the
            % data is pooled for the entire movie - note you checked the nanmean
            % (Default is to take the mean across rows)
            % so currently taking the cellular mean and the cellular std.
            diffMetric = ip.Results.diffMetric;
            switch diffMetric
                case 'perCellZ';
                   % mean = cellfun(@(x) nanmean(nanmean(x)),toPlot.(params{iParam}).dataMat);
                   % Calculate the perNeuriteStat value per movie 
                   % Calculate the mean of the distribution 
                    values = cellfun(@(x) perNeuriteStat(x),toPlot.(params{iParam}).dataMat,'uniformoutput',0); 
                    meanForZ = cellfun(@(x) nanmean(x),values); 
                    % Calculate the standard deviation in the values 
                    stdForZ = cellfun(@(x) nanstd(x),values);
                    
                    % currently have a slight problem with some of the
                    % orientations parameter given complex numbers need to troubleshoot 
                    % quick fix is to just make real here. 
                    meanForZ = real(meanForZ);
                    stdForZ = real(stdForZ);
                    
                    % Standardize the Difference Metric by the Standard
                    % Deviation of the Control Distribution 
                    valuesz = (meanForZ-meanForZ(1))./stdForZ(1);
                    toPlot.(params{iParam}).zScore = real(valuesz);
                                
                case 'pooledZ';
                    mean = cellfun(@(x) nanmean(x(:)),toPlot.(params{iParam}).dataMat);
                    std = cellfun(@(x) nanstd(x(:)),toPlot.(params{iParam}).dataMat);
                    mean = real(mean);
                    std = real(std);
                    valuesz = (mean-mean(1))./std(1);
                    toPlot.(params{iParam}).zScore = real(valuesz);
                    
            end
            %arrayfun(@(x)
        end
        
        
        if ip.Results.clearUnusedFields
            paramsAll = fieldnames(toPlot);
           
                paramsAll = paramsAll(~strcmpi(paramsAll,'info')) ; % keep the info
                paramsRemove = paramsAll(~paramSelectLogic); 
            toPlot = rmfield(toPlot,paramsRemove);
         
        end
        time = clock; 
        %save('compMat.mat','compMat');
        save([outDir filesep 'toPlotHierarchicalCluster'],'toPlot','time');
    end
    % collect data
    forHeatMap = arrayfun(@(x) toPlot.(params{x}).zScore,1:numel(params),'uniformoutput',0);
    forHeatMap = vertcat(forHeatMap{:});
    forHeatMap = forHeatMap(:,2:numel(toPlot.info.names));
    close all
    acutePertPlots = false;
    if acutePertPlots == true
        names = toPlot.info.names;
    else
        names = toPlot.info.names(2:end);
    end
    
    %colorMap = HeatMap(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
      %  'colormap','redbluecmap','ColumnLabelsRotate',45,'DisplayRange',70,'Symmetric',1);
   % save(['colorMap_' ip.Results.diffMetric],'colorMap');
    %   saveas(gcf,'SummaryFigure.png');
    %   saveas(gcf,'SummaryFigure.fig');
    %   saveas(gcf,'SummaryFigure.eps','psc2');
    
    data = mat2dataset(forHeatMap,'ObsNames',params,'varNames',names); 
%     
%     filename  = ['dataToCluster_' ip.Results.diffMetric' '_' ip.Results.perNeuriteStat '.csv']; 
    export(data,'file',[outDir filesep 'dataToCluster.csv']);
    
    if ~isempty(ip.Results.cutOffPValue)
      
        plotDir = [outDir filesep 'IndividualPlotsBeforeCluster']; 
        if ~isdir(plotDir)
            mkdir(plotDir)
        end 
        toPlotMeas =  GCAGroupAnalysisGroupPlots(toPlot,'Measurements',params,...
            'Interactive',false,'OutputDirectory',plotDir,... 
            'plotType','perCell','perNeuriteStat',ip.Results.perNeuriteStat,'FontSize',40,...
            'makePlot',false);
        pValuesCell = arrayfun(@(x) toPlotMeas.(params{x}).pValues.perCell,1:numel(params),'uniformoutput',0); 
        pValues = vertcat(pValuesCell{:});
        pValues = pValues(:,2:end);  % dont' include control 
        filterMat = pValues < ip.Results.cutOffPValue;  
        if ip.Results.filterNonSigFeatures
            test = sum(filterMat,2); 
            filterMat(test==0,:) = []; % filter if all perturbations insignificant. 
            forHeatMap(test==0,:)= []; 
            params(test==0) =[];
        end 
        forHeatMap = forHeatMap.*filterMat; 
    end 
     
    % if filter out non-sig rows for clustering.
    
    
   
    
    
    ClusterObj = clustergram(forHeatMap,'ColorMap','redbluecmap','RowLabels',params,'ColumnLabels',names,'Dendrogram',{'default','default'});
    sortedLabels = get(ClusterObj,'ColumnLabels');
    IDSort = cellfun(@(x) find(strcmpi(x,toPlot.info.names)),sortedLabels);
    IDSort = IDSort-1; % take out the control
    
    sortedLabelsMeas = get(ClusterObj,'RowLabels'); 
    IDRow = cellfun(@(x) find(strcmpi(x,params)),sortedLabelsMeas); 
 
    forHeatMapSort = forHeatMap(IDRow,IDSort); 
    
    %% Add a small bit here to save a imagesc of the sorted Z Scores 
    % (with the clusterObj it is non-intuative how the scaling is
    % performed): also the .eps file of the imagesc images are typically in
    % weird output that doesn't translate position well in illustrator
    % (seems to be saved as some repeated background image) 
    % Therefore it is good to just print the image the approx size you want
    % using imagesc. 
    if ip.Results.plotImagesc
        
        setFigure(172,250,'on')
        imagesc(forHeatMapSort);
        % make sure 'ydir' is normal (typically 'reverse' by default for imagesc)
        set(gca,'ydir','normal');
        
        cmap = get(ClusterObj,'Colormap');
        set(gcf,'Colormap',cmap);
        
        set(gca,'clim',[-3,3])
        helperScreen2png([outDir filesep 'JustHeatMap.png']); 
    end
    %% 
    dataSort = mat2dataset(forHeatMapSort,'ObsNames',sortedLabelsMeas,'varNames',sortedLabels); 
    export(dataSort,'file',[outDir filesep 'dataPostCluster.csv']); 
    
    time = clock; 
    save([outDir filesep 'HierarchCluster.mat'],'ClusterObj','dataSort','data','time');
    
    
     if ip.Results.makePlots 
        % find the variable names that are greater than cut
        % off 
        plotDir = [outDir filesep 'GroupPlots_' ip.Results.diffMetric 'GreaterThan' num2str(ip.Results.plotCutOff)];
        if ~isdir(plotDir)
            mkdir(plotDir)
        end 

        % find only the measurements that are of a
        % significance level that the user wants to target 

        sigValues = abs(forHeatMap)>ip.Results.plotCutOff ; 
        measToPlot = params(sum(sigValues,2)>0); 

        % get the sorted labels from the cluster gram so
        % can plot the groups in that order (simply more
        % intuitive) 
        switch ip.Results.diffMetric
            case 'perCellZ'
                plotType = 'perCell';
            case 'pooledZ'
                plotType = 'pooled';
        end
            
        
        GCAGroupAnalysisGroupPlots(toPlot,'order', sortedLabels,... 
            'Interactive',false,'OutputDirectory',plotDir,... 
            'plotType',plotType,'measurements',measToPlot,'perNeuriteStat',ip.Results.perNeuriteStat,'FontSize',40);
        time = clock; 
        save([plotDir filesep 'timeStamp.mat'],'time');  
     end 
    
    
% end  % split movie
end

