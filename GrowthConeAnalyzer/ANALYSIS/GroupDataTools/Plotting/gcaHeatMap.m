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


%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');
ip.addParameter('Interactive',true,@(x) islogical(x));
ip.addParameter('OutputDirectory',pwd,@(x) ischar(x));
ip.addParameter('yLimOff',false,@(x) islogical(x));
ip.addParameter('diffMetric','perCellZ');
ip.addParameter('splitMovie',false);
ip.addParameter('calcZ',true); % adds a z-score to the structure

ip.parse(toPlot,varargin{:});
%%

params = fieldnames(toPlot);
params = params(~strcmpi('info',params)) ;
saveDir = ip.Results.OutputDirectory;
if ip.Results.Interactive
    paramSelect  = listSelectGUI(params,[],'move');
    
    params  = params(paramSelect);
end

if ip.Results.calcZ
    
    if ip.Results.splitMovie
        for iParam = 1:length(params)
            
            if strcmpi(ip.Results.diffMetric,'ks');
                for iGroup = 1:numel(toPlot.info.names);
                    
                    % for now 20151015 lets just make the before and after z score
                    dataMat = toPlot.(params{iParam}).dataMat{iGroup};
                    
                    %ks = arrayfun(@(x) distribTest(dataMat(:,x),dataMat(:,x-1)),2:2:size(dataMat,2));
                    [~,~,ks] = arrayfun(@(x) kstest2(dataMat(:,x),dataMat(:,x-1)),2:2:size(dataMat,2));
                    % get a directionality : for now just just use the
                    % delta in the mean ? 
                    
                    signs =  arrayfun(@(x) sign(nanmedian(dataMat(:,x))-nanmedian(dataMat(:,x-1))),2:2:size(dataMat,2));   %
                    
                    ks = signs.*ks; 
                    
                    toPlot.(params{iParam}).ks{iGroup} = ks;
                    clear ks
                end
            end % strcmpi
            
            
            if strcmpi(ip.Results.diffMetric,'medPermTest');
                for iGroup = 1:numel(toPlot.info.names);
                    
                    % for now 20151015 lets just make the before and after z score
                    dataMat = toPlot.(params{iParam}).dataMat{iGroup};
                    
                    [~,diffMed] = arrayfun(@(x) permTest(dataMat(:,x),dataMat(:,x-1),'CmpFunction',@median),2:2:size(dataMat,2));
                   
                    toPlot.(params{iParam}).permTestMed{iGroup} = diffMed;
                    clear ks
                end
                
                forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).permTestMed{:}),1:numel(params),'uniformoutput',0); 
                
                
            end % strcmpi
            
            
            
            
            
            
            
            if strcmpi(ip.Results.diffMetric,'perCellZ');
                for iGroup = 1:numel(toPlot.info.names);
                    
                    % for now 20151015 lets just make the before and after z score
                    dataMat = toPlot.(params{iParam}).dataMat{iGroup};
                    
                    zValuesPerNeurite = arrayfun(@(x) (nanmean(dataMat(:,x)) - nanmean(dataMat(:,x-1)))./nanstd(dataMat(:,x-1)),2:2:size(dataMat,2));
                    zValuesInd.(params{iParam}){iGroup} = zValuesPerNeurite; 
                    toPlot.(params{iParam}).zPerNeurite{iGroup} = zValuesPerNeurite;
                    clear zValuesPerNeurite
                end
                
                     forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).zPerNeurite{:}),1:numel(params),'uniformoutput',0); 
                
                
            end % strcmpi
            
        end
        
  forHeatMap = arrayfun(@(x) horzcat(toPlot.(params{x}).(ip.Results.diffMetric){:}),1:numel(params),'uniformoutput',0); 
        forHeatMap =vertcat(forHeatMap{:});
        projList =  vertcat(toPlot.info.projList{:});
        names = projList(:,2);
        
        % clustergram(forHeatMap,'ColorMap','redbluecmap','RowLabels',params,'ColumnLabels',names);
        
        
        %          colorMap = HeatMapMine(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
        %       'colormap','redbluecmap','ColumnLabelsRotate',45,'Symmetric',1);
        colorMap = HeatMap(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
            'colormap','redbluecmap','ColumnLabelsRotate',45,'DisplayRange',4,'Symmetric',1);
        %
%         save('zValuesInd.mat','zValuesInd'); 
    else % not splitMovie
        
        
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
                    mean = cellfun(@(x) nanmean(nanmean(x)),toPlot.(params{iParam}).dataMat);
                    std = cellfun(@(x) nanstd(nanstd(x)),toPlot.(params{iParam}).dataMat);
                    mean = real(mean);
                    std = real(std);
                    
                    valuesz = (mean-mean(1))./std(1);
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
        
        %save('compMat.mat','compMat');
        save('toPlotWithStats','toPlot');
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
    
    colorMap = HeatMap(forHeatMap,'RowLabels',params,'ColumnLabels',names,...
        'colormap','redbluecmap','ColumnLabelsRotate',45,'DisplayRange',70,'Symmetric',1);
    save(['colorMap_' ip.Results.diffMetric],'colorMap');
    %   saveas(gcf,'SummaryFigure.png');
    %   saveas(gcf,'SummaryFigure.fig');
    %   saveas(gcf,'SummaryFigure.eps','psc2');
    
    
    clustergram(forHeatMap,'ColorMap','redbluecmap','RowLabels',params,'ColumnLabels',names);
    
end  % split movie
end

