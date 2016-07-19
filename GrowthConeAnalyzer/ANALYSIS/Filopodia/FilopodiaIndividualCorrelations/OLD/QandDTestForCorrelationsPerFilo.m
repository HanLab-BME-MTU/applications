function [ output_args ] = QandDTestForCorrelationsPerFilo( toPlot ,varargin)
%QandDPlotLengthVsCurvInd
% Plots the length versus the curvature of individual filopodia and
% calculates the slope of the relationship 

ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

% ip.addParameter('curvLimits',[0,1.5]);% min/max values of the curvature 
% ip.addParameter('lengthLimits',[0,10]);  


ip.addParameter('Interactive',true,@(x) islogical(x)); 
ip.addParameter('Measurements',[]); % cell of measurements to be plotted 


ip.addParameter('getRepScatter',true); 

ip.addParameter('makeScatterPlots',true,@(x) islogical(x)); 

ip.addParameter('collectRepScatters',true); 
ip.addParameter('collectDirectory', [pwd filesep 'RepresentativeScatters']); 


ip.parse(varargin{:});

%%


if ~isempty(ip.Results.Measurements);
    params = ip.Results.Measurements;
else
    params = fieldnames(toPlot);
    params = params(~strcmpi('info',params)) ;
    
    if ip.Results.Interactive
        paramSelect  = listSelectGUI(params,[],'move');
        params  = params(paramSelect);
    end
end

nameFieldCompareR = ['r_' params{1} '_vs_' params{2}]; 
nameFieldCompareP = ['p_' params{1} '_vs_' params{2}]; 

nGroups = numel(toPlot.info.names);
nProjsCon = size(toPlot.info.projList{1},1); 
cMapCon = [toPlot.info.colorShades{1} ; toPlot.info.colorShades{1}];
outDir = pwd;

% 
if ip.Results.getRepScatter 
    if ~isdir(ip.Results.collectDirectory)
        mkdir(ip.Results.collectDirectory); 
    end    
end 

% get minimum and maximum 
    dataAllX = cellfun(@(x) x(:),toPlot.(params{1}).dataMat,'uniformoutput',0); 
    dataAllX = vertcat(dataAllX{:}); 
    xMin = min(dataAllX); 
    xMax = max(dataAllX); 
    
    dataAllY = cellfun(@(x) x(:), toPlot.(params{2}).dataMat,'uniformoutput',0); 
    dataAllY = vertcat(dataAllY{:}); 
    yMin = min(dataAllY); 
    yMax = max(dataAllY); 
    
%     dataAllX = horzcat(toPlot.params{1}.dataMat{:}); 
%     xMax =  max(dataAll(:)); 
%     xMin = min(dataAll(:)); 
   
   
    % get the axis coordinates
%     xMax = max(dataAll(:,1));
%     xMin = min(currentMatAll(:,1));
%     
%     yMax = max(currentMatAll(:,2));
%     yMin = min(currentMatAll(:,2));


%% Start 
for iGroup = 1:nGroups
    
    %if ip.Results.makeScatterPlots 
        figureDir = [outDir filesep 'ScatterPlots_' toPlot.info.names{iGroup}];
        
        if ~isdir(figureDir)
            mkdir(figureDir)
        end
    %end
    n = size(toPlot.(params{1}).dataMat{iGroup},2); 
    
    [currentMat] = arrayfun(@(x) [toPlot.(params{1}).dataMat{iGroup}(:,x), toPlot.(params{2}).dataMat{iGroup}(:,x)],1:n,'uniformoutput',0);
    
  
    IDs = toPlot.info.projList{iGroup}(:,2);
    
    IDs = cellfun(@(x) strrep(x,'_',' '),IDs,'uniformoutput',0);
    
  
    
    colors = toPlot.info.colorShades{iGroup}; 
    cMap = [colors; colors]; 
    
    for i = 1:numel(currentMat)
        
        LC = currentMat{i};
      
        % convert to form pixels to 1/um
      
      
        % test for nan 
        idxNoNan = ~isnan(LC(:,1)) & ~isnan(LC(:,2)); 
        LC = LC(idxNoNan,:); 
        [r(1,i),p(1,i)] = corr(LC(:,1),LC(:,2),'type','Spearman','rows','pairwise');
        
        
        
        if ip.Results.makeScatterPlots
            
            setAxis('off')
            
            scatter(LC(:,1),LC(:,2),50,toPlot.info.color{iGroup},'filled');
            
            
            ylabel(params{2});
            xlabel(params{1});
            
            
            values(i,1) = nanmean(LC(:,1));
            values(i,2) = nanmean(LC(:,2));
            hold on
            scatter(values(i,1),values(i,2),500,'+');
            values(i,3) = r(1,i);
            values(i,4) = p(1,i);
            title({[IDs{i} ' : ' 'Spearman r= ' num2str(r(1,i),3) ]; ...
                ['PValue = ' num2str(values(i,4),3) params{1} ' = ' ...
                num2str(values(i,1),3) params{2} ' = ' num2str(values(i,2),3)]});
            axis([xMin,xMax,yMin,yMax]);
            saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'Ind' params{1} 'vs' params{2} '.fig'])
            saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'Ind' params{1} 'vs' params{2} '.eps'],'psc2')
            saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'Ind' params{1} 'vs' params{2} '.png']);
            close gcf
        end % if makeScatterPlots
    end
    toPlot.(nameFieldCompareR).dataMat{iGroup} = r; 
    toPlot.(nameFieldCompareP).dataMat{iGroup} = p; 
    percentP = sum(p<0.05)./length(p);
    toPlot.(nameFieldCompareP).fractionCellsSig{iGroup} = percentP; 
    
%     toPlot.(nameFieldCompareP).
    
    
    % find the project closest to the mean r and move that file 
    if ip.Results.getRepScatter 
       delt =  abs(r - nanmean(r));
       idx =  delt == min(delt); 
       from  = [figureDir filesep toPlot.info.projList{iGroup}{idx,2} 'Ind' params{1} 'vs' params{2} '.png']; 
       to = [ip.Results.collectDirectory filesep toPlot.info.projList{iGroup}{idx,2} 'Ind' params{1} 'vs' params{2} '.png']; 
       copyfile(from,to);  
    end 
    
    clear r p
   
    clear values 

end
save([outDir filesep 'toPlotCompare' params{1} 'vs' params{2} '.mat'],'toPlot'); 
