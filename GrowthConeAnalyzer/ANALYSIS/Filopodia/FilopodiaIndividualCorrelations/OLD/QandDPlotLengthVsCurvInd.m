function [ output_args ] = QandDPlotLengthVsCurvInd( toPlot ,varargin)
%QandDPlotLengthVsCurvInd
% Plots the length versus the curvature of individual filopodia and
% calculates the slope of the relationship 

ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('curvLimits',[0,1.5]);% min/max values of the curvature 
ip.addParameter('lengthLimits',[0,10]);  

ip.parse(varargin{:});

%%


nGroups = numel(toPlot.info.names);
nProjsCon = size(toPlot.info.projList{1},1); 
cMapCon = [toPlot.info.colorShades{1} ; toPlot.info.colorShades{1}];
outDir = pwd;
compareDir = [outDir filesep 'FitCompareToControl']; 
if ~isdir(compareDir) 
    mkdir(compareDir) 
end 
for iGroup = 1:nGroups
    
    
    figureDir = [outDir filesep 'ScatterPlots_' toPlot.info.names{iGroup}];
    
    if ~isdir(figureDir)
        mkdir(figureDir)
    end
    
    n = size(toPlot.filoCurvature.dataMat{iGroup},2);
    nProjs = size(toPlot.info.projList{iGroup},1); 
    
    
    [lengthCurv] = arrayfun(@(x) [toPlot.filoLengthToVeil.dataMat{iGroup}(:,x), toPlot.filoCurvature.dataMat{iGroup}(:,x)],1:n,'uniformoutput',0);
    
    
    
    setAxis('off')
    hold on
    
    IDs = toPlot.info.projList{iGroup}(:,2);
    
    IDs = cellfun(@(x) strrep(x,'_',' '),IDs,'uniformoutput',0);
    
    
    lengthCurvAll = vertcat(lengthCurv{:});
    lengthCurvAll(:,2)= lengthCurvAll(:,2)./0.216 ; % convert to pixels 1/um
    
    idx1 = (lengthCurvAll(:,2) > ip.Results.curvLimits(1) & lengthCurvAll(:,2) < ip.Results.curvLimits(2));
    idx2 = (lengthCurvAll(:,1) > ip.Results.lengthLimits(1) & lengthCurvAll(:,1) < ip.Results.lengthLimits(2));
    
    lengthCurvAll = lengthCurvAll(idx1 & idx2,:);
    
    xMax = max(lengthCurvAll(:,1));
    xMin = min(lengthCurvAll(:,1));
    
    yMax = max(lengthCurvAll(:,2));
    yMin = min(lengthCurvAll(:,2));
    
    % convert
%     yMax = yMax/0.216; % curvature in 1/um - fix eventually earlier
%     yMin = yMin/0.216;
%     
    colors = toPlot.info.colorShades{iGroup}; 
    cMap = [colors; colors]; 
    
    for i = 1:numel(lengthCurv)
        
        LC = lengthCurv{i};
        length = LC(:,1);
        curv = LC(:,2);
%         curv = LC(:,2)./0.216 ; % convert to pixels 1/um
        % make the scatter
        % convert to form pixels to 1/um
        setAxis('on')
        
        scatter(length,curv,50,toPlot.info.color{iGroup},'filled');
        lengthClean = length(~isnan(curv)); 
        curvClean = curv(~isnan(curv));
        
        [fits,errors] = polyfit(lengthClean,curvClean,2); 
        xfit = linspace(0,max(lengthClean)); 
        yfit =  polyval(fits,linspace(0,max(lengthClean))); 
         hold on 
        plot(xfit,yfit,'color',cMap(i,:)); 
        
        %store plotting info
        toPlot.fitLenVsCurvPlots{iGroup}{i} = [xfit' yfit']; 
        toPlot.fitLenVsCurv{iGroup}{i} = fits; 
        toPlot.fitLenVsCurvError{iGroup}{i} = errors; 
        
        
        [r,p] = corr(length,curv,'type','Spearman','rows','pairwise');
        
        ylabel('Curvature (um-1)');
        xlabel('Length (um)');
        
        axis([xMin,xMax,yMin,yMax]);
        
        %     values(i,1) = nanmean(LC(:,1)); % length
        %     values(i,2) = nanmean(LC(:,2)); % curv
        %     values(i,3) = c(1,2); % correlation
        %     values(i,4) =
        values(i,1) = nanmean(LC(:,1));
        values(i,2) = nanmean(LC(:,2)./0.216);
        hold on
        scatter(values(i,1),values(i,2),500,'+');
        values(i,3) = r;
        values(i,4) = p;
        title({[IDs{i} ' : ' 'Spearman r= ' num2str(r,3) ]; ...
            ['PValue = ' num2str(values(i,4),3) ' Mean Length = ' ...
            num2str(values(i,1),3) 'Mean Curv = ' num2str(values(i,2),3)]});
        
        saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'IndLenVsIndCurv.fig'])
        saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'IndLenVsIndCurv.eps'],'psc2')
        saveas(gcf,[figureDir filesep toPlot.info.projList{iGroup}{i,2} 'IndLenVsIndCurv.png']);
        close gcf
        
    end
    % Combine the fits 
    setAxis('on') 
    hold on
    % Always plot control 
    arrayfun(@(x) plot(toPlot.fitLenVsCurvPlots{1}{x}(:,1),toPlot.fitLenVsCurvPlots{1}{x}(:,2) ...
        ,'color',cMapCon(x,:)),1:nProjsCon); 
    hold on 
    if iGroup >1 
         arrayfun(@(x) plot(toPlot.fitLenVsCurvPlots{iGroup}{x}(:,1),toPlot.fitLenVsCurvPlots{iGroup}{x}(:,2) ...
        ,'color',cMap(x,:),'Linewidth',2),1:nProjs);
    end 
    
    xlabel('Filopodia Length (um)'); 
    ylabel({'Filopodia Maximum' ;' Curvature (um-1)'}); 
    title(toPlot.info.names{iGroup}); 
    axis([0 30 0 5]);  
    saveas(gcf,[compareDir filesep 'FitsPerGroup_' toPlot.info.names{iGroup} '.fig']); 
    saveas(gcf,[compareDir filesep 'FitsPerGroup_' toPlot.info.names{iGroup} '.png']);
    close gcf
  
   % save([outDir filesep 'values'],'values');
end
save([outDir filesep 'toPlotCompareLandC.mat'],'toPlot'); 
