function [ output_args ] = QandDPlotLengthVsCurvInd( toPlot ,varargin)
%QandDPlotLengthVsCurvInd
% Plots the length versus the curvature of individual filopodia and
% calculates the slope of the relationship NOTE This needs to be fixed
% 20160706 NEED TO make sure to filter out the data used in the fit you
% don't really filter it here! 

ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('curvLimits',[0,1.5]);% min/max values of the curvature 
ip.addParameter('lengthLimits',[0,10]);  
ip.addParameter('fitType','exp1'); 

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
        length2 = LC(:,1);
        curv = LC(:,2);
%         curv = LC(:,2)./0.216 ; % convert to pixels 1/um
        % make the scatter
        % convert to form pixels to 1/um
        setAxis('on')
        
        %scatter(length,curv,50,toPlot.info.color{iGroup},'filled');
        lengthClean = length2(~isnan(curv)); 
        curvClean = curv(~isnan(curv));
        hold on        
        %         [fits,errors] = polyfit(lengthClean,curvClean,2);
        
        
        ft = fittype(ip.Results.fitType); 
%         options = fitoptions(ft); 
%         options.lower = [1 -inf ];
%         options.upper = [1, inf]; 
        fitC = fit(lengthClean,curvClean,ft);

        x =linspace(xMin,xMax)';
        y = fitC(x); 
        ci = predint(fitC,x);
     
        dist(:,1) = abs(y-ci(:,1));
        dist(:,2) = abs(y-ci(:,2)); 
        [l,p] = boundedline(x,y,dist);
        set(p,'FaceColor',cMap(i,:)); 
        set(l,'Color',cMap(i,:)); 
        hold on 
      
        scatter(length2,curv,50,toPlot.info.color{iGroup},'filled');
        outlinebounds(l,p);
        plot(x,y,'color',cMap(i,:)); 
%         xfit = linspace(0,max(lengthClean)); 
%         yfit =  polyval(fits,linspace(0,max(lengthClean))); 
%          hold on 
%         plot(xfit,yfit,'color',cMap(i,:)); 
        
        %store plotting infofit
        conInt = confint(fitC); 
        coefs = coeffvalues(fitC); 
        toPlot.fitLenVsCurvPlots{iGroup}{i} = [x' y']; 
        toPlot.fitLenVsCurvConInt{iGroup}{i} = conInt; 
        toPlot.fitLenVsCurv{iGroup}{i} = coefs; 
%         toPlot.fitLenVsCurvError{iGroup}{i} = errors; 
        
        
       % [r,p] = corr(length,curv,'type','Spearman','rows','pairwise');
        
        ylabel('Curvature (um-1)');
        xlabel('Length (um)');
        
        axis([xMin,xMax,yMin,yMax]);
        
        %     values(i,1) = nanmean(LC(:,1)); % length
        %     values(i,2) = nanmean(LC(:,2)); % curv
        %     values(i,3) = c(1,2); % correlation
        %     values(i,4) =
%         values(i,1) = nanmean(LC(:,1));
%         values(i,2) = nanmean(LC(:,2)./0.216);
        hold on
%         scatter(values(i,1),values(i,2),500,'+');
%         values(i,3) = r;
%         values(i,4) = p;
%         title({[IDs{i} ' : ' 'Spearman r= ' num2str(r,3) ]; ...
%             ['PValue = ' num2str(values(i,4),3) ' Mean Length = ' ...
%             num2str(values(i,1),3) 'Mean Curv = ' num2str(values(i,2),3)]});
%  title({IDs{i} ; [' y = ' num2str(coefs(1),2) 'exp' num2str(coefs(2),2) 'x']})
        
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
end % iGroups 
save([outDir filesep 'toPlotCompareLandC.mat'],'toPlot'); 

% if makeplot 

 colors = toPlot.info.color;
 names = toPlot.info.names;
 colorShades = toPlot.info.colorShades;

 % get the number of parameters in the fit 
  n = fittype(ip.Results.fitType); 
  nParams = numel(coeffnames(n)); 
 
for iParam = 1:nParams 
 
for iGroup = 1: numel(toPlot.info.names)
    
    values = cellfun(@(x) x(1,iParam), toPlot.fitLenVsCurv{iGroup});
    values = values'; 
    
    valuesCollect{iGroup,1} = values; 
    
   lower{iGroup,1} = cellfun(@(x,y) abs(x(1,1)-y(2,1)),toPlot.fitLenVsCurv{iGroup},toPlot.fitLenVsCurvConInt{iGroup});
   upper{iGroup,1} = cellfun(@(x,y) abs(x(1,1)-y(2,1)),toPlot.fitLenVsCurv{iGroup},toPlot.fitLenVsCurvConInt{iGroup}); 
end 


coeffAll{iParam} = vertcat(valuesCollect{:}); 

setAxis('on',0.75,20);
dataMat = reformatDataCell(valuesCollect);
h  = notBoxPlot(dataMat);
hold on
 arrayfun(@(i) errorbar(h(i).data.XData,h(i).data.YData,lower{i},upper{i},'Linestyle','none','color',colors{i}),1:numel(names)); 
  %set the colors
  arrayfun(@(i) set(h(i).data,'markerFaceColor', colors{i}),1:numel(names));
  arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(names));
  arrayfun(@(i) set(h(i).mu,'color',colors{i}),1:numel(names));
  arrayfun(@(i) set(h(i).semPtch,'faceColor',colorShades{i}(4,:)),1:numel(names));
  arrayfun(@(i) set(h(i).sdPtch,'faceColor',colorShades{i}(1,:)),1:numel(names));
  groupNames = toPlot.info.names;

ylabel({'Param' num2str(iParam)});
%ylabel(toPlot.(params{iParam}).yLabel);
set(gca,'XTick',1:numel(groupNames));
set(gca,'XTickLabel',groupNames,'FontSize',20);

h1 = get(gcf,'CurrentAxes');
yLim = h1.YLim; % use whatever they used
axis([0.5 numel(names) + 0.5 yLim(1) yLim(2)]);
nameFitType = ip.Results.fitType; 
nameFitType = strrep(nameFitType,'*','_Mult_'); 

saveas(gcf,['CompareFitParam' num2str(iParam) '_' nameFitType '.fig']); 
saveas(gcf,['CompareFitParam' num2str(iParam) '_' nameFitType '.eps'],'psc2');
saveas(gcf,['CompareFitParam' num2str(iParam) '_' nameFitType '.png']); 
close gcf 
end % iParam 
%% Scan 
% minMaxCoeff = cellfun(@(x) [min(x),max(x)],coeffAll,'uniformoutput',0); 
% 
% medianCoeff = cellfun(@(x) [median(x)],coeffAll,'uniformoutput',0); 
% 
% % example of how function varies 
% x =linspace(0,10);
% a = linspace(minMaxCoeff{1}(1),minMaxCoeff{1}(2)); 
% b = medianCoeff{2} ; 
% cMap = jet(length(a)); 
% 
%  % plot holding b constant 
%  
%      y = arrayfun(@(i) a(i)*exp(b(1)*x),1:length(a),'uniformoutput',0);
%      setAxis('on')
%      hold on
%      arrayfun(@(i) plot(x,y{i},'color',cMap(i,:)),1:length(a),'uniformoutput',0);  
%      
%      title({['A from' num2str(a(1),2) 'to' num2str(a(end),2)] ; ['B' num2str(b(1),2) ]}); 
%   saveas(gcf,'scan1.png'); 
%   
%  close gcf
%  %y = arrayfun(@(i,j) a(i)*(exp(b(j)*x)),a,b,'uniformoutput',0); 
% %%
%  
% % example of how function varies 
% b = linspace(minMaxCoeff{2}(1),minMaxCoeff{2}(2)); 
% a = medianCoeff{1} ; 
% cMap = jet(length(b)); 
% 
%  % plot holding b constant 
%  
%      y = arrayfun(@(i) a(1)*exp(b(i)*x),1:length(b),'uniformoutput',0);
%      setAxis('on')
%      hold on 
%      arrayfun(@(i) plot(x,y{i},'color',cMap(i,:)),1:length(b),'uniformoutput',0);  
%      
%      title({['B from' num2str(b(1),2) 'to' num2str(b(end),2)] ; [' A ' num2str(a(1),2)]}); 
% saveas(gcf,scan2.png); 
end



