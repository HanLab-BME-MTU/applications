function [ output_args ] = GCAGroupAnalysisCalculateDBIndices( toPlot,varargin)
% GCACalculateDBIndices 
 

dataAll = toPlot.compiled.dataMatZ_nanmeanNoNaN;   
grouping = toPlot.info.groupIDEachFrameNoNaN;

dataAll = dataAll'; % make r = to features and c = to obs

% get all control obs
controlCompare= dataAll(:,grouping==1); 
%controlCompare= dataAll(:,grouping==2); 
% nCond = 8; 
groupsToCompare = [1,2,3,4,5,6,7,8]; 
%groupsToCompare = [2,8]; 
%controlCompare{1} = dataFinalAll{1};
outDir = pwd; 
%  for iCompare = 1:numel(controlCompare)
%        
       
        [DB,pValues,~,CIs] = arrayfun(@(x) bootDiscrimMetrics(dataAll(:,grouping==groupsToCompare(x)),controlCompare,'nrep',10000),1:length(groupsToCompare),'uniformoutput',0);
       
        DB = horzcat(DB{:});
        pValues = vertcat(pValues{:});
        
        setAxis('on');
        [DBSort,idxSort] = sort(DB);
        DBSort = DBSort(2:end);
        pValues = pValues(idxSort,:);
        pValues = pValues(1:end,:);
        
        
        colorsGroup = toPlot.info.color;
        colorsGroup = colorsGroup(groupsToCompare); 
        colorsGroup = colorsGroup(idxSort);
        colorsGroup = colorsGroup(2:end);
        CIsSort = CIs(idxSort); 
        CIsSort = CIsSort(2:end);
        CIsSort = vertcat(CIsSort{:}); 
       
        setAxis('on');
        arrayfun(@(x) scatter(x,DBSort(x),200,colorsGroup{x},'filled'),1:length(DBSort))
        hold on 
        errorbar(1:length(DBSort),DBSort,abs(DBSort-CIsSort(:,1)'),abs(DBSort-CIsSort(:,2)'),'color','k'); 
        arrayfun(@(x) errorbar(x,DBSort(x),abs(DBSort(x)-CIsSort(x,1)),abs(DBSort(x)-CIsSort(x,2)),'color',colorsGroup{x}),1:length(DBSort)); 
        hold on 
        
        namesGroup = toPlot.info.names;
        namesSort = namesGroup(idxSort);
        
        ylabel({'Separation Statistic-'; 'Inverse Davies Bouldin Index'});
        set(gca,'XTick',1:numel(namesSort)-1);
        set(gca,'XTickLabel',namesSort);
        set(gca,'XTickLabelRotation',45); 
        saveas(gcf , [outDir filesep  'InverseDBIndexWithControl' '.png']);
        saveas(gcf , [outDir filesep  'InverseDBIndexWithControl.eps'],'psc2');
        saveas(gcf , [outDir filesep  'InverseDBIndexWithControl.fig']);
        save([outDir filesep 'discrimResults'],'pValues','DB','DBSort','CIsSort');
        close gcf
%%        
% make an example plot in 2D
% idxSort = 7;
% colorsGroup{1} = toPlot.info.color{7}; 
% colorCL = toPlot.info.colorShades{1}(1,:);  
% colorSL = toPlot.info.colorShades{7}(1,:); 
% %colorLine = brewermap(1,'oranges'); 
% colorLine = 'm'; 
% 
% varNames =  fieldnames(toPlot.dataInd);
% idx1=  find(strcmpi(varNames,'filoLengthToVeil'));
% idx2 = find(strcmpi(varNames,'veilStemThickness'));
% 
% % idx1 = find(strcmpi(varNames,'filoDensityAlongVeil'));
% % idx2 = find(strcmpi(varNames,'filoIntensityToVeil')); 
% 
% [DB,pValues,plotInfo ] = bootDiscrimMetrics(dataAll([idx1,idx2],grouping==idxSort(end)),dataAll([idx1,idx2],grouping == 1));
% 
% setAxis('on')
% hold on
% 
% 
% 
% % scatter all the control data
% scatter(dataAll(idx1,grouping==1),dataAll(idx2,grouping == 1),50,...
%     'k','filled','MarkerEdgeColor','w');
% 
% % scatter all the condition data
% scatter(dataAll(idx1,grouping==idxSort(end)),dataAll(idx2,grouping == idxSort(end)),50,...
%     colorsGroup{end},'filled','MarkerEdgeColor','w');
% 
% % draw a line between the center of the two clusters- this is d.
% line([plotInfo.meanS1(1,1),plotInfo.meanS2(1,1)], ...
%     [plotInfo.meanS1(1,2),plotInfo.meanS2(1,2)],'color',colorLine,'Linewidth',5)
% 
% % draw an example of the intra cluster distance
% sdAll1 = plotInfo.intraClustDistS1;
% sdAll2 = plotInfo.intraClustDistS2;
% 
% 
% ex1= find(min(sdAll1 - plotInfo.c1)); % example of mean distance from cluster.
% 
% data1x = dataAll(idx1,grouping==idxSort(end)); 
% data1y = dataAll(idx2,grouping==idxSort(end)); 
% 
% line([data1x(ex1),plotInfo.meanS1(1,1)],[data1y(ex1),plotInfo.meanS1(1,2)], ...
%     'color',colorSL,'Linewidth',5);
% 
% ex2 = find(min(sdAll2-plotInfo.c2));
% dataCx = dataAll(idx1,grouping==1); 
% dataCy = dataAll(idx2,grouping==1); 
% line([dataCx(ex2),plotInfo.meanS2(1,1)],[dataCy(ex2),plotInfo.meanS2(1,2)], ...
%     'color',colorCL,'Linewidth',5);
% 
% 
% % plot the means 
% scatter(plotInfo.meanS1(1,1),plotInfo.meanS1(1,2), 500,colorSL,'d','filled','MarkerEdgeColor',colorsGroup{end},'Linewidth',2); 
% scatter(plotInfo.meanS2(1,1),plotInfo.meanS2(1,2), 500,colorCL,'d','filled','MarkerEdgeColor','k','Linewidth',2); 
% 
% 
% scatter(data1x(ex1),data1y(ex1),500,colorsGroup{end},'filled','MarkerEdgeColor','w','Linewidth',2);
% 
% 
% 
% scatter(dataCx(ex2),dataCy(ex2),500,'k','filled','MarkerEdgeColor','w','LineWidth',2);
% 
% xlabel({'Example Feature 1:' ; 'Filopodia Length To Veil (Normalized) '}); 
% ylabel({'Example Feature 2:' ; 'VeilStemThickness (Normalized)'});
% 
% title({['SepStat =' num2str(DB,3) ' d = ' num2str(plotInfo.d,3)] ; [ 'stdPert = '  num2str(plotInfo.c1,3) ...
%     'stdControl = ' num2str(plotInfo.c2,3)]}); 
%     
%      
%       saveas(gcf,[outDir filesep 'ExampleOfSeparationStatCalc.fig']);
%       saveas(gcf,[outDir filesep 'ExampleOfSeparationStatCalc.eps'],'psc2'); 
%       saveas(gcf,[outDir filesep 'ExampleOfSeparationStatCalc.png']); 
%         
    end 




% end

