function [ output_args ] = GCAAnalysisWriteGroupReport(toPlotGroup,figureDirectory)
warning('off','all');
axisFont = 14;
labelFont = 20;
titleFont = 30;




projListC = toPlotGroup.info.projList{1}(:,1);
% Collect all the current parameters.
% perform all the collection and reformatting based on the group
% toPlotGroup = reformatExpressionLevelForPlotting(toPlotGroup,0);
% toPlotGroup = reformatLengthForPlotting(toPlotGroup,0);
% toPlotGroup = reformatOrientForPlotting(toPlotGroup,0);
% toPlotGroup = reformatDensityForPlotting(toPlotGroup,0);
% toPlotGroup = GCAReformatVeilStats(toPlotGroup,0);
% toPlotGroup = reformatBodyShapeForPlotting(toPlotGroup,0);
%write a quick test adding a data mat of normally distributed values of
%different parmeters
%toPlotGroup = extraAddTestBoxPlot(toPlotGroup) ;

%% Get all parameters 

toPlotGroup = GCACollectGroupData(toPlotGroup); 
toPlotGroup = GCAReformatVeilStats(toPlotGroup,0);

[neuriteLengths,deltas] = gcaCollectMultNeuriteLengthsAndDeltas(projListC);
nProjs = size((toPlotGroup.info.projList{1}),1);
cMap = linspecer(nProjs);
[~, idxSortGrowth] = sort(deltas,'ascend');
%% Control Plots: Plot Global Expression Level Estimate: Test For Correlations with Filopodia Length, Neurite Outgrowth, and Protrusion Velocity
%     fsFigure(1);
%     projListAll = vertcat(toPlotGroup.info.projList{:});
%
%     %toPlotGroup = reformatExpressionLevelForPlotting(toPlotGroup,0);
%     names = projListAll(:,2);
%     names = strrep(names,'_',' ');
%     dataMat = toPlotGroup.expressLevel{1};
%     h12 = boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
%     %h1=  boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
%     %set(h1(:),'Linewidth',1);
%     values = dataMat(:);
%     values = values(~isnan(values));
%     minValue = min(values);
%     maxValue = max(values);
%
%     axis([0.5,size(dataMat,2)+0.5, minValue-5, maxValue+5]);
%     set(gca,'FontSize',axisFont);
%     title(['Group: ' toPlotGroup.info.names{1}  ],'FontName','Arial','FontSize',titleFont);
%     ylabel({'Mean Background Subtracted Intensity' ; '(AU)'}','FontName','Arial','FontSize',labelFont);
%
%     % make directory
%
%     cDir = [figureDirectory filesep 'ExpressionLevel'];
%     if ~isdir(cDir);
%         mkdir(cDir);
%     end
%
%     saveas(gcf,[cDir filesep 'ExpressionLevelBoxPlot.fig']);
%     %saveas(gcf,[cDir filesep 'ExpressionLevelBoxPlot.eps'],'psc2');
%     %saveas(gcf,[cDir filesep 'ExpressionLevelBoxPlot.png']);
%
%     snapnow
%     close gcf
%
%     %%% Test and Plot the correlation coefficients
%       % Vs Outgrowth
%       fsFigure(1)
%       set(gca,'FontSize',axisFont);
%       x = nanmedian(toPlotGroup.expressLevel{1},1);
%
%       toTest  = [x' deltas];
%       [r,p] = corrcoef(toTest);
%       title({['Correlation Coefficient' num2str(r(1,2),3)]; ['p-value' num2str(p(1,2),3)]},'FontSize',20);
%       xlabel({'Mean Background Subtracted Intensity' ; '(AU)'},'FontSize',labelFont,'FontName','Arial');
%       ylabel('Net Outgrowth In 10 min','FontSize',labelFont,'FontName','Arial');
%       scatter(x,deltas,100,'k','filled');
%       snapnow
%       if ~isempty(figureDirectory);
%          saveDir =  [figureDirectory filesep 'TestExpressionLevelEffect' ];
%          if ~isdir(saveDir)
%              mkdir(saveDir)
%          end
%          saveas(gcf,[saveDir filesep 'ExpressionEstVsNeuriteOutgrowth.fig']);
%       end
%          close gcf
%
%     % ...Vs Other Parameters
%     paramsToTest{1} = 'Filo';
%     paramsToTest{2} = 'protrusionAnalysis_mednVeloc';
%     labels{1} = 'Filopodia Length (um)';
%     labels{2} = 'Veil Protrusion Velocity (nm/sec)';
%
% for iParam = 1:numel(paramsToTest)
%       fsFigure(1) ;
%       set(gca,'FontSize',axisFont);
%       % for now just collect the mean values might be more correct to collect
%       % the median in the future (mean actually makes the most sense at this
%       % point though)...
%       y = nanmean(toPlotGroup.(paramsToTest{iParam}){1},1);
%       scatter(x',y',100,'k');
%       toTest  = [x' deltas];
%       [r,p] = corrcoef(toTest);
%       title({['Correlation Coefficient' num2str(r(1,2),3)]; ['p-value' num2str(p(1,2),3)]},'FontSize',20);
%
%       ylabel(labels{iParam},'FontSize',labelFont,'FontName','Arial');
%       xlabel({'Mean Background Subtracted Intensity' ; '(AU)'},'FontSize',labelFont,'FontName','Arial');
%       if ~isempty(figureDirectory);
%           saveas(gcf,[saveDir filesep 'ExpressionEstVs' paramsToTest{iParam} '.fig']);
%       end
%       snapnow
%       close gcf
% end
%% Examine Neurite Outgrowth
% % Blue Indicates Lower Outgrowth Cluster and Red Indicates Higher Outgrowth
% % Cluster Based Via Simple K-Means Clustering of Net Outgrowth (Simplest
% % mechanism of clustering for the time being).
% fsFigure(1);
% subplot(1,2,1);
%
% hold on
%
[idx,cCenters]  = kmeans(deltas,2,'replicates',20);
indexMin = find(cCenters == min(cCenters));
indexMax = find(cCenters == max(cCenters));
%  % make it such that 2 is always high and 1 is always the low cluster
sortedIdx = zeros(length(idx),1);
sortedIdx(idx==indexMin) = 1;
sortedIdx(idx==indexMax) = 2;
%
%
%
%
%  neuriteLengthsLow = neuriteLengths(sortedIdx==1);
%  neuriteLengthsHigh = neuriteLengths(sortedIdx==2);
%
% cellfun(@(x)  plotNeuriteOutgrowthInTime(x,'b',1,5,0,[],[],1,1),neuriteLengthsLow);
% hold on
% cellfun(@(x) plotNeuriteOutgrowthInTime(x,'r',1,5,0,[],[],1,1),neuriteLengthsHigh);
%
% subplot(1,2,2);
% % Potentially Make an Optional Modality to Cluster the Data (can potentially have different options...)
% % if runClust == 1
% % sort the indices such that indexMax == 2 and indexMin ==1
%
%
% groupLow = deltas(sortedIdx==1);
% groupHi = deltas(sortedIdx==2);
% scatter(ones(length(groupHi),1),groupHi,100,'r','filled');
% hold on
% scatter(ones(length(groupLow),1),groupLow,100,'b','filled');
% text(1,cCenters(indexMin),['Low Outgrowth = ' num2str(cCenters(indexMin,1),3) ]);
% text(1,cCenters(indexMax),['High Outgrowth = ' num2str(cCenters(indexMax,1),3)]);
%
%
% %scatter(ones(length(deltas),1),deltas,100 ,'k','filled','MarkerEdge','w');
%
% axis([0.5,1.5,-10,25]);
% ylabel({'Delta Neurite Outgrowth (um)' ; 'in 10 mins' },'FontSize',20,'FontName','Arial');
%
% set(gca,'XTick',1);
% set(gca,'XTickLabel',[ toPlotGroup.info.names{1}  ]);
% set(gca,'FontName','Arial','FontSize',18);
% snapnow
%
% if ~isempty(figureDirectory)
%       saveas(gcf,[figureDirectory filesep 'neuriteOutgrowthCompiled.fig']);
% end
%
% close gcf
% Outgrowth Plots Color-Coded by Identifier
fsFigure(1);
subplot(1,2,1);
% get colormap


hold on
%nProjs = size(toPlotGroup.info.projList{1},1);

%cMap = linspecer(nProjs);
neuriteLengths = neuriteLengths(idxSortGrowth);
hold on
arrayfun(@(i) plotNeuriteOutgrowthInTime(neuriteLengths{i},cMap(i,:),1,5,0,[],[],1,1),1:nProjs);


subplot(1,2,2);
% Potentially Make an Optional Modality to Cluster the Data (can potentially have different options...)
% if runClust == 1
% sort the indices such that indexMax == 2 and indexMin ==1
%hold on
%arrayfun(@(i) plotNeuriteOutgrowthInTime(neuriteLengths{i},cMap(i,:),1,5,0,[],[],1,1,'vel'),1:nProjs);

% text(1,cCenters(indexMin),['Low Outgrowth = ' num2str(cCenters(indexMin,1),3) ]);
% text(1,cCenters(indexMax),['High Outgrowth = ' num2str(cCenters(indexMax,1),3)]);

deltasSort = deltas(idxSortGrowth)./10;
hold on
arrayfun(@(i) scatter(1,deltasSort(i),300 ,cMap(i,:),'filled','MarkerEdge','w'),1:length(deltasSort));




snapnow

if ~isempty(figureDirectory)
    saveas(gcf,[figureDirectory filesep 'neuriteOutgrowthCompiledColorByNeurite.fig']);
end

close gcf
%% Visualize Per Cell Parameter Distributions
% Before you can pool the data one needs to test if each of the cell
% measurements is exhibiting a similar distribution:
% It also helps to understand the typical distribution observed for each
% parameter: this will inevitably help pick which statistical tests are
% most appropriate: Currently each per movie distribution
% is tested against the null hypothesis that it comes from a
% ...normal,extremevalue,log normal,weibull
% using the anderson darling test. If the per movie data fails to reject the null
% hypothesis for one of the stat tests a text indicating this is shown.
% The null hypothesis here is that the data sampled is from
% the respective test distribution.

projListAll = vertcat(toPlotGroup.info.projList{:}) ;
names = projListAll(:,2);
names = strrep(names,'_',' ');
%% Test adTest Screen
fsFigure(1);
toPlotGroup = extraAddTestBoxPlot(toPlotGroup); 

dataMat = toPlotGroup.testDist{1};
gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,labelFont,titleFont);
saveas(gcf,[figureDirectory filesep 'test.fig']); 
snapnow
close gcf

%% Filopodia Length (um) Per Cell Distributions (Traditional Definition: 'External')
fsFigure(1);
dataMat = toPlotGroup.filoLengthToVeil{1};



gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,labelFont,titleFont);
snapnow
if ~isempty([figureDirectory]);
    saveas(gcf,[figureDirectory filesep 'Filopodia_Length_Distributions.fig']) ;
end
close gcf
%% Filopodia Orientation Relative to the Veil (Degrees 0-90) Per Cell Distributions
% fsFigure(1);
% dataMat = toPlotGroup.FiloOrient{1};
% dataMat = dataMat(:,1:end-1); % need to rerun 21 quck fix for now
% namesC = names(1:end-1);
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,namesC,axisFont,labelFont,titleFont);
% if ~isempty(figureDirectory);
%    saveas(gcf,[figureDirectory filesep 'Filopodia_Orientation_Distributions.fig']);
% end
% snapnow
% close gcf
%% Filodia Density (Number of Filopodia Per 10 um of Veil
% Note This Includes Filopodia Attached to Veil Exhibiting Branches
% but not the Branch Structures Themselves
% fsFigure(1);
% dataMat = toPlotGroup.FiloDensity{1};
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,labelFont,titleFont);
% snapnow
% close gcf
%% Veil Protrusion Velocity Per Cell Distributions
% fsFigure(1);
% dataMat = toPlotGroup.protrusionAnalysis_mednVeloc{1};
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,labelFont,titleFont);
% if ~isempty(figureDirectory)
%     saveas(gcf,[figureDirectory filesep 'Filopodia_Protrusion_Velocity_Distributions.fig']) ;
% end
%
% snapnow
% close gcf
%% Veil Retraction Velocity Per Cell Distributions
% fsFigure(1);
% dataMat = toPlotGroup.retractionAnalysis_mednVeloc{1};
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,titleFont,labelFont);
% if ~isempty(figureDirectory)
%     saveas(gcf,[figureDirectory filesep 'Filopodia_Retraction_Velocity_Distributions.fig']);
% end
% snapnow
% close gcf
%% Veil Protrusion Persistence Time Per Cell Distributions
% fsFigure(1);
% dataMat = toPlotGroup.protrusionAnalysis_persTime{1};
% 
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,titleFont,labelFont);
% 
% if ~isempty(figureDirectory)
%     saveas(gcf,[figureDirectory filesep 'Filopodia_Protrusion_Persistence_Distributions.fig']) ;
% end
% 
% 
% snapnow
% close gcf
%% Veil Retraction Persistence Time PerCell Distributions
% fsFigure(1);
% dataMat = toPlotGroup.retractionAnalysis_persTime{1};
% gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,names,axisFont,labelFont,titleFont);
% snapnow
% close gcf
%% Screen for Correlations With Neurite Outgrowth (still a bit rough...)
[results,dataSetArray] = GCAAnalysisScreenCorrelations(deltas,toPlotGroup);

% save the dataSetarray 
export(dataSetArray,'File',[figureDirectory filesep  'DataSetMeanWholeMovie.csv'],'Delimiter',',');

mdl = fitlm(dataSetArray); 
save([figureDirectory filesep 'ModelObjectWholeMovie.mat'],'mdl'); 
save([figureDirectory filesep 'DataSetArray.mat'],'dataSetArray'); 

corrplotMine(dataSetArray,'testR','on'); 
saveas(gcf,[figureDirectory filesep 'CorrelationsAll.fig']); 
close gcf

% results of the correlations will be a simple %
%if ~isempty(resultsCorrScreen.Hit{1,1});
% Get the Indices of the Hits
idxHit = find(arrayfun(@(x) results(x).p(1,2) <0.05,1:length(results)));
idxBorderline = find(arrayfun(@(x) (results(x).p(1,2) >0.05 & results(x).p(1,2) < 0.1 ),1:length(results)));
idxNoCorr = find(arrayfun(@(x) (results(x).p(1,2) > 0.1),1:length(results)));
%% Potentially Significant Correlation With Outgrowth

sigCorr = [figureDirectory filesep 'SigCorr']; 
if ~isdir(sigCorr) 
    mkdir(sigCorr); 
end 

% Make the Correlation Plots
for iParam =1: length(idxHit)
    
    fsFigure(1);
    %  subplot(2,2,1:2);
    resultsC = results(idxHit(iParam));
    hold on
    sorty = resultsC.values(idxSortGrowth,1); % put outgrowth on the y axis
    
    sortx = resultsC.values(idxSortGrowth,2);
    arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity in 10 min','FontSize',labelFont+7,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', labelFont+7,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',axisFont+7);
    
    
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[sigCorr filesep resultsC.name 'Corr.png']); 
    close gcf
    fsFigure(1);
    %subplot(2,2,3:4);
    % sort by outgrowth
    
    
    dataMat = toPlotGroup.(resultsC.name){1};
    if strcmpi(resultsC.name,'retractionAnalysis_persTime') || strcmpi(resultsC.name,'protrusionAnalysis_persTime')
        nMovies = size(dataMat,2);
        % get the 75th percentile of each movie for the persistence time
        % /retraction time % set anything lower to this value to zero.
        prctile75 = arrayfun(@(i) prctile(dataMat(:,i),75),1:nMovies);
        toRemove = arrayfun(@(i) dataMat(:,i)<prctile75(i),1:nMovies,'uniformoutput',0);
        toRemoveMat = horzcat(toRemove{:});
        %dataMat75 = dataMat;
        dataMat(toRemoveMat)=NaN; % set to NaN;
        
    end
    
    
    
    dataMatSort = dataMat(:,idxSortGrowth);
    % sort the cluster indices (I know this is confusing should rename)
    % Make the boxplot of the distrubtions again sorted by outgrowth
    % and colored by cluster.
    growthClusterIdx =  sortedIdx(idxSortGrowth);
    namesSorted = names(idxSortGrowth);
    % h12 = boxplot(dataMatSort,'color',['b','r'],'colorGroup',growthClusterIdx,...
    %  'notch','on','outlierSize',1,'Labels',namesSorted,'labelorientation','inline');
    
    h12 = boxplot(dataMatSort,'color',cMap,...
        'notch','on','outlierSize',1,'Labels',namesSorted,'labelorientation','inline');
    
    values = dataMatSort(:);
    values = values(~isnan(values));
    minValue = min(values);
    maxValue = max(values);
    
    axis([0.5,size(dataMatSort,2)+0.5, minValue, maxValue]);
    set(gca,'FontSize',axisFont);
    %title(['Group: ' toPlotGroup.info.names{1}  ],'FontName','Arial','FontSize',titleFont);
    saveas(gcf,[resultsC.name 'boxplot.fig']);
    saveas(gcf,[resultsC.name 'boxplot.eps'],'psc2');
    snapnow
    close gcf
    
end

%% BorderLine Correlations With Outgrowh
bordDir = [figureDirectory filesep 'BorderLineCorr']; 
if ~isdir(bordDir) 
    mkdir(bordDir); 
end 

for iParam =1: length(idxBorderline)
    
    fsFigure(1);
    % subplot(2,2,1:2);
    resultsC = results(idxBorderline(iParam));
    
    sorty = resultsC.values(idxSortGrowth,1);
    
    sortx = resultsC.values(idxSortGrowth,2);
    hold on
    arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity in 10 min','FontSize',labelFont+7,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', labelFont+7,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',axisFont+7);
    
    
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[bordDir filesep resultsC.name 'Corr.png']); 
    
    
    
    %  scatter(resultsC.values(:,1),resultsC.values(:,2),100,'k','filled');
    % title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',14);
    %xlabel('Neurite Net Velocity in 10 min','FontSize',labelFont-2,'FontName','Arial');
    %ylabel(strrep(resultsC.name,'_',' '),'FontSize', labelFont-2,'FontName','Arial');
    % subplot(2,2,3:4);
    % sort by outgrowth
    %[~, idxSortGrowth] = sort(deltas,'ascend');
    close gcf
    fsFigure(1);
    dataMat = toPlotGroup.(resultsC.name){1};
    dataMatSort = dataMat(:,idxSortGrowth);
    % sort the cluster indices (I know this is confusing should rename)
    % Make the boxplot of the distrubtions again sorted by outgrowth
    % and colored by cluster.
    growthClusterIdx =  sortedIdx(idxSortGrowth);
    namesSorted = names(idxSortGrowth);
    h12 = boxplot(dataMatSort,'color',['b','r'],'colorGroup',growthClusterIdx,...
        'notch','on','outlierSize',1,'Labels',namesSorted,'labelorientation','inline');
    
    values = dataMatSort(:);
    values = values(~isnan(values));
    minValue = min(values);
    maxValue = max(values);
    
    axis([0.5,size(dataMatSort,2)+0.5, minValue, maxValue]);
    set(gca,'FontSize',axisFont);
    %title(['Group: ' toPlotGroup.info.names{1}  ],'FontName','Arial','FontSize',titleFont);
    
    
    snapnow
    
    close gcf
    
end
%%  No Correlation with Outgrowth

noCorrDir = [figureDirectory filesep 'noCorrelation']; 
if ~isdir(noCorrDir) 
    mkdir(noCorrDir)
end 

for iParam = 1:length(idxNoCorr)
    fsFigure(1);
    % subplot(2,2,1:2);
    resultsC = results(idxNoCorr(iParam));
    
    sorty = resultsC.values(idxSortGrowth,1);
    
    sortx = resultsC.values(idxSortGrowth,2);
    hold on
    arrayfun(@(i) scatter(sortx(i),sorty(i),500,cMap(i,:),'filled'),1:nProjs);
    title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',20);
    ylabel('Neurite Net Velocity in 10 min','FontSize',labelFont+7,'FontName','Arial');
    xlabel(strrep(resultsC.name,'_',' '),'FontSize', labelFont+7,'FontName','Arial');
    set(gca,'FontName','Arial','FontSize',axisFont+7);
    
    
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.fig']);
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.eps'],'psc2');
    saveas(gcf,[noCorrDir filesep resultsC.name 'Corr.png']); 
    
    close gcf
    fsFigure(1);
    %scatter(resultsC.values(:,1),resultsC.values(:,2),100,'k','filled');
    %title({['Correlation Coefficient' num2str(resultsC.r(1,2),3)]; ['p-value' num2str(resultsC.p(1,2),3)]},'FontSize',14);
    %xlabel('Neurite Net Velocity (um/min)','FontSize',labelFont,'FontName','Arial');
    %ylabel(strrep(resultsC.name,'_',' '),'FontSize', labelFont,'FontName','Arial');
    
    % subplot(2,2,3:4);
    % sort by outgrowth
    %[~, idxSortGrowth] = sort(deltas,'ascend');
    
    dataMat = toPlotGroup.(resultsC.name){1};
    dataMatSort = dataMat(:,idxSortGrowth);
    % sort the cluster indices (I know this is confusing should rename)
    % Make the boxplot of the distrubtions again sorted by outgrowth
    % and colored by cluster.
    growthClusterIdx =  sortedIdx(idxSortGrowth);
    namesSorted = names(idxSortGrowth);
    h12 = boxplot(dataMatSort,'color',['b','r'],'colorGroup',growthClusterIdx,...
        'notch','on','outlierSize',1,'Labels',namesSorted,'labelorientation','inline');
    
    values = dataMatSort(:);
    values = values(~isnan(values));
    minValue = min(values);
    maxValue = max(values);
    
    axis([0.5,size(dataMatSort,2)+0.5, minValue, maxValue]);
    set(gca,'FontSize',axisFont);
    
    
    %subplot(2,2,3:4);
    % sort by outgrowth
    %[~, idxSortGrowth] = sort(deltas,'ascend');
    
    dataMat = toPlotGroup.(resultsC.name){1};
    dataMatSort = dataMat(:,idxSortGrowth);
    % sort the cluster indices (I know this is confusing should rename)
    % Make the boxplot of the distrubtions again sorted by outgrowth
    % and colored by cluster.
    growthClusterIdx =  sortedIdx(idxSortGrowth);
    namesSorted = names(idxSortGrowth);
    h12 = boxplot(dataMatSort,'color',['b','r'],'colorGroup',growthClusterIdx,...
        'notch','on','outlierSize',1,'Labels',namesSorted,'labelorientation','inline');
    
    values = dataMatSort(:);
    values = values(~isnan(values));
    minValue = min(values);
    maxValue = max(values);
    
    axis([0.5,size(dataMatSort,2)+0.5, minValue, maxValue]);
    set(gca,'FontSize',axisFont);
    
    snapnow
    close gcf
end


end

