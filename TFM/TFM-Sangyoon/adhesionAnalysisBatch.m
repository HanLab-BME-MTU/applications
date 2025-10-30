%% open necessary MLs
clear

[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,~,MLdirect]=chooseSelectedFolders;
nameList=groupNames';
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
% rootAnalysis = pathAnalysisAll{1};
summaryPath = [rootAnalysis '/AnalysisSummary_AdhesionStatic' specificName];
ii=0;
while exist(summaryPath, 'dir')
    ii=ii+1;
    summaryPath = [rootAnalysis '/AnalysisSummary_AdhesionStatic' specificName num2str(ii)];
end
figPath = [summaryPath '/Figs'];
mkdir(figPath)
dataPath = [summaryPath '/Data'];
mkdir(dataPath)
save([rootAnalysis filesep 'selectedFolders' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end

%% Run detectMovieNascentAdhesion for each list
N=zeros(numConditions,1);
NAdensity = cell(numConditions,1);
FAarea = cell(numConditions,1);
FAlength = cell(numConditions,1);
FAdensity = cell(numConditions,1);
FAdensityPeri = cell(numConditions,1);
FAdensityInside = cell(numConditions,1);
NAstructGroup= cell(numConditions,1);
FAstructGroup= cell(numConditions,1);
forceAllGroup= cell(numConditions,1);
ampTheOtherAllGroup= cell(numConditions,1);
bandNA = 5; % 2 um from the edge

sampleMovie = MLAll(1).movies_{1};
numChan = numel(sampleMovie.channels_);
if numChan==1
    iPax = numChan;
else
    iPax = input('Enter adhesion channel number (1 or 2 ..): ');
end

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);

    curNAdensity = zeros(N(ii),1);
    curFAarea = cell(N(ii),1);
    curFAlength = cell(N(ii),1);
    curFAdensity = zeros(N(ii),1);
    curFAdensityPeri = zeros(N(ii),1);
    curFAdensityInside = zeros(N(ii),1);
    curML = MLAll(ii);
    
    % Combine data per each condition (1,2,3,4 for 3.4, 18, 100, 100Y, respectively)
    for k=1:N(ii)
        [curNAstruct,curFAstruct] = detectMovieNascentAdhesion(curML.movies_{k},bandNA,iPax);
        if numel(curFAstruct)>1
            curNAdensity(k) = mean(arrayfun(@(x) x.NAdensity,curNAstruct));
            curFAarea{k} = cell2mat(arrayfun(@(x) x.area,curFAstruct,'unif',false));
            curFAlength{k} = cell2mat(arrayfun(@(x) x.length,curFAstruct,'unif',false));
            curFAdensity(k) = mean(arrayfun(@(x) x.FAdensity,curFAstruct));
            curFAdensityPeri(k) = mean(arrayfun(@(x) x.FAdensityPeri,curFAstruct));
            curFAdensityInside(k) = mean(arrayfun(@(x) x.FAdensityInside,curFAstruct));
            NAstruct(k) = curNAstruct(end);
            try
                FAstruct(k) = curFAstruct(end);
            catch
                % This time they have in dissimilar structure
                FAstruct(k)=FAstruct(k-1); % I assume k is already more than 1
                fieldList = fieldnames(curFAstruct(end));
                for pp=1:numel(fieldList)
                    FAstruct(k).(fieldList{pp})=curFAstruct(end).(fieldList{pp});
                end
            end
        else %numel(curFAstruct)==1
            curNAdensity(k) = curNAstruct.NAdensity;
            curFAarea{k} = curFAstruct.area;
            curFAlength{k} = curFAstruct.length;
            curFAdensity(k) = curFAstruct.FAdensity;
            NAstruct(k) = curNAstruct;
            if isfield(FAstruct,'ampTheOther') && ~isfield(curFAstruct,'ampTheOther')
                FAstruct = rmfield(FAstruct,'ampTheOther');
            end
            FAstruct = orderfields(FAstruct, curFAstruct);
            FAstruct(k) = curFAstruct;
            % FAdensity at cell periphery
            curFAdensityPeri(k) = mean(arrayfun(@(x) x.FAdensityPeri,curFAstruct));
            % FAdensity at cell center
            curFAdensityInside(k) = mean(arrayfun(@(x) x.FAdensityInside,curFAstruct));
        end
        % Reading potential forces
        iTFM = curML.movies_{k}.getPackageIndex('TFMPackage');
        if ~isempty(iTFM)
            outputFilePath = [curML.movies_{k}.outputDirectory_ filesep 'Adhesion Quantification'];
            curDataPath= [outputFilePath filesep 'data'];
            forceAll(k) = load([curDataPath filesep 'forcePerAdhesionType.mat'],'forceFA', 'forceFC', 'forceNA', 'forceBGinCell', 'forceBGoutCell');
%         else
%             forceAll(k)=[];
        end
        % Reading the other channels's amplitudes
        nChan = numel(curML.movies_{k}.channels_);
        if nChan>2
            outputFilePath = [curML.movies_{k}.outputDirectory_ filesep 'Adhesion Quantification'];
            curDataPath= [outputFilePath filesep 'data'];
            ampTheOtherAll(k) = load([curDataPath filesep 'ampTheOtherPerAdhesionType.mat']);
%         else
%             ampTheOtherAll(k).ampTheOtherFA=NaN;
%             ampTheOtherAll(k).ampTheOtherFC=NaN;
%             ampTheOtherAll(k).ampTheOtherNA=NaN;
%             ampTheOtherAll(k).ampTheOtherGroup=NaN;
        end
    end
    NAdensity{ii}=curNAdensity;
    FAarea{ii}=curFAarea;
    FAlength{ii}=curFAlength;
    FAdensity{ii}=curFAdensity;
    FAdensityPeri{ii}=curFAdensityPeri;
    FAdensityInside{ii}=curFAdensityInside;
    NAstructGroup{ii}=NAstruct;
    FAstructGroup{ii}=FAstruct;
    if exist('forceAll','var')
        forceAllGroup{ii}=forceAll;
        clear forceAll
    end
    if exist('forceAll','var')
        ampTheOtherAllGroup{ii}=ampTheOtherAll;
        clear ampTheOtherAll
    end
    clear NAstruct FAstruct 
end
%% convert
pixSize = MLAll(1).getMovie(1).pixelSize_; % nm/pix
convertArea = (pixSize/1000)^2;
convertL = pixSize/1000;
%% FA area
try
    FAareaCell=cellfun(@(x) cell2mat(x),FAarea,'Unif',false);
catch
    FAareaCell=cellfun(@(x) cell2mat(x'),FAarea,'Unif',false);
end
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
barPlotCellArray(FAareaCell,nameList,1,1)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (\mum^2)')
title('FA area - bar plot')
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAarea'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAarea'),'-v7.3')

try
    FAareaCellConverted = cellfun(@(x) x, FAareaCell,'unif',false);
    tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
catch
    FAareaCellConverted = cellfun(@(x) x', FAareaCell,'unif',false);
    tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
end
writetable(tableFAarea,strcat(dataPath,'/FAarea.csv'))
%% FA area - boxplot
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
boxPlotCellArray(FAareaCell,nameList,1,false,true)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (\mum^2)')
title('FA area')
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAareaBoxPlot'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAareaBoxPlot'),'-v7.3')
%% FA density
h2=figure; 
barPlotCellArray(FAdensity,nameList)

title('FA density')
ylabel('FA density (#/\mum^2)')
hgexport(h2,strcat(figPath,'/FAdensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensity'),'-v7.3')

tableFAdensity=table(FAdensity,'RowNames',nameList);
writetable(tableFAdensity,strcat(dataPath,'/FAdensity.csv'))
%% FA density at Periphery
h2=figure; 
barPlotCellArray(FAdensityPeri,nameList)

title({'FA density in the cell periphery'; ['(up to ' num2str(bandNA) ' um from the cell edge)']})
ylabel('FA density (#/\mum^2)')
hgexport(h2,strcat(figPath,'/FAdensityPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensityPeri'),'-v7.3')
tableFAdensityPeri=table(FAdensityPeri,'RowNames',nameList);
writetable(tableFAdensityPeri,strcat(dataPath,'/FAdensityPeri.csv'))
%% FA area - cell periphery
faAreaPeri = cellfun(@(x) cell2mat(arrayfun(@(y) y.FAareaPeri,x,'unif',false)'),FAstructGroup','unif',false);
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
boxPlotCellArray(faAreaPeri',nameList,1,false,true)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (\mum^2)')
title({'FA area in the cell periphery'; ['(up to ' num2str(bandNA) ' um from the cell edge)']})
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAareaPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAareaPeri'),'-v7.3')
%% FA area periphery bar plot
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
barPlotCellArray(faAreaPeri',nameList,1,true)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (\mum^2)')
title({'FA area in the cell periphery'; ['(up to ' num2str(bandNA) ' um from the cell edge)']})
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAareaPeriBar'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAareaPeriBar'),'-v7.3')
%% FA length - cell periphery
faLengthPeri = cellfun(@(x) cell2mat(arrayfun(@(y) y.FAlengthPeri,x,'unif',false)'),FAstructGroup','unif',false);
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
boxPlotCellArray(faLengthPeri',nameList,1,false,true)
ylabel('FA length (\mum)')
title({'FA length in the cell periphery'; ['(up to ' num2str(bandNA) ' um from the cell edge)']})
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAlengthPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAlengthPeri'),'-v7.3')
%% FA density Inside
h2=figure; 
barPlotCellArray(FAdensityInside,nameList)

title('FA density inside a cell (from the center to the 5 um from the edge)')
ylabel('FA density (#/\mum^2)')
hgexport(h2,strcat(figPath,'/FAdensityInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensityInside'),'-v7.3')
tableFAdensityInside=table(FAdensityInside,'RowNames',nameList);
writetable(tableFAdensityInside,strcat(dataPath,'/FAdensityInside.csv'))
%% NA density
h3=figure; 
barPlotCellArray(NAdensity,nameList,1)

title('NA density')
ylabel('NA density (#/um^2)')
hgexport(h3,strcat(figPath,'/NAdensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h3,strcat(figPath,'/NAdensity'),'-v7.3')
tableNAdensity=table(NAdensity,'RowNames',nameList);
writetable(tableNAdensity,strcat(dataPath,'/NAdensity.csv'))
%% NA density- boxplot
h3=figure; 
boxPlotCellArray(NAdensity,nameList,1,false,true);

title('NA density')
ylabel('NA density (#/um^2)')
hgexport(h3,strcat(figPath,'/NAdensityBox'),hgexport('factorystyle'),'Format','eps')
hgsave(h3,strcat(figPath,'/NAdensityBox'),'-v7.3')
%% FA length
try
    FAlengthCell=cellfun(@(x) cell2mat(x),FAlength,'Unif',false);
catch
    FAlengthCell=cellfun(@(x) cell2mat(x'),FAlength,'Unif',false);
end
h4=figure; 
barPlotCellArray(FAlengthCell,nameList,1)
title('FA length')
ylabel('FA length (\mum)')
hgexport(h4,strcat(figPath,'/FAlength'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAlength'),'-v7.3')
%% FA length - boxplot
h1=figure; 
boxPlotCellArray(FAlengthCell,nameList,1,false,true)
title('FA length')
ylabel('FA length (um)')
hgexport(h1,strcat(figPath,'/FAlengthBoxPlot'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAlengthBoxPlot'),'-v7.3')

FAlenthCellConverted = cellfun(@(x) x, FAlengthCell,'unif',false);
tableFAlength=table(FAlenthCellConverted,'RowNames',nameList);
writetable(tableFAlength,strcat(dataPath,'/FAlength.csv'))
%% top 10 percentile
percLT=10;
perc=percLT/100;
FAlengthCellSmall=cell(numel(FAlengthCell),1);
FAareaCellSmall=cell(numel(FAareaCell),1);
FAareaPeriSmall=cell(numel(faAreaPeri),1);
for ii=1:numel(FAlengthCell)
    FAlengthCellSmall{ii,1} = ...
    quantile(FAlengthCell{ii},(1-perc)+(perc-0.0001)*rand(1,round((perc-0.0001)*sum(~isnan(FAlengthCell{ii})))));
    FAareaCellSmall{ii,1} = ...
    quantile(FAareaCell{ii},(1-perc)+(perc-0.0001)*rand(1,round((perc-0.0001)*sum(~isnan(FAlengthCell{ii})))));
    FAareaPeriSmall{ii,1} = ...
    quantile(faAreaPeri{ii},(1-perc)+(perc-0.0001)*rand(1,round((perc-0.0001)*sum(~isnan(FAlengthCell{ii})))));
end
%% FA length - boxplot -top 10 percentile
h1=figure;
plotSuccess=boxPlotCellArray(FAlengthCellSmall,nameList,1,false,true);
if plotSuccess
    ylabel(['Focal adhesion length (\mum)'])
    title(['Focal adhesion length (top ' num2str(percLT) ' percentile)'])
    ylim auto
    hgexport(h1,[figPath filesep 'faLengthTop' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faLengthTop' num2str(percLT)],'-v7.3')
    print(h1,[figPath filesep 'faLengthTop' num2str(percLT)],'-dtiff')
    faAreaTop20=table(FAlengthCellSmall,'RowNames',nameList);
    writetable(faAreaTop20,[dataPath filesep 'faLengthTop' num2str(percLT) '.csv'],'WriteRowNames',true)    
end
%% FA area - boxplot -top 10 percentile
h1=figure;
plotSuccess=boxPlotCellArray(FAareaCellSmall,nameList,1,false,true);
if plotSuccess
    ylabel(['Focal adhesion area (\mum^2)'])
    title(['Focal adhesion area (top ' num2str(percLT) ' percentile)'])
    ylim auto
    hgexport(h1,[figPath filesep 'faAreaTop' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faAreaTop' num2str(percLT)],'-v7.3')
    print(h1,[figPath filesep 'faAreaTop' num2str(percLT)],'-dtiff')
    faAreaTop10=table(FAlengthCellSmall,'RowNames',nameList);
    writetable(faAreaTop10,[dataPath filesep 'faAreaTop' num2str(percLT) '.csv'],'WriteRowNames',true)    
end
%% FA area - barplot -top 10 percentile
h1=figure;
barPlotCellArray(FAareaCellSmall,nameList,1,true);
ylabel('Focal adhesion area (\mum^2)')
title(['Focal adhesion area (top ' num2str(percLT) ' percentile)'])
ylim auto
hgexport(h1,[figPath filesep 'faAreaTopBar' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'faAreaTopBar' num2str(percLT)],'-v7.3')
print(h1,[figPath filesep 'faAreaTopBar' num2str(percLT)],'-dtiff')
%% FA area peri- boxplot -top 10 percentile
h1=figure;
boxPlotCellArray(FAareaPeriSmall,nameList,1,false,true);
ylabel(['Focal adhesion area (\mum^2)'])
title(['Focal adhesion area at periphery (top ' num2str(percLT) ' percentile)'])
ylim auto
hgexport(h1,[figPath filesep 'faAreaPeriTopBox' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'faAreaPeriTopBox' num2str(percLT)],'-v7.3')
print(h1,[figPath filesep 'faAreaPeriTopBox' num2str(percLT)],'-dtiff')
faAreaPeriTopBox10=table(FAareaPeriSmall,'RowNames',nameList);
writetable(faAreaPeriTopBox10,[dataPath filesep 'faAreaPeriTopBox' num2str(percLT) '.csv'],'WriteRowNames',true)    
%% FA area peri- barplot -top 10 percentile
h1=figure;
boxPlotCellArray(FAareaPeriSmall,nameList,1,false,true); %Now we don't need to convert again.
ylabel(['Focal adhesion area (\mum^2)'])
title(['Focal adhesion area at periphery (top ' num2str(percLT) ' percentile)'])
ylim auto
hgexport(h1,[figPath filesep 'faAreaPeriTopBar' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'faAreaPeriTopBar' num2str(percLT)],'-v7.3')
print(h1,[figPath filesep 'faAreaPeriTopBar' num2str(percLT)],'-dtiff')
faAreaPeriTopBar10=table(FAareaPeriSmall,'RowNames',nameList);
writetable(faAreaPeriTopBar10,[dataPath filesep 'faAreaPeriTopBar' num2str(percLT) '.csv'],'WriteRowNames',true)    
%% Ratio of FA over NA
FAtoNAratio = cell(numConditions,1);
for ii=1:numConditions
    FAtoNAratio{ii} = [];
    for k=1:N(ii)
        FAtoNAratio{ii} = [FAtoNAratio{ii}; FAstructGroup{ii}(k).numberFA/NAstructGroup{ii}(k).numberNA];
    end
end
%% Plotting Ratio of FA over NA
h4=figure; 
try
    barPlotCellArray(FAtoNAratio,nameList);
catch
    text(0,0.2,'Not plottable')
end
title('Ratio of FA over NA')
ylabel('Ratio of FA over NA (ratio)')
hgexport(h4,strcat(figPath,'/FAtoNARatio'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAtoNARatio'),'-v7.3')

tableFAtoNAratio=table(FAtoNAratio,'RowNames',nameList);
writetable(tableFAtoNAratio,strcat(dataPath,'/FAtoNAratio.csv'))
%% Overal adhesion area per cell area 
FAareaToCellArea = cell(numConditions,1);
for ii=1:numConditions
    FAareaToCellArea{ii} = [];
    for k=1:N(ii)
        FAareaToCellArea{ii} = [FAareaToCellArea{ii}; FAstructGroup{ii}(k).meanFAarea*FAstructGroup{ii}(k).numberFA/FAstructGroup{ii}(k).cellArea];
    end
end
%% Overal adhesion area per cell area plotting
h4=figure; 
barPlotCellArray(FAareaToCellArea,nameList)
title('FA area over cell area')
ylabel('FA area over cell area (\mum^2)')
hgexport(h4,strcat(figPath,'/FAoverCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAoverCell'),'-v7.3')

tableFAareaToCellArea=table(FAareaToCellArea,'RowNames',nameList);
writetable(tableFAareaToCellArea,strcat(dataPath,'/FAareaToCellArea.csv'))
%% Overal cell area 
cellArea = cell(numConditions,1);
for ii=1:numConditions
    cellArea{ii} = [];
    for k=1:N(ii)
        cellArea{ii} = [cellArea{ii}; FAstructGroup{ii}(k).cellArea];
    end
end
h4=figure; 
barPlotCellArray(cellArea,nameList)
title('Cell area')
ylabel('Cell area (\mum^2)')
hgexport(h4,strcat(figPath,'/CellArea'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/CellArea'),'-v7.3')

tableFAareaToCellArea=table(FAareaToCellArea,'RowNames',nameList);
writetable(tableFAareaToCellArea,strcat(dataPath,'/FAareaToCellArea.csv'))
%% FA quantity
% numFA = cellfun(@(x,y) x.*y,FAdensity,cellArea,'unif',false);
numFA = cellfun(@(x) cell2mat(arrayfun(@(y) y.numberFA, x,'unif',false)'),FAstructGroup','unif',false)';
h2=figure; 
barPlotCellArray(numFA,nameList)

title('FA quantity per cell')
ylabel('FA quantity (#/cell)')
hgexport(h2,strcat(figPath,'/FAquantity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAquantity'),'-v7.3')

tableFAquantity=table(numFA,'RowNames',nameList);
writetable(tableFAquantity,strcat(dataPath,'/FAquantity.csv'))
%% NA quantity
% numNA = cellfun(@(x,y) x.*y,NAdensity,cellArea,'unif',false);
numNA = cellfun(@(x) cell2mat(arrayfun(@(y) y.numberNA, x,'unif',false)'),NAstructGroup','unif',false)';
% bandAreaNA = cellfun(@(x) cell2mat(arrayfun(@(y) y.bandArea, x,'unif',false)'),NAstructGroup','unif',false);
h2=figure; 
barPlotCellArray(numNA,nameList)
% boxPlotCellArray(numNA,nameList,1,0,1,1)

% barPlotCellArray(bandAreaNA',nameList)

title('NA quantity per cell')
ylabel('NA quantity (#/cell)')
hgexport(h2,strcat(figPath,'/NAquantity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/NAquantity'),'-v7.3')

tableNAquantity=table(numNA,'RowNames',nameList);
writetable(tableNAquantity,strcat(dataPath,'/NAquantity.csv'))
%% Adhesion eccentricity
eccCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ecc, x,'unif',false)'),FAstructGroup','unif',false);
h4=figure; 
boxPlotCellArray(eccCell,nameList',1,0,1)
title('Adhesion eccentricity')
ylabel('Eccentricity (1)')
hgexport(h4,strcat(figPath,'/ecc'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/ecc'),'-v7.3')

tableEcc=table(eccCell','RowNames',nameList);
writetable(tableEcc,strcat(dataPath,'/ecc.csv'))
%% Ecc - top 10 percentile
percLT=20;
perc=percLT/100;
eccCellSmall=cell(numel(eccCell),1);
for ii=1:numel(eccCell)
    percInd = (1-0.9*perc)+(0.4*perc)*randn(1,round((perc-0.0001)*sum(~isnan(eccCell{ii}))));
    percInd(percInd>1) = 1;
    eccCellSmall{ii,1} = ...
    quantile(eccCell{ii},percInd);
end
h4=figure; 
boxPlotCellArray(eccCellSmall,nameList,1,0,1), ylim auto
title(['Adhesion eccentricity top ' num2str(percLT) ' percentile'])
ylabel('Eccentricity (1)')
hgexport(h4,strcat(figPath,'/eccT', num2str(percLT), 'P'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/eccT', num2str(percLT), 'P'),'-v7.3')

%% Eccentricity barplot
h4=figure; 
barPlotCellArray(eccCell,nameList',1)
title('Adhesion eccentricity')
ylabel('Eccentricity (1)')
hgexport(h4,strcat(figPath,'/eccBar'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/eccBar'))

%% Adhesion width
widthCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.width, x,'unif',false)'),FAstructGroup','unif',false);
h4=figure; 
boxPlotCellArray(widthCell,nameList',1,0,1)
title('Adhesion width')
ylabel('width (\mum)')
hgexport(h4,strcat(figPath,'/FAwidth'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAwidth'),'-v7.3')

tablewidth=table(widthCell','RowNames',nameList);
writetable(tablewidth,strcat(dataPath,'/width.csv'))
%% Adhesion width - bottom 30 percentile
percLT=30;
perc=percLT/100;
widthCellSmall=cell(numel(widthCell),1);
for ii=1:numel(widthCell)
    widthCellSmall{ii,1} = ...
    quantile(widthCell{ii},(perc)+0.2*(perc)*randn(1,round((perc-0.0001)*sum(~isnan(widthCell{ii})))));
end
h4=figure; 
boxPlotCellArray(widthCellSmall,nameList,1,0,1); ylim auto
title(['Adhesion width bottom ' num2str(percLT) ' percentile'])
ylabel('width (\mum)')
hgexport(h4,strcat(figPath,'/FAwidthB',num2str(percLT),'P'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAwidthB',num2str(percLT),'P'),'-v7.3')

%% Width barplot
h5=figure; 
barPlotCellArray(widthCell,nameList',convertL)
title('Adhesion width')
ylabel('width (\mum)')
hgexport(h5,strcat(figPath,'/FAwidthBar'),hgexport('factorystyle'),'Format','eps')
hgsave(h5,strcat(figPath,'/FAwidthBar'),'-v7.3')
%% Force per focal adhesion
try
    forcesFACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceFA, x)),forceAllGroup','unif',false);
catch
    forcesFACell = cellfun(@(x) cell2mat(arrayfun(@(y) cell2mat(y.forceFA), x,'unif',false)),forceAllGroup','unif',false);
end
h4=figure; 
boxPlotCellArray(forcesFACell,nameList',1,0,1)
title('Traction per FA')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceFA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceFA'),'-v7.3')

%% Force per focal complex
try
    forcesFCCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceFC, x)),forceAllGroup','unif',false);
catch
    forcesFCCell = cellfun(@(x) cell2mat(arrayfun(@(y) cell2mat(y.forceFC), x,'unif',false)),forceAllGroup','unif',false);
end
h4=figure; 
boxPlotCellArray(forcesFCCell,nameList',1,0,1)
title('Traction per FC')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceFC'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceFC'),'-v7.3')
%% Force per nascent adhesion
try
    forcesNACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceNA, x)),forceAllGroup','unif',false);
catch
    forcesNACell = cellfun(@(x) cell2mat(arrayfun(@(y) cell2mat(y.forceNA), x, 'unif', false)),forceAllGroup','unif',false);
end
h4=figure; 
boxPlotCellArray(forcesNACell,nameList',1,0,1)
title('Traction per NA')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceNA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceNA'),'-v7.3')
%% the amplitude of the other channel (than TFM and than the adhesion) for FAs
if isfield(FAstructGroup{1}(1),'ampTheOther')
    try
        ampTheOtherFACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherFA, x)),ampTheOtherAllGroup','unif',false);
    catch
        indCell = cellfun(@(x) cell2mat(arrayfun(@(y) iscell(y.ampTheOtherFA), x, 'unif', false)),ampTheOtherAllGroup','unif',false);
        ampTheOtherFACell = cellfun(@(x,z) cell2mat(arrayfun(@(y) cell2mat(y.ampTheOtherFA), x(z), 'unif', false)),ampTheOtherAllGroup',indCell,'unif',false);
    end
    h4=figure; 
    boxPlotCellArray(ampTheOtherFACell,nameList',1,0,1)
    title('The amplitude (background-subtracted) of the third channel per FA')
    ylabel('The amplitude (a.u.)')
    hgexport(h4,strcat(figPath,'/ampTheOtherFA'),hgexport('factorystyle'),'Format','eps')
    hgsave(h4,strcat(figPath,'/ampTheOtherFA'))
end
%% the amplitude of the other channel (than TFM and than the adhesion) for FCs
if isfield(FAstructGroup{1}(1),'ampTheOther')
    try
        ampTheOtherFCCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherFC, x)),ampTheOtherAllGroup','unif',false);
    catch
        indCell = cellfun(@(x) cell2mat(arrayfun(@(y) iscell(y.ampTheOtherFC), x, 'unif', false)),ampTheOtherAllGroup','unif',false);
        ampTheOtherFCCell = cellfun(@(x,z) cell2mat(arrayfun(@(y) cell2mat(y.ampTheOtherFC), x(z), 'unif', false)),ampTheOtherAllGroup',indCell,'unif',false);
    end
    h4=figure; 
    boxPlotCellArray(ampTheOtherFCCell,nameList',1,0,1)
    title('The amplitude (background-subtracted) of the third channel per FC')
    ylabel('The amplitude (a.u.)')
    hgexport(h4,strcat(figPath,'/ampTheOtherFC'),hgexport('factorystyle'),'Format','eps')
    hgsave(h4,strcat(figPath,'/ampTheOtherFC'))
end
%% the amplitude of the other channel (than TFM and than the adhesion) for NAs
if isfield(FAstructGroup{1}(1),'ampTheOther')
    try
        ampTheOtherNACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherNA, x)),ampTheOtherAllGroup','unif',false);
    catch
        indCell = cellfun(@(x) cell2mat(arrayfun(@(y) iscell(y.ampTheOtherNA), x, 'unif', false)),ampTheOtherAllGroup','unif',false);
        ampTheOtherNACell = cellfun(@(x,z) cell2mat(arrayfun(@(y) cell2mat(y.ampTheOtherNA), x(z), 'unif', false)),ampTheOtherAllGroup',indCell,'unif',false);
    end
    h4=figure; 
    boxPlotCellArray(ampTheOtherNACell,nameList',1,0,1)
    title('The amplitude (background-subtracted) of the third channel per NA')
    ylabel('The amplitude (a.u.)')
    hgexport(h4,strcat(figPath,'/ampTheOtherNA'),hgexport('factorystyle'),'Format','eps')
    hgsave(h4,strcat(figPath,'/ampTheOtherNA'))
end
%% Plotting the other channel amplitudes
if isfield(FAstructGroup{1}(1),'ampTheOther')
    ampTheOther = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOther, x','unif',false)'),FAstructGroup','unif',false);
    h4=figure; 
    boxPlotCellArray(ampTheOther,nameList',1,0,1)
    title('Amplitude of the other channel')
    ylabel('Fluorescence amplitude (a.u.)')
    hgexport(h4,strcat(figPath,'/ampTheOther'),hgexport('factorystyle'),'Format','eps')
    hgsave(h4,strcat(figPath,'/ampTheOther'),'-v7.3')

    tableAmpTheOther=table(ampTheOther','RowNames',nameList);
    writetable(tableAmpTheOther,strcat(dataPath,'/ampTheOther.csv'))
end
%% the name for all together
nameListAdh={'NA','FC','FA'};
nameListAdhComb=cell(numel(eccCell)*3,1);

for ii=1:numel(eccCell)
    p=ii-1;
    for jj=1:3
        nameListAdhComb{3*p+jj,1} = [nameList{ii} '-' nameListAdh{jj}];
    end
end

%% the amplitude of the other channel all together
if exist('ampTheOtherNACell','var')
    amplitudeGroupAll=cell(numel(ampTheOtherNACell)*3,1);
    amplitudeGroup={ampTheOtherNACell, ampTheOtherFCCell, ampTheOtherFACell};
    for ii=1:numel(ampTheOtherNACell)
        p=ii-1;
        for jj=1:3
            amplitudeGroupAll{3*p+jj,1} = amplitudeGroup{jj}{ii};
        end
    end
    h1=figure; 
    plotSuccess=boxPlotCellArray(amplitudeGroupAll,nameListAdhComb,1,false,true);
    if plotSuccess
        ylabel('Fluorescence amplitude (a.u.)')
        title('Background-subtracted amplitude of the other channel')
        hgexport(h1,[figPath filesep 'amplitudeTheOtherAll'],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'amplitudeTheOtherAll'])
        print(h1,[figPath filesep 'amplitudeTheOtherAll'],'-dtiff')
    end
end
%% amplitude the other bar plot
if isfield(FAstructGroup{1}(1),'ampTheOther') && ~isempty(amplitudeGroupAll{1})
    h1=figure; 
    barPlotCellArray(amplitudeGroupAll,nameListAdhComb);
    ylabel('Fluorescence amplitude (a.u.)')
    title('Background-subtracted amplitude of the other channel')
    hgexport(h1,[figPath filesep 'amplitudeTheOtherAllBar'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'amplitudeTheOtherAllBar'])
    print(h1,[figPath filesep 'amplitudeTheOtherAllBar'],'-dtiff')
end
%% the force in adhesions all together
forceGroup = {forcesNACell, forcesFCCell,forcesFACell};
for ii=1:numel(forcesNACell)
    p=ii-1;
    for jj=1:3
        forceGroupAll{3*p+jj,1} = forceGroup{jj}{ii};
    end
end
h1=figure; 
boxPlotCellArray(forceGroupAll,nameListAdhComb,1,false,true);
title('Traction at all adhesions')
ylabel('Traction (Pa)')
hgexport(h1,strcat(figPath,'/forceAll'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/forceAll'),'-v7.3')

%% saving
close all
save([dataPath filesep 'adhesionData.mat'],'-v7.3');
