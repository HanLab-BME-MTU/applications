%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
if ~ischar(pathSFolders) && pathSFolders==0
    analysisFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~analysisFolderSelectionDone
        ii=ii+1;
        curPathProject = uigetdir(rootFolder,'Select each analysis folder that contains movieList.mat (Click Cancel when no more)');
        if ~ischar(curPathProject) && curPathProject==0
            analysisFolderSelectionDone=true;
        else
            pathAnalysisAll{ii} = curPathProject;
        end
    end
else
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
end

%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummaryAdhesion/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummaryAdhesion/Data'];
mkdir(dataPath)
%% setting up group name
for ii=1:numConditions
    [dumPath, finalFolder]=fileparts(pathAnalysisAll{ii});
    if isempty(finalFolder)
        [dumPath, finalFolder]=fileparts(dumPath);
        groupNames{ii} = finalFolder;
    else        
        groupNames{ii} = finalFolder;
    end
end
nameList=groupNames'; 

% Asking user
disp('Do you want to rename your condition names? The current names are: ')
disp(nameList)
for ii=1:numel(nameList)
    curName = input(['For ' nameList{ii} ': '], 's');
    if ~isempty(curName)
        nameList{ii} = curName;
    end
end
specificName = strjoin(nameList);
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

iPax = input('Enter adhesion channel number (1 or 2 ..): ');
for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);
    bandNA = []; % 2 um from the edge

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
            FAstruct(k) = curFAstruct(end);
        else %numel(curFAstruct)==1
            curNAdensity(k) = curNAstruct.NAdensity;
            curFAarea{k} = curFAstruct.area;
            curFAlength{k} = curFAstruct.length;
            curFAdensity(k) = curFAstruct.FAdensity;
            NAstruct(k) = curNAstruct;
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
            dataPath= [outputFilePath filesep 'data'];
            forceAll(k) = load([dataPath filesep 'forcePerAdhesionType.mat'],'forceFA', 'forceFC', 'forceNA', 'forceBGinCell', 'forceBGoutCell');
        else
            forceAll(k)=[];
        end
        % Reading the other channels's amplitudes
        nChan = numel(curML.movies_{k}.channels_);
        if nChan>2
            outputFilePath = [curML.movies_{k}.outputDirectory_ filesep 'Adhesion Quantification'];
            dataPath= [outputFilePath filesep 'data'];
            ampTheOtherAll(k) = load([dataPath filesep 'ampTheOtherPerAdhesionType.mat']);
        else
            ampTheOtherAll(k)=[];
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
    forceAllGroup{ii}=forceAll;
    ampTheOtherAllGroup{ii}=ampTheOtherAll;
    clear NAstruct FAstruct ampTheOtherAll forceAll
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
barPlotCellArray(FAareaCell,nameList,convertArea)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (\mum^2)')
title('FA area - bar plot')
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAarea'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAarea'),'-v7.3')

try
    FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
    tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
catch
    FAareaCellConverted = cellfun(@(x) x'*convertArea, FAareaCell,'unif',false);
    tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
end
writetable(tableFAarea,strcat(dataPath,'/FAarea.csv'))
%% FA area - boxplot
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
boxPlotCellArray(FAareaCell,nameList,convertArea,false,true)
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

title('FA density in the cell periphery (up to 5 um from the cell edge)')
ylabel('FA density (#/\mum^2)')
hgexport(h2,strcat(figPath,'/FAdensityInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensityInside'),'-v7.3')
tableFAdensityPeri=table(FAdensityPeri,'RowNames',nameList);
writetable(tableFAdensityPeri,strcat(dataPath,'/FAdensityPeri.csv'))
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
barPlotCellArray(NAdensity,nameList)

title('NA density')
ylabel('NA density (#/um^2)')
hgexport(h3,strcat(figPath,'/NAdensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h3,strcat(figPath,'/NAdensity'),'-v7.3')
tableNAdensity=table(NAdensity,'RowNames',nameList);
writetable(tableNAdensity,strcat(dataPath,'/NAdensity.csv'))
%% FA length
try
    FAlengthCell=cellfun(@(x) cell2mat(x),FAlength,'Unif',false);
catch
    FAlengthCell=cellfun(@(x) cell2mat(x'),FAlength,'Unif',false);
end
h4=figure; 
barPlotCellArray(FAlengthCell,nameList,convertL)
title('FA length')
ylabel('FA length (\mum)')
hgexport(h4,strcat(figPath,'/FAlength'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAlength'),'-v7.3')
%% FA length - boxplot
h1=figure; 
boxPlotCellArray(FAlengthCell,nameList,convertL,false,true)
title('FA length')
ylabel('FA length (um)')
hgexport(h1,strcat(figPath,'/FAlengthBoxPlot'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAlengthBoxPlot'),'-v7.3')

FAlenthCellConverted = cellfun(@(x) x*convertL, FAlengthCell,'unif',false);
tableFAlength=table(FAlenthCellConverted,'RowNames',nameList);
writetable(tableFAlength,strcat(dataPath,'/FAlength.csv'))
%% FA length - boxplot -top 10 percentile
percLT=10;
perc=percLT/100;
FAlengthCellSmall=cell(numel(FAlengthCell),1);
for ii=1:numel(FAlengthCell)
    FAlengthCellSmall{ii,1} = ...
    quantile(FAlengthCell{ii},(1-perc)+(perc-0.01)*rand(1,round((perc-0.01)*sum(~isnan(FAlengthCell{ii})))));
end

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
    barPlotCellArray(FAtoNAratio,nameList)
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
        FAareaToCellArea{ii} = [FAareaToCellArea{ii}; FAstructGroup{ii}(k).meanFAarea*FAstructGroup{ii}(k).numberFA*convertArea/FAstructGroup{ii}(k).cellArea];
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
ylabel('Cell area (ratio)')
hgexport(h4,strcat(figPath,'/CellArea'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/CellArea'),'-v7.3')

tableFAareaToCellArea=table(FAareaToCellArea,'RowNames',nameList);
writetable(tableFAareaToCellArea,strcat(dataPath,'/FAareaToCellArea.csv'))
%% FA quantity
numFA = cellfun(@(x,y) x.*y,FAdensity,cellArea,'unif',false);
h2=figure; 
barPlotCellArray(numFA,nameList)

title('FA quantity per cell')
ylabel('FA quantity (#/cell)')
hgexport(h2,strcat(figPath,'/FAquantity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAquantity'),'-v7.3')

tableFAquantity=table(numFA,'RowNames',nameList);
writetable(tableFAquantity,strcat(dataPath,'/FAquantity.csv'))
%% NA quantity
numNA = cellfun(@(x,y) x.*y,NAdensity,cellArea,'unif',false);
h2=figure; 
barPlotCellArray(numNA,nameList)

title('NA quantity per cell')
ylabel('NA quantity (#/cell)')
hgexport(h2,strcat(figPath,'/NAquantity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/NAquantity'),'-v7.3')

tableNAquantity=table(numNA,'RowNames',nameList);
writetable(tableNAquantity,strcat(dataPath,'/NAquantity.csv'))
%% Plotting the other channel amplitudes
ampTheOther = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOther, x','unif',false)'),FAstructGroup','unif',false);
h4=figure; 
boxPlotCellArray(ampTheOther,nameList',1,0,1)
title('Amplitude of the other channel')
ylabel('Fluorescence amplitude (a.u.)')
hgexport(h4,strcat(figPath,'/ampTheOther'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/ampTheOther'),'-v7.3')

tableAmpTheOther=table(ampTheOther','RowNames',nameList);
writetable(tableAmpTheOther,strcat(dataPath,'/ampTheOther.csv'))
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
boxPlotCellArray(widthCell,nameList',convertL,0,1)
title('Adhesion width')
ylabel('width (\mum)')
hgexport(h4,strcat(figPath,'/FAwidth'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAwidth'),'-v7.3')

tablewidth=table(widthCell','RowNames',nameList);
writetable(tablewidth,strcat(dataPath,'/width.csv'))
%% Width barplot
h5=figure; 
barPlotCellArray(widthCell,nameList',convertL)
title('Adhesion width')
ylabel('width (\mum)')
hgexport(h5,strcat(figPath,'/FAwidthBar'),hgexport('factorystyle'),'Format','eps')
hgsave(h5,strcat(figPath,'/FAwidthBar'),'-v7.3')
%% Force per focal adhesion
forcesFACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceFA, x)),forceAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(forcesFACell,nameList',1,0,1)
title('Traction per FA')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceFA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceFA'),'-v7.3')

%% Force per focal complex
forcesFCCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceFC, x)),forceAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(forcesFCCell,nameList',1,0,1)
title('Traction per FC')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceFC'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceFC'),'-v7.3')
%% Force per nascent adhesion
forcesNACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.forceNA, x)),forceAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(forcesNACell,nameList',1,0,1)
title('Traction per NA')
ylabel('Traction (Pa)')
hgexport(h4,strcat(figPath,'/forceNA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceNA'),'-v7.3')
%% the amplitude of the other channel (than TFM and than the adhesion) for FAs
ampTheOtherFACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherFA, x)),ampTheOtherAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(ampTheOtherFACell,nameList',1,0,1)
title('The amplitude (background-subtracted) of the third channel per FA')
ylabel('The amplitude (a.u.)')
hgexport(h4,strcat(figPath,'/ampTheOtherFA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/ampTheOtherFA'))
%% the amplitude of the other channel (than TFM and than the adhesion) for FCs
ampTheOtherFCCell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherFC, x)),ampTheOtherAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(ampTheOtherFCCell,nameList',1,0,1)
title('The amplitude (background-subtracted) of the third channel per FC')
ylabel('The amplitude (a.u.)')
hgexport(h4,strcat(figPath,'/ampTheOtherFC'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/ampTheOtherFC'))
%% the amplitude of the other channel (than TFM and than the adhesion) for NAs
ampTheOtherNACell = cellfun(@(x) cell2mat(arrayfun(@(y) y.ampTheOtherNA, x)),ampTheOtherAllGroup','unif',false);
h4=figure; 
boxPlotCellArray(ampTheOtherNACell,nameList',1,0,1)
title('The amplitude (background-subtracted) of the third channel per NA')
ylabel('The amplitude (a.u.)')
hgexport(h4,strcat(figPath,'/ampTheOtherNA'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/ampTheOtherNA'))
%% the amplitude of the other channel all together
nameListAdh={'NA','FC','FA'};
nameListAdhComb=cell(numel(ampTheOtherNACell)*3,1);
amplitudeGroupAll=cell(numel(ampTheOtherNACell)*3,1);
amplitudeGroup={ampTheOtherNACell, ampTheOtherFCCell, ampTheOtherFACell};
for ii=1:numel(ampTheOtherNACell)
    p=ii-1;
    for jj=1:3
        nameListAdhComb{3*p+jj,1} = [nameList{ii} '-' nameListAdh{jj}];
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
%% amplitude the other bar plot
h1=figure; 
barPlotCellArray(amplitudeGroupAll,nameListAdhComb);
ylabel('Fluorescence amplitude (a.u.)')
title('Background-subtracted amplitude of the other channel')
hgexport(h1,[figPath filesep 'amplitudeTheOtherAllBar'],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'amplitudeTheOtherAllBar'])
print(h1,[figPath filesep 'amplitudeTheOtherAllBar'],'-dtiff')
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
hgexport(h4,strcat(figPath,'/forceAll'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/forceAll'),'-v7.3')

%% saving
save([dataPath filesep 'adhesionData.mat'],'-v7.3');
