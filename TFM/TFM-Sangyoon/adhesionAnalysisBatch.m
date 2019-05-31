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
    end
    NAdensity{ii}=curNAdensity;
    FAarea{ii}=curFAarea;
    FAlength{ii}=curFAlength;
    FAdensity{ii}=curFAdensity;
    FAdensityPeri{ii}=curFAdensityPeri;
    FAdensityInside{ii}=curFAdensityInside;
    NAstructGroup{ii}=NAstruct;
    FAstructGroup{ii}=FAstruct;
    clear NAstruct FAstruct
end
%% convert
pixSize = MLAll(1).getMovie(1).pixelSize_; % nm/pix
convertArea = (pixSize/1000)^2;
convertL = pixSize/1000;
%% setting up group name
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
nameList=groupNames'; 
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
ylabel('FA area (um^2)')
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAarea'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAarea'),'-v7.3')

FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
writetable(tableFAarea,strcat(dataPath,'/FAarea.csv'))
%% FA area - boxplot
h1=figure; 
% faAreaConverted=cellfun(@(x) x*convertArea,FAarea,'uniformoutput',false);
boxPlotCellArray(FAareaCell,nameList,convertArea,false,true)
% ylim([0 max(mean(FAareaCtrl),mean(FAareaGamma))*convertArea*1.2])
ylabel('FA area (um^2)')
% set(gca,'XTickLabel',{'Control' 'PIP5K-\gamma'})
hgexport(h1,strcat(figPath,'/FAareaBoxPlot'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAareaBoxPlot'),'-v7.3')

%% FA density
h2=figure; 
barPlotCellArray(FAdensity,nameList)

title('FA density')
ylabel('FA density (#/um^2)')
hgexport(h2,strcat(figPath,'/FAdensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensity'),'-v7.3')

tableFAdensity=table(FAdensity,'RowNames',nameList);
writetable(tableFAdensity,strcat(dataPath,'/FAdensity.csv'))
%% FA density at Periphery
h2=figure; 
barPlotCellArray(FAdensityPeri,nameList)

title('FA density in the cell periphery (up to 5 um from the cell edge)')
ylabel('FA density (#/um^2)')
hgexport(h2,strcat(figPath,'/FAdensityInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(figPath,'/FAdensityInside'),'-v7.3')
tableFAdensityPeri=table(FAdensityPeri,'RowNames',nameList);
writetable(tableFAdensityPeri,strcat(dataPath,'/FAdensityPeri.csv'))
%% FA density Inside
h2=figure; 
barPlotCellArray(FAdensityInside,nameList)

title('FA density inside a cell (from the center to the 5 um from the edge)')
ylabel('FA density (#/um^2)')
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
    FAlenthCell=cellfun(@(x) cell2mat(x),FAlength,'Unif',false);
catch
    FAlenthCell=cellfun(@(x) cell2mat(x'),FAlength,'Unif',false);
end
h4=figure; 
barPlotCellArray(FAlenthCell,nameList,convertL)
title('FA length')
ylabel('FA length (um)')
hgexport(h4,strcat(figPath,'/FAlength'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAlength'),'-v7.3')
%% FA length - boxplot
h1=figure; 
boxPlotCellArray(FAlenthCell,nameList,convertL,false,true)
title('FA length')
ylabel('FA length (um)')
hgexport(h1,strcat(figPath,'/FAlengthBoxPlot'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FAlengthBoxPlot'),'-v7.3')

FAlenthCellConverted = cellfun(@(x) x*convertL, FAlenthCell,'unif',false);
tableFAlength=table(FAlenthCellConverted,'RowNames',nameList);
writetable(tableFAlength,strcat(dataPath,'/FAlength.csv'))
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
ylabel('FA area over cell area (ratio)')
hgexport(h4,strcat(figPath,'/FAoverCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h4,strcat(figPath,'/FAoverCell'),'-v7.3')

tableFAareaToCellArea=table(FAareaToCellArea,'RowNames',nameList);
writetable(tableFAareaToCellArea,strcat(dataPath,'/FAareaToCellArea.csv'))
%% Plotting the other channel amplitudes - Mehdi: fill here

%% saving
save([dataPath filesep 'adhesionData.mat'],'-v7.3');
