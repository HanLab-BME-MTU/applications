% TractionBatch.m is a script that collects Traction data analyzed and
% sampled through TFM and Windowing packages. Is a rewrite of
% flowSpeedBatch for use with traction instead of speed.
% Sangyoon Han Nov 2021
% Etienne Michels Jul 2022
%% open necessary MLs
clear

[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,~,MLdirect]=chooseSelectedFolders;
% Asking user
disp('The current names are: ')
nameList = groupNames';
disp(nameList)
namesOK = input('Do you want to rename your condition names? (y/n)','s');
if strcmp(namesOK, 'y')
    for ii=1:numel(nameList)
        curName = input(['For ' nameList{ii} ': '], 's');
        if ~isempty(curName)
            nameList{ii} = curName;
        end
    end
    specificName = strjoin(nameList, '_');
end
%% Just in case when there is a larger condition.
largerCondition = input('Do you want to add a larger condition name, e.g., WT or Blebbi etc (y/n)?','s');
if strcmp(largerCondition, 'y')
    largerConditionName = input(['Enter the condition name:'], 's');
    specificName = [largerConditionName specificName];
end

%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
% rootAnalysis = pathAnalysisAll{1};
summaryPath = [rootAnalysis '/ForceFieldSummary' specificName];
ii=0;
while exist(summaryPath, 'dir')
    ii=ii+1;
    summaryPath = [rootAnalysis '/ForceFieldSummary' specificName num2str(ii)];
end
figPath = [summaryPath '/Figs'];
mkdir(figPath)
dataPath = [summaryPath '/Data'];
mkdir(dataPath)
save([rootAnalysis filesep 'selectedFoldersForFlow_' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end

%% Run quantifyMovieFlowTraction for each list
N=zeros(numConditions,1);
TractionL1Group = cell(numConditions,1);
TractionL2Group = cell(numConditions,1);
TractionL3Group = cell(numConditions,1);
TractionL4Group = cell(numConditions,1);
TractionL5Group = cell(numConditions,1);
TractionAllGroup = cell(numConditions,1);


sampleMovie = MLAll(1).movies_{1};

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);

    curNAdensity = zeros(N(ii),1);
    TractionL1 = cell(N(ii),1);
    TractionL2 = cell(N(ii),1);
    TractionL3 = cell(N(ii),1);
    TractionL4 = cell(N(ii),1);
    TractionL5 = cell(N(ii),1);
    TractionAll = cell(N(ii),1);
    curML = MLAll(ii);
    
    % Combine data per each condition (1,2,3,4 for 3.4, 18, 100, 100Y, respectively)
    for k=1:N(ii)
        tMapUnshifted = quantifyTractionWindows(curML.movies_{k});
        % output will be in TractionCell(windows:layers:frame) format
        % Will need to separately store them
        TractionL1{k} = squeeze(tMapUnshifted(:,1,:));
        TractionL2{k} = squeeze(tMapUnshifted(:,2,:));
        TractionL3{k} = squeeze(tMapUnshifted(:,3,:));
        TractionL4{k} = squeeze(tMapUnshifted(:,4,:));
        TractionL5{k} = squeeze(tMapUnshifted(:,5,:));
        TractionAll{k} = tMapUnshifted;

    end
    TractionL1Group{ii} = TractionL1;
    TractionL2Group{ii} = TractionL2;
    TractionL3Group{ii} = TractionL3;
    TractionL4Group{ii} = TractionL4;
    TractionL5Group{ii} = TractionL5;
    TractionAllGroup{ii} = TractionAll;
    
    clear TractionL1 TractionL2 TractionL3 TractionL4 TractionL5 TractionAll  
end

%% Plotting TractionL1
TractionL1Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionL1Group,'unif',false);
% TractionL1Cell = cellfun(@(x) x(:), TractionL1Cell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionL1Cell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at first layer')
hgexport(h1,strcat(figPath,'/TractionL1'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL1'),'-v7.3')
% FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
% tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
% writetable(tableFAarea,strcat(dataPath,'/FAarea.csv'))
%% Scatter plot
% name should be only numeric for scatter plot
try
    xValues = cellfun(@(x) str2double(x), nameList);
    xLabel = input(['Label and unit? e.g. Stiffness (kPa): '], 's');
catch
    disp('Your initial x labels were not numerical texts. Run this code again with renaming the labels.')
end
%% scatter plot - TractionL1
h1=figure; 
errorBarPlotCellArray(TractionL1Cell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at first layer')
hgexport(h1,strcat(figPath,'/TractionL1Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL1Scatter'),'-v7.3')
%% Plotting TractionL1 - only top 10 percentile per movie per frame
TractionL1CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    TractionL1Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(TractionL1Group{ii}{k},2),1);
        for p=1:size(TractionL1Group{ii}{k},2)
            curL = length(TractionL1Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(TractionL1Group{ii}{k}(:,k),round(curL/10));
        end
        TractionL1Top10{k,1} = cell2mat(top10ForFrame);
    end
    TractionL1CellTop10{ii} = cell2mat(TractionL1Top10);
end
h1=figure; 
boxPlotCellArray(TractionL1CellTop10,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at first layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL1Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL1Top10'),'-v7.3')
% FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
%% scatter plot - TractionL1Top10
h1=figure; 
errorBarPlotCellArray(TractionL1CellTop10,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at first layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL1T10Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL1T10Scatter'),'-v7.3')
%% Plotting TractionL2
TractionL2Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionL2Group,'unif',false);
% TractionL2Cell = cellfun(@(x) x(:), TractionL2Cell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionL2Cell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at second layer')
hgexport(h1,strcat(figPath,'/TractionL2'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL2'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(TractionL2Cell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at second layer')
hgexport(h1,strcat(figPath,'/TractionL2Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL2Scatter'),'-v7.3')
%% Plotting TractionL2 - only top 10 percentile per movie per frame
TractionL2CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    TractionL2Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(TractionL2Group{ii}{k},2),1);
        for p=1:size(TractionL2Group{ii}{k},2)
            curL = length(TractionL2Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(TractionL2Group{ii}{k}(:,k),round(curL/10));
        end
        TractionL2Top10{k,1} = cell2mat(top10ForFrame);
    end
    TractionL2CellTop10{ii} = cell2mat(TractionL2Top10);
end
h1=figure; 
boxPlotCellArray(TractionL2CellTop10,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at second layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL2Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL2Top10'),'-v7.3')
%% Plotting TractionL3
TractionL3Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionL3Group,'unif',false);
% TractionL3Cell = cellfun(@(x) x(:), TractionL3Cell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionL3Cell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at third layer')
hgexport(h1,strcat(figPath,'/TractionL3'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL3'),'-v7.3')
%% scatter plot - TractionL3
h1=figure; 
errorBarPlotCellArray(TractionL3Cell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at the third layer')
hgexport(h1,strcat(figPath,'/TractionL3Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL3Scatter'),'-v7.3')
%% Plotting TractionL3 - only top 10 percentile per movie per frame
TractionL3CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    TractionL3Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(TractionL3Group{ii}{k},2),1);
        for p=1:size(TractionL3Group{ii}{k},2)
            curL = length(TractionL3Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(TractionL3Group{ii}{k}(:,k),round(curL/10));
        end
        TractionL3Top10{k,1} = cell2mat(top10ForFrame);
    end
    TractionL3CellTop10{ii} = cell2mat(TractionL3Top10);
end
h1=figure; 
boxPlotCellArray(TractionL3CellTop10,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at third layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL3Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL3Top10'),'-v7.3')
%% Plotting TractionL4
TractionL4Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionL4Group,'unif',false);
% TractionL4Cell = cellfun(@(x) x(:), TractionL4Cell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionL4Cell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at fourth layer')
hgexport(h1,strcat(figPath,'/TractionL4'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL4'),'-v7.3')
%% scatter plot - TractionL4Cell
h1=figure; 
errorBarPlotCellArray(TractionL4Cell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at the fourth layer')
hgexport(h1,strcat(figPath,'/TractionL4Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL4Scatter'),'-v7.3')
%% Plotting TractionL4 - only top 10 percentile per movie per frame
TractionL4CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    TractionL4Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(TractionL4Group{ii}{k},2),1);
        for p=1:size(TractionL4Group{ii}{k},2)
            curL = length(TractionL4Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(TractionL4Group{ii}{k}(:,k),round(curL/10));
        end
        TractionL4Top10{k,1} = cell2mat(top10ForFrame);
    end
    TractionL4CellTop10{ii} = cell2mat(TractionL4Top10);
end
h1=figure; 
boxPlotCellArray(TractionL4CellTop10,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at fourth layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL4Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL4Top10'),'-v7.3')
%% Plotting TractionL5
TractionL5Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionL5Group,'unif',false);
% TractionL5Cell = cellfun(@(x) x(:), TractionL5Cell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionL5Cell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at fifth layer')
hgexport(h1,strcat(figPath,'/TractionL5'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL5'),'-v7.3')
%% scatter plot - TractionL5Cell
h1=figure; 
errorBarPlotCellArray(TractionL5Cell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at the fifth layer')
hgexport(h1,strcat(figPath,'/TractionL5Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL5Scatter'),'-v7.3')
%% Plotting TractionL5 - only top 10 percentile per movie per frame
TractionL5CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    TractionL5Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(TractionL5Group{ii}{k},2),1);
        for p=1:size(TractionL5Group{ii}{k},2)
            curL = length(TractionL5Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(TractionL5Group{ii}{k}(:,k),round(curL/10));
        end
        TractionL5Top10{k,1} = cell2mat(top10ForFrame);
    end
    TractionL5CellTop10{ii} = cell2mat(TractionL5Top10);
end
h1=figure; 
boxPlotCellArray(TractionL5CellTop10,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at fifth layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/TractionL5Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionL5Top10'),'-v7.3')
%% Plotting TractionAll
TractionAllCell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), TractionAllGroup,'unif',false);
% TractionAllCell = cellfun(@(x) x(:), TractionAllCell,'unif',false);
h1=figure; 
boxPlotCellArray(TractionAllCell,nameList,1,1,1)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at all layers')
hgexport(h1,strcat(figPath,'/TractionAll'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionAll'),'-v7.3')
%% scatter plot - TractionAllCell
h1=figure; 
errorBarPlotCellArray(TractionAllCell,xValues,1);
xlabel(xLabel)
ylabel('Traction Map Unshifted')
title('Traction Map Unshifted at all layers')
hgexport(h1,strcat(figPath,'/TractionAllScatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/TractionAllScatter'),'-v7.3')
%% saving
close all
save([dataPath filesep 'UnshiftedTractionMapDataAll.mat'],'-v7.3');
disp(['Figures are stored in ' figPath '.'])
disp(['Data are stored in ' dataPath '.'])