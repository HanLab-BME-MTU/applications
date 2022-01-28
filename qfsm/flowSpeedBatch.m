% flowSpeedBatch.m is a script that collects speed data analyzed and
% sampled through QFSM and Windowing packages via quantifyMovieFlowSpeed
% function. It works similarly with adhesionAnalysisBatch and
% strainEnergyBatch.
% Sangyoon Han Nov 2021
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
summaryPath = [rootAnalysis '/FlowSpeedSummary' specificName];
ii=0;
while exist(summaryPath, 'dir')
    ii=ii+1;
    summaryPath = [rootAnalysis '/FlowSpeedSummary' specificName num2str(ii)];
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

%% Run quantifyMovieFlowSpeed for each list
N=zeros(numConditions,1);
SpeedL1Group = cell(numConditions,1);
SpeedL2Group = cell(numConditions,1);
SpeedL3Group = cell(numConditions,1);
SpeedL4Group = cell(numConditions,1);
SpeedL5Group = cell(numConditions,1);
SpeedAllGroup = cell(numConditions,1);

sampleMovie = MLAll(1).movies_{1};

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);

    curNAdensity = zeros(N(ii),1);
    SpeedL1 = cell(N(ii),1);
    SpeedL2 = cell(N(ii),1);
    SpeedL3 = cell(N(ii),1);
    SpeedL4 = cell(N(ii),1);
    SpeedL5 = cell(N(ii),1);
    SpeedAll = cell(N(ii),1);
    curML = MLAll(ii);
    
    % Combine data per each condition (1,2,3,4 for 3.4, 18, 100, 100Y, respectively)
    for k=1:N(ii)
        speedCell = quantifyMovieFlowSpeed(curML.movies_{k});
        % output will be in speedCell(windows:layers:frame) format
        % Will need to separately store them
        SpeedL1{k} = squeeze(speedCell(:,1,:));
        SpeedL2{k} = squeeze(speedCell(:,2,:));
        SpeedL3{k} = squeeze(speedCell(:,3,:));
        SpeedL4{k} = squeeze(speedCell(:,4,:));
        SpeedL5{k} = squeeze(speedCell(:,5,:));
        SpeedAll{k} = speedCell;
    end
    SpeedL1Group{ii} = SpeedL1;
    SpeedL2Group{ii} = SpeedL2;
    SpeedL3Group{ii} = SpeedL3;
    SpeedL4Group{ii} = SpeedL4;
    SpeedL5Group{ii} = SpeedL5;
    SpeedAllGroup{ii} = SpeedAll;
    clear SpeedL1 SpeedL2 SpeedL3 SpeedL4 SpeedL5 SpeedAll  
end

%% Plotting SpeedL1
speedL1Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedL1Group,'unif',false);
% speedL1Cell = cellfun(@(x) x(:), speedL1Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL1Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at first layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL1'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL1'),'-v7.3')
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
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL1Cell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at first layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL1Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL1Scatter'),'-v7.3')
%% Plotting SpeedL1 - only top 10 percentile per movie per frame
speedL1CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    speedL1Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(SpeedL1Group{ii}{k},2),1);
        for p=1:size(SpeedL1Group{ii}{k},2)
            curL = length(SpeedL1Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(SpeedL1Group{ii}{k}(:,k),round(curL/10));
        end
        speedL1Top10{k,1} = cell2mat(top10ForFrame);
    end
    speedL1CellTop10{ii} = cell2mat(speedL1Top10);
end
h1=figure; 
boxPlotCellArray(speedL1CellTop10,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at first layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL1Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL1Top10'),'-v7.3')
% FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL1CellTop10,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at first layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL1T10Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL1T10Scatter'),'-v7.3')
%% Plotting SpeedL2
speedL2Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedL2Group,'unif',false);
% speedL2Cell = cellfun(@(x) x(:), speedL2Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL2Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at second layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL2'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL2'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL1Cell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at second layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL2Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL2Scatter'),'-v7.3')
%% Plotting SpeedL2 - only top 10 percentile per movie per frame
speedL2CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    speedL2Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(SpeedL2Group{ii}{k},2),1);
        for p=1:size(SpeedL2Group{ii}{k},2)
            curL = length(SpeedL2Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(SpeedL2Group{ii}{k}(:,k),round(curL/10));
        end
        speedL2Top10{k,1} = cell2mat(top10ForFrame);
    end
    speedL2CellTop10{ii} = cell2mat(speedL2Top10);
end
h1=figure; 
boxPlotCellArray(speedL2CellTop10,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at second layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL2Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL2Top10'),'-v7.3')
%% Plotting SpeedL3
speedL3Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedL3Group,'unif',false);
% speedL3Cell = cellfun(@(x) x(:), speedL3Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL3Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at third layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL3'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL3'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL3Cell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at the third layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL3Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL3Scatter'),'-v7.3')
%% Plotting SpeedL3 - only top 10 percentile per movie per frame
speedL3CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    speedL3Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(SpeedL3Group{ii}{k},2),1);
        for p=1:size(SpeedL3Group{ii}{k},2)
            curL = length(SpeedL3Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(SpeedL3Group{ii}{k}(:,k),round(curL/10));
        end
        speedL3Top10{k,1} = cell2mat(top10ForFrame);
    end
    speedL3CellTop10{ii} = cell2mat(speedL3Top10);
end
h1=figure; 
boxPlotCellArray(speedL3CellTop10,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at third layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL3Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL3Top10'),'-v7.3')
%% Plotting SpeedL4
speedL4Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedL4Group,'unif',false);
% speedL4Cell = cellfun(@(x) x(:), speedL4Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL4Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at fourth layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL4'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL4'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL4Cell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at the fourth layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL4Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL4Scatter'),'-v7.3')
%% Plotting SpeedL4 - only top 10 percentile per movie per frame
speedL4CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    speedL4Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(SpeedL4Group{ii}{k},2),1);
        for p=1:size(SpeedL4Group{ii}{k},2)
            curL = length(SpeedL4Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(SpeedL4Group{ii}{k}(:,k),round(curL/10));
        end
        speedL4Top10{k,1} = cell2mat(top10ForFrame);
    end
    speedL4CellTop10{ii} = cell2mat(speedL4Top10);
end
h1=figure; 
boxPlotCellArray(speedL4CellTop10,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at fourth layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL4Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL4Top10'),'-v7.3')
%% Plotting SpeedL5
speedL5Cell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedL5Group,'unif',false);
% speedL5Cell = cellfun(@(x) x(:), speedL5Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL5Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at fifth layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL5'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL5'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedL5Cell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at the fifth layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL5Scatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL5Scatter'),'-v7.3')
%% Plotting SpeedL5 - only top 10 percentile per movie per frame
speedL5CellTop10 = cell(numConditions,1);
for ii=1:numConditions
    speedL5Top10 = cell(N(ii),1);
    for k=1:N(ii)
        top10ForFrame=cell(size(SpeedL5Group{ii}{k},2),1);
        for p=1:size(SpeedL5Group{ii}{k},2)
            curL = length(SpeedL5Group{ii}{k}(:,k));
            top10ForFrame{p} = maxk(SpeedL5Group{ii}{k}(:,k),round(curL/10));
        end
        speedL5Top10{k,1} = cell2mat(top10ForFrame);
    end
    speedL5CellTop10{ii} = cell2mat(speedL5Top10);
end
h1=figure; 
boxPlotCellArray(speedL5CellTop10,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at fifth layer, top 10 percentile')
hgexport(h1,strcat(figPath,'/FlowSpeedL5Top10'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL5Top10'),'-v7.3')
%% Plotting SpeedAll
speedAllCell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), SpeedAllGroup,'unif',false);
% speedAllCell = cellfun(@(x) x(:), speedAllCell,'unif',false);
h1=figure; 
boxPlotCellArray(speedAllCell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at all layers')
hgexport(h1,strcat(figPath,'/FlowSpeedAll'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedAll'),'-v7.3')
%% scatter plot
h1=figure; 
errorBarPlotCellArray(speedAllCell,xValues,1);
xlabel(xLabel)
ylabel('Flow speed (nm/min)')
title('Flow speed at all layers')
hgexport(h1,strcat(figPath,'/FlowSpeedAllScatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedAllScatter'),'-v7.3')
%% saving
close all
save([dataPath filesep 'flowDataAll.mat'],'-v7.3');
disp(['Figures are stored in ' figPath '.'])
disp(['Data are stored in ' dataPath '.'])