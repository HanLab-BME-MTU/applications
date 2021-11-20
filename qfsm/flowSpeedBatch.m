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
    specificName = strjoin(nameList);
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
        speedCell = quantifyMovieFlowSpeed(MD);
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
speedL1Cell = cellfun(@(x) cell2mat(x), SpeedL1Group,'unif',false);
speedL1Cell = cellfun(@(x) x(:), speedL1Cell,'unif',false);
h1=figure; 
boxPlotCellArray(speedL1Cell,nameList,1,1,1)
ylabel('Flow speed (nm/min)')
title('Flow speed at first layer')
hgexport(h1,strcat(figPath,'/FlowSpeedL1'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/FlowSpeedL1'),'-v7.3')
% FAareaCellConverted = cellfun(@(x) x*convertArea, FAareaCell,'unif',false);
% tableFAarea=table(FAareaCellConverted,'RowNames',nameList);
% writetable(tableFAarea,strcat(dataPath,'/FAarea.csv'))
%% Plotting SpeedL1 - only top 10 percentile per movie per frame
speedL1CellTop10 = cell(1,numConditions);
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
