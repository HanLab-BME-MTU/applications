function [] = protrusivityBatch
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
    % saving pathAnalysisAll...
    pathSFolders = pathAnalysisAll{1};
    disp(['Saving these as selectedFolders.mat in ' pathSFolders '...'])
    save([pathSFolders filesep 'selectedFolders.mat'],'pathAnalysisAll');
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
rootAnalysis = (pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummaryProtrusion/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummaryProtrusion/Data'];
mkdir(dataPath)
%% setting up group name
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
nameList=groupNames; 
                
%% Running
protVelCondition = cell(numConditions,1);
retrVelCondition = cell(numConditions,1);
for ii=1:numConditions
    curML=MLAll(ii);
    curNumMDs(ii) = numel(curML.movieDataFile_);
    curMovies = curML.movies_;
    for k=1:curNumMDs(ii)
        % get the tracksNA
        curMovie=curMovies{k};
        % get the strain energy
        iCurProtSamp = curMovie.getProcessIndex('ProtrusionSamplingProcess');
        curProtSampProc = curMovie.getProcess(iCurProtSamp);
        % Load sampling
        curProtSamp=curProtSampProc.loadChannelOutput;
        curNFrames = curMovie.nFrames_;
        % Go through persistent time analysis per time series
        numWin = size(curProtSamp.avgNormal,1);
        deltaT = curMovie.timeInterval_;
        iP=0; iR=0; clear overallProt overallTimeProt overallRetr overallTimeRetr
        motionStateWin = zeros(numWin,1); % 1: pro, -1:ret, 0: quiescent
        for jj=1:numWin
            [prot,retr,motionTS]=getPersistenceTime(curProtSamp.avgNormal(jj,:),deltaT);
            if ~isnan(prot.Veloc(1))
                advance = prot.persTime'*prot.Veloc-retr.persTime'*retr.Veloc;
                posQuiescentMovement= curNFrames*prot.limit*deltaT/10;
                if advance>posQuiescentMovement
                    iP=iP+1;
                    overallProt(jj)=advance;
                    overallTimeProt(jj)=(sum(prot.persTime)+sum(retr.persTime));
                    avgTimeProt(jj)=mean(prot.persTime);
                    overallRetr(jj)=NaN;
                    overallTimeRetr(jj)=NaN;
                    avgTimeRetr(jj)=NaN;
                    motionStateWin(jj)=1;
                elseif advance<-posQuiescentMovement
                    iR=iR+1;
                    overallRetr(jj)=advance;
                    overallTimeRetr(jj)=(sum(prot.persTime)+sum(retr.persTime));
                    avgTimeRetr(jj)=mean(retr.persTime);
                    overallProt(jj)=NaN;
                    overallTimeProt(jj)=NaN;
                    avgTimeProt(jj)=NaN;
                    motionStateWin(jj)=-1;
                end
            end
        end
        fractionProt(k)=sum(motionStateWin==1)/length(motionStateWin);
        fractionQuiet(k)=sum(motionStateWin==0)/length(motionStateWin);
        fractionRet(k)=sum(motionStateWin==-1)/length(motionStateWin);
        speedProtCell{k}=(overallProt)./(overallTimeProt);
        speedRetrCell{k}=(overallRetr)./(overallTimeRetr);
        distProtruded{k}=overallProt;
        distRetracted{k}=overallRetr;
        persTimeProt{k}=avgTimeProt;
        persTimeRetr{k}=avgTimeRetr;
    end
%     protVelCondition{ii} = cell2mat(protVelGroup);
%     retrVelCondition{ii} = cell2mat(retrVelGroup);
    protSpeedGroup{ii}=speedProtCell;
    retrSpeedGroup{ii}=speedRetrCell;
    protDistGroup{ii}=distProtruded;
    retrDistGroup{ii}=distRetracted;
    persTimeProtGroup{ii}=persTimeProt;
    persTimeRetrGroup{ii}=persTimeRetr;
    fractionProtGroup{ii}=fractionProt;
    fractionQuietGroup{ii}=fractionQuiet;
    fractionRetGroup{ii}=fractionRet;
end
%% Ploting - fractionProtGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
% top20pProt = cellfun(@(x) cellfun(@(y) quantile(y,.75:0.001:.99),x, 'unif',false),protSpeedGroup,'unif',false);
% fractionProtGroupArray = cellfun(@(y) cell2mat(y),fractionProtGroup,'unif',false);
boxPlotCellArray(fractionProtGroup,nameList,1,false,true,true)
% barPlotCellArray(fractionProtGroup,nameList,1)
ylabel('fractionProtGroupArray')
title('fractionProtGroupArray')
hgexport(h1,strcat(figPath,'/fractionProtGroup'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/fractionProtGroup'),'-v7.3')
print(h1,strcat(figPath,'/fractionProtGroup.tif'),'-dtiff')

tableFractionProt=table(fractionProtGroup','RowNames',nameList');
writetable(tableFractionProt,strcat(dataPath,'/tableFractionProt.csv'),'WriteRowNames',true)
%% Ploting - fractionRetGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
% top20pProt = cellfun(@(x) cellfun(@(y) quantile(y,.75:0.001:.99),x, 'unif',false),protSpeedGroup,'unif',false);
% fractionProtGroupArray = cellfun(@(y) cell2mat(y),fractionProtGroup,'unif',false);
boxPlotCellArray(fractionRetGroup,nameList,1,false,true)
% barPlotCellArray(fractionProtGroup,nameList,1)
ylabel('fractionRetGroup')
title('fractionRetGroup')
hgexport(h1,strcat(figPath,'/fractionRetGroup'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/fractionRetGroup'),'-v7.3')
print(h1,strcat(figPath,'/fractionRetGroup.tif'),'-dtiff')

tableFractionRet=table(fractionRetGroup','RowNames',nameList');
writetable(tableFractionRet,strcat(dataPath,'/fractionRetGroup.csv'),'WriteRowNames',true)
%% Ploting - fractionQuietGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
% top20pProt = cellfun(@(x) cellfun(@(y) quantile(y,.75:0.001:.99),x, 'unif',false),protSpeedGroup,'unif',false);
% fractionProtGroupArray = cellfun(@(y) cell2mat(y),fractionProtGroup,'unif',false);
boxPlotCellArray(fractionQuietGroup,nameList,1,false,true,true)
% barPlotCellArray(fractionProtGroup,nameList,1)
ylabel('fractionQuietGroup')
title('fractionQuietGroup')
hgexport(h1,strcat(figPath,'/fractionQuietGroup'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/fractionQuietGroup'),'-v7.3')
print(h1,strcat(figPath,'/fractionQuietGroup.tif'),'-dtiff')

tableFractionQuiet=table(fractionQuietGroup','RowNames',nameList');
writetable(tableFractionQuiet,strcat(dataPath,'/fractionQuietGroup.csv'),'WriteRowNames',true)
%% Ploting - protSpeedGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
top20pProt = cellfun(@(x) cellfun(@(y) quantile(y,.25:0.001:.75),x, 'unif',false),protSpeedGroup,'unif',false);
protSpeedGroupArray = cellfun(@(y) cell2mat(y),top20pProt,'unif',false);
boxPlotCellArray(protSpeedGroupArray,nameList,1,false,true)
% barPlotCellArray(protSpeedGroupArray,nameList,1)
ylabel('Protrusion velocity (um/min)')
title('Protrusion velocity')
hgexport(h1,strcat(figPath,'/protVelCondition'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/protVelCondition'),'-v7.3')
print(h1,strcat(figPath,'/protVelCondition.tif'),'-dtiff')

tableProtSpeed=table(protSpeedGroupArray','RowNames',nameList');
writetable(tableProtSpeed,strcat(dataPath,'/protSpeedGroup.csv'),'WriteRowNames',true)
%% Ploting - retrSpeedGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
top20pRetr = cellfun(@(x) cellfun(@(y) quantile(y,.05:0.001:.75),x, 'unif',false),retrSpeedGroup,'unif',false);
retrSpeedGroupArray = cellfun(@(y) cell2mat(y),top20pRetr,'unif',false);
boxPlotCellArray(retrSpeedGroupArray,nameList,1,false,true)
% barPlotCellArray(retrSpeedGroupArray,nameList,1)
ylabel('Retraction velocity (um/min)')
title('Retraction velocity')
hgexport(h1,strcat(figPath,'/retrVelCondition'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/retrVelCondition'),'-v7.3')
print(h1,strcat(figPath,'/retrVelCondition.tif'),'-dtiff')

tableRetrSpeed=table(retrSpeedGroupArray','RowNames',nameList');
writetable(tableRetrSpeed,strcat(dataPath,'/retrSpeedGroup.csv'),'WriteRowNames',true)
%% Ploting - persTimeProtGroup
h1=figure; 
% barPlotCellArray(protVelCondition,nameList',1)
top20pProtTime = cellfun(@(x) cellfun(@(y) quantile(y,.50:0.001:.99),x, 'unif',false),persTimeProtGroup,'unif',false);
protTimeGroupArray = cellfun(@(y) cell2mat(y),top20pProtTime,'unif',false);
boxPlotCellArray(protTimeGroupArray,nameList,1,false,true)
% barPlotCellArray(protTimeGroupArray,nameList,1/60)
ylabel('Protrusion persistent time (min)')
title('Protrusion persistent time')
hgexport(h1,strcat(figPath,'/protTime'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/protTime'),'-v7.3')
print(h1,strcat(figPath,'/protTime.tif'),'-dtiff')

tableProtTime=table(protTimeGroupArray','RowNames',nameList');
writetable(tableProtTime,strcat(dataPath,'/protTime.csv'),'WriteRowNames',true)
%% Ploting - protDistGroup
h1=figure; 
top20pProtDist = cellfun(@(x) cellfun(@(y) quantile(y,.75:0.001:.99),x, 'unif',false),protDistGroup,'unif',false);
protDistGroupArray = cellfun(@(y) cell2mat(y),top20pProtDist,'unif',false);
boxPlotCellArray(protDistGroupArray,nameList,1,false,true)
ylabel('Protrusion distance (pix)')
title('Protrusion distance total')
hgexport(h1,strcat(figPath,'/protDist'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/protDist'),'-v7.3')
print(h1,strcat(figPath,'/protDist.tif'),'-dtiff')

tableProtDist=table(protDistGroupArray','RowNames',nameList');
writetable(tableProtDist,strcat(dataPath,'/protDist.csv'),'WriteRowNames',true)

%% Ploting - protVelCondition
% h1=figure; 
% % barPlotCellArray(protVelCondition,nameList,1)
% boxPlotCellArray(retrVelCondition,nameList',1,false,false)
% ylabel('Retraction velocity (um/min)')
% title('Retraction velocity')
% hgexport(h1,strcat(figPath,'/retrVelCondition'),hgexport('factorystyle'),'Format','eps')
% hgsave(h1,strcat(figPath,'/retrVelCondition'),'-v7.3')
% print(h1,strcat(figPath,'/retrVelCondition.tif'),'-dtiff')
% 
% tableRet=table(retrVelCondition,'RowNames',nameList);
% writetable(tableRet,strcat(dataPath,'/retrVelCondition.csv'),'WriteRowNames',true)

end
%% nested function callSteps
% function dataSet = callSteps(ML,exclude,plotYes,labels)
% 
% minLen   = 179;
% outLevel = 8;
% scale    = true;
% % missObs  = 10;
% 
% formatEdgeVelocity(ML,'minLength',minLen,'outLevel',outLevel,'scale',scale);
% 
% %exclude    = cell(1);
% 
% 
% excludeWindowsFromAnalysis(ML,'excBorder',3,'excludeW',exclude);
% 
% % interval = {1:59,60:100,101:179};
% %interval = {1:179};
% [cellData,dataSet] = edgeVelocityQuantification(ML);
% 
% if plotYes
%     nCell = numel(ML.movies_);
%     for iCell = 1:nCell
%         
%         figure
%         set(gcf,'Name',labels{iCell})
%         set(gcf,'number','off')
%         plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,cellData(iCell).protrusionAnalysis.total.percentage,'LineWidth',2)
%         hold on
%         plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,cellData(iCell).retractionAnalysis.total.percentage,'r','LineWidth',2)
%         plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,ones(179,1) - (cellData(iCell).protrusionAnalysis.total.percentage + cellData(iCell).retractionAnalysis.total.percentage),'g','LineWidth',2)
%         xlabel('Time [min]')
%         ylabel('Percentage of the cell Edge')
%         title(labels{iCell})
%         legend({'Protrusion','Retraction','Quiescent'})
%     end
% end
% 
% end