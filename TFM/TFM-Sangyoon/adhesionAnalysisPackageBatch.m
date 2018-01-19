%% Asking initializing
initYes=input('Do you want to start from scratch?(1/0):');
if initYes
    %% First Get all project folders
    dataFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~dataFolderSelectionDone
        ii=ii+1;
        curPathProject = uigetdir(rootFolder,'Select data folder (Type zero when no more)');
        [rootFolder,finalFolder] = fileparts(curPathProject);
        if strcmp(finalFolder,'0')
            dataFolderSelectionDone=true;
        else
            pathProject{ii} = curPathProject;
            groupNames{ii} = finalFolder;
        end
    end
    %% Now you choose specific dv files of interest per folder
    numFolder = numel(pathProject);
    orgPath=pwd;
    for k=1:numFolder
        cd(pathProject{k})
        [fileImgDV{k}, pathImgDV{k}] = uigetfile('*.dv','Select dv files of your interest. Use Ctrl for multiselection',...
        'MultiSelect','on');
    end
    cd(orgPath)
    %% Get the analysis folder (root) and make each analysis folder
    rootAnalysis = uigetdir(pwd,'Select the root analysis folder');
    % populate groupNames to make absolute paths with root analysis folder
    pathAnalysis = cellfun(@(x) [rootAnalysis filesep x],groupNames,'unif',false);
    %% Make new movieLists for each condition (interactively)
    pathDataAll = pathProject;
    pathAnalysisAll = pathAnalysis;
    %% Print information about pathAnalysis and pathProject
    groupNamesCat = strjoin(groupNames);
    save([rootAnalysis filesep 'selectedFolders' groupNamesCat '.mat'], 'rootAnalysis', 'pathDataAll','pathAnalysisAll')
    %% Set up MLs with dvs
    for k=1:numFolder
        curDataPath = pathDataAll{k};
        curAnalysisPath = pathAnalysisAll{k};
        numCells = numel(fileImgDV{k});
        curCellDir = cellfun(@(x) [curDataPath filesep x],fileImgDV{k},'unif',false);
        curCellAnalDir = cellfun(@(x) [curAnalysisPath filesep x],fileImgDV{k},'unif',false);
        for ii=1:numCells
            curRawPath=curCellDir{ii}; %[curCellDir(ii).folder filesep curCellDir(ii).name];
            outputDir=curCellAnalDir{ii};
            cellMD(ii)=bfImport(curRawPath,true, 'outputDirectory', outputDir);
            cellMD(ii).sanityCheck;
            % Save the movie
            cellMD(ii).save
        end
        ML = MovieList(cellMD,curAnalysisPath);
        ML.setPath(curAnalysisPath);
        ML.setFilename('movieList.mat');
        ML.sanityCheck;
        ML.save
        clear cellMD
    end
else
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
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummary/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary/Data'];
mkdir(dataPath)
%% running through node (MDCS)
ClusterInfo.setEmailAddress('sjhan@mtu.edu');
% ClusterInfo.getWallTime
% create job script m-file 
for k=1:numConditions
    curML = MLAll(k);
    curNumCells = numel(curML.movieDataFile_);
    for ii=1:curNumCells
        curMD = curML.getMovie(ii);
        [~,curProjectName] = fileparts(curMD.movieDataPath_);
%         ClusterInfo.setQueueName('super')
        ClusterInfo.setQueueName('256GB')
        ClusterInfo.setNNode(1);
        disp(['Currently working on ' curProjectName])
%         jobs{k,ii}=batch(@faPackageRun,0,{curMD.getFullPath},'Pool',15, 'profile','nucleus_r2017a');
        faPackageRun(curMD.getFullPath);
    end
end
disp('Done submitting all adhesion package jobs!')
