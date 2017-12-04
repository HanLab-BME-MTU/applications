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
        ClusterInfo.setQueueName('super')
        ClusterInfo.setNNode(1);
        disp(['Currently working on ' curProjectName])
        jobs{k,ii}=batch(@faPackageRun,0,{curMD.getFullPath},'Pool',15, 'profile','nucleus_r2017a');
    end
end
disp('Done submitting all adhesion package jobs!')
