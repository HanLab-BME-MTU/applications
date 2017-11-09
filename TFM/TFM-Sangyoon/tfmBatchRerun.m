% This function is used when you changed the parameters in the TFM package 
%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.');
selectedFolders=load([pathSFolders filesep fileSFolders]);
pathAnalysisAll=selectedFolders.pathAnalysisAll;
numConditions = numel(pathAnalysisAll);
%% Load movieLists for each condition
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
numML=numel(MLAll);
%% running through node (MDCS)
ClusterInfo.setEmailAddress('sjhan@mtu.edu');
% ClusterInfo.getWallTime
% create job script m-file 
for k=1:numML
    curML = MLAll(k);
    curNumCells = numel(curML.movieDataFile_);
    for ii=1:curNumCells
        curMD = curML.getMovie(ii);
        disp(['Currently working on ' curMD.movieDataPath_])
        iPack=  curMD.getPackageIndex('TFMPackage');
        tfmPack = curMD.getPackage(iPack);
        params = tfmPack.getProcess(4).funParams_;
        bemOrFTTC=params.method;
        if strcmp(bemOrFTTC,'FastBEM') || curMD.nFrames_>10
            [~,curProjectName] = fileparts(curMD.movieDataPath_);
    %         curProjectName = [curProjectName(end-21:end-20) curProjectName(end-2:end)];
    %         ClusterInfo.setProjectName(curProjectName); 
            curRoiMask=curMD.roiMask;
            if isempty(curRoiMask)
                curRoiMask=imread(curMD.roiMaskPath_);
            end
            boundROI=bwboundaries(curRoiMask);
            curW=max(boundROI{1}(:,2))-min(boundROI{1}(:,2));
            curH=max(boundROI{1}(:,1))-min(boundROI{1}(:,1));
            if curW*curH>700000 && strcmp(params.method,'FastBEM')
                ClusterInfo.setQueueName('256GB')
            else
                ClusterInfo.setQueueName('super')
            end
            ClusterInfo.setNNode(1);
    %         ClusterInfo.setPrivateKeyFile
            disp(['Currently working on ' curProjectName])
            jobs{k,ii}=batch(@tfmRun,0,{curMD.getFullPath},'Pool',15, 'profile','nucleus_r2017a');
        elseif strcmp(bemOrFTTC,'FTTC') && curMD.nFrames_<=10
            iPack=  curMD.getPackageIndex('TFMPackage');
            %% Run the stage drift correction
            if ~curMD.getPackage(iPack).getProcess(1).success_
                curMD.getPackage(iPack).getProcess(1).run();
            end
            if ~curMD.getPackage(iPack).getProcess(2).success_
                curMD.getPackage(iPack).getProcess(2).run();
            end
            if ~curMD.getPackage(iPack).getProcess(3).success_
                curMD.getPackage(iPack).getProcess(3).run();
            end
            if ~curMD.getPackage(iPack).getProcess(4).success_
                curMD.getPackage(iPack).getProcess(4).run();
            end
            if ~curMD.getPackage(iPack).getProcess(5).success_
                curMD.getPackage(iPack).getProcess(5).run();
            end
        else
            disp('You must choose either FastBEM or FTTC for bemOrFTTC!')
        end
    end
end