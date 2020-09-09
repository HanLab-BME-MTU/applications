% setupTFMPackageDefault.m sets up TFM package for chosen set of MLs.
% It is recommended that you run setupMovieDataFromOMXDV.m first. 

%% Read selectedFolders.mat
clear refDirTifAllCell
try
    [pathAnalysisAll, MLNames, groupNames,usedSelectedFoldersMat,specificName,refDirTifAllCell] = chooseSelectedFolders;
catch
    [pathAnalysisAll, MLNames, groupNames,usedSelectedFoldersMat] = chooseSelectedFolders;
end

%% Load movieLists for each condition
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end

%% TFM parameter
prompt = {'Enter YoungModulus (in Pa):','Enter solMethodBEM(e.g.QR or 1NormReg:',...
    'Enter useLcurve (true/false):','Enter FastBEM or FTTC',...
    'Enter lcorner or optimal:','Enter reg parameter:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'5000','QR','true','FastBEM','lcorner','1e-6'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
YoungModulus = str2double(answer{1}); % Pa
solMethodBEM=answer{2};
useLcurve = str2num(answer{3});
bemOrFTTC = answer{4};
lcornerOptimal = answer{5};
regParam = str2double(answer{6});
% save
save([rootAnalysis filesep 'ForceVariables' groupNamesCat '.mat'], 'YoungModulus',...
    'solMethodBEM','useLcurve','bemOrFTTC','lcornerOptimal','regParam')
%% Setting up each movie for TFM
MLgroup=MLAll;
numML=numel(MLgroup);
for k=1:numML
    curML = MLgroup(k);
    curNumCells = numel(curML.movieDataFile_);

    curDataPath = pathDataAll{k};
    curRefDir=refDirTifAllCell{k}; %dir([curDataPath filesep '*Ref*.tif']);
    for ii=1:curNumCells
        curMD = curML.getMovie(ii);
        %% Create TFM package and retrieve package index
        iPack=  curMD.getPackageIndex('TFMPackage');
        if isempty(iPack)
            curMD.addPackage(TFMPackage(curMD));
            iPack=  curMD.getPackageIndex('TFMPackage');
        end
        tfmPack = curMD.getPackage(iPack);
        %% Create first process
        if isempty(tfmPack.processes_{1})
            tfmPack.createDefaultProcess(1)
        end
        curSDCProc = tfmPack.getProcess(1);
        if 1 %curSDCProc.updated_ %&& curSDCProc.procChanged_
            params = curSDCProc.funParams_;
            refImgPath = curRefDir{ii}; %[curRefDir(ii).folder filesep curRefDir(ii).name];
            params.referenceFramePath=refImgPath;
            params.usfac=100;
            curSDCProc.setPara(params);
        end
        %% ROI for cell to reduce memory burden
        if isempty(curMD.roiMask)
            sampleCellImgFirst = (curMD.channels_(2).loadImage(1));
            sampleCellImgD1=double(sampleCellImgFirst);
            sampleCellImgNorm1=(sampleCellImgD1-min(sampleCellImgD1(:)))/(max(sampleCellImgD1(:))-min(sampleCellImgD1(:)));
            nFrames=curMD.nFrames_;
            sampleCellImgLast = (curMD.channels_(2).loadImage(nFrames));
            sampleCellImgDLast=double(sampleCellImgLast);
            sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
            sampleBeadImg = curMD.channels_(1).loadImage(1);
            sampleBeadImgD=double(sampleBeadImg);
            sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
            compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
            compImg(:,:,3)=sampleBeadImgNorm;
            compImg(:,:,2)=sampleCellImgNorm2;
            compImg(:,:,1)=sampleCellImgNorm1;
            h12=figure; imshow(compImg), hold on
            disp(['Draw rectangle for ROI for ' curMD.movieDataPath_ '.'])

            h=imrect;
            ROI_rect = wait(h);
            roiMask=createMask(h);

            % Save it as ROI mask associated with MD
            roiPath=[curMD.outputDirectory_ filesep 'roiMask.tif'];
            imwrite(roiMask,roiPath);
            curMD.setROIMaskPath(roiPath);
            % maskArray = imread(MD.roiMaskPath_);
            curMD.roiMask=roiMask;
            close(h12)
%         else
%             sampleCellImgFirst = (curMD.channels_(2).loadImage(1));
%             sampleCellImgD1=double(sampleCellImgFirst);
%             sampleCellImgNorm1=(sampleCellImgD1-min(sampleCellImgD1(:)))/(max(sampleCellImgD1(:))-min(sampleCellImgD1(:)));
%             nFrames=curMD.nFrames_;
%             sampleCellImgLast = (curMD.channels_(2).loadImage(nFrames));
%             sampleCellImgDLast=double(sampleCellImgLast);
%             sampleCellImgNorm2=(sampleCellImgDLast-min(sampleCellImgDLast(:)))/(max(sampleCellImgDLast(:))-min(sampleCellImgDLast(:)));
%             sampleBeadImg = curMD.channels_(1).loadImage(1);
%             sampleBeadImgD=double(sampleBeadImg);
%             sampleBeadImgNorm=(sampleBeadImgD-min(sampleBeadImgD(:)))/(max(sampleBeadImgD(:))-min(sampleBeadImgD(:)));
%             compImg = zeros(size(sampleBeadImg,1),size(sampleBeadImg,2),3);
%             compImg(:,:,3)=sampleBeadImgNorm;
%             compImg(:,:,2)=sampleCellImgNorm2;
%             compImg(:,:,1)=sampleCellImgNorm1;
%             h12=figure; imshow(compImg), hold on
%             disp(['Draw rectangle for ROI for ' curMD.movieDataPath_ '.'])
%             roiMask=curMD.roiMask;
%             boundROI=bwboundaries(roiMask);
%             hold on, 
%             try
%                 plot(boundROI{1}(:,2),boundROI{1}(:,1),'w')
%                 h=imrect;
%             catch                
%                 h=imrect;
%             end
%             ROI_rect = wait(h);
%             roiMask=createMask(h);
% 
%             % Save it as ROI mask associated with MD
%             roiPath=[curMD.outputDirectory_ filesep 'roiMask.tif'];
%             imwrite(roiMask,roiPath);
%             curMD.setROIMaskPath(roiPath);
%             % maskArray = imread(MD.roiMaskPath_);
%             curMD.roiMask=roiMask;
%             close(h12)
        end
        %% displacement 
        if isempty(tfmPack.processes_{2})
            tfmPack.createDefaultProcess(2)
        end
        if 1 %tfmPack.getProcess(2).updated_
            params = tfmPack.getProcess(2).funParams_;
%             lastImgFileName=curMD.channels_(1).getImageFileNames(end);
            params.referenceFramePath = refImgPath; %[curMD.channels_(1).channelPath_ filesep lastImgFileName{1}]; % last frame
            params.alpha = 0.05;
            params.minCorLength = 19;
            params.addNonLocMaxBeads=true;
            params.maxFlowSpeed = 40; %for 16 kPa
            params.highRes = true;
            params.useGrid = false; % PIV suite
            params.mode = 'accurate';
            params.noFlowOutwardOnBorder=1;
            params.trackSuccessively=false;
            tfmPack.getProcess(2).setPara(params);
            curMD.save;
        end
        %% Create third process and run
        if isempty(tfmPack.processes_{3})
            tfmPack.createDefaultProcess(3)
        end
        params = tfmPack.getProcess(3).funParams_;
        params.outlierThreshold = 2;
        params.fillVectors=true;
        tfmPack.getProcess(3).setPara(params);
        %% Create force reconstruction process and run
        if isempty(tfmPack.processes_{4})
            tfmPack.createDefaultProcess(4)
        end
        params = tfmPack.getProcess(4).funParams_;

        params.YoungModulus = YoungModulus;
        params.regParam = regParam;
        params.method = bemOrFTTC;
        % params.solMethodBEM = '1NormReg';
        params.solMethodBEM = solMethodBEM;%'1NormRegLaplacian';
        params.useLcurve = useLcurve;
        params.basisClassTblPath = ['/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/TFM basis functions/' 'basisClass' strrep(num2str(YoungModulus/1000),'.','_') 'kPa' num2str(5) 'pix.mat'];

        tfmPack.getProcess(4).setPara(params);
        curMD.save;
        %% Create strain energy calculation process and run
        if isempty(tfmPack.processes_{5})
            tfmPack.createDefaultProcess(5)
        end
        params = tfmPack.getProcess(5).funParams_;

        tfmPack.getProcess(5).setPara(params);
        curMD.save;
        clear curMD
    end
end

close all

%% running through node (MDCS)
% ClusterInfo.setEmailAddress('sjhan@mtu.edu');
% ClusterInfo.getWallTime
% create job script m-file 
wayOfRun = input('Run by this computer(1), or run by job submission (2), or stop (0 or empty)?:');
if ~isempty(wayOfRun) && wayOfRun>0 && wayOfRun<3
    for k=1:numML
        curML = MLgroup(k);
        curNumCells = numel(curML.movieDataFile_);
        for ii=1:curNumCells
            curMD = curML.getMovie(ii);
            disp(['Currently working on ' curMD.movieDataPath_])
            if wayOfRun==2 %strcmp(bemOrFTTC,'FastBEM') || curMD.nFrames_>10
                [~,curProjectName] = fileparts(curMD.movieDataPath_);
        %         curProjectName = [curProjectName(end-21:end-20) curProjectName(end-2:end)];
        %         ClusterInfo.setProjectName(curProjectName); 
                curRoiMask=curMD.roiMask;
                boundROI=bwboundaries(curRoiMask);
                curW=max(boundROI{1}(:,2))-min(boundROI{1}(:,2));
                curH=max(boundROI{1}(:,1))-min(boundROI{1}(:,1));
                iPack=  curMD.getPackageIndex('TFMPackage');
                tfmPack = curMD.getPackage(iPack);
                params = tfmPack.getProcess(4).funParams_;
                if curW*curH>700000 && strcmp(params.method,'FastBEM')
                    ClusterInfo.setQueueName('256GB')
                else
                    ClusterInfo.setQueueName('super')
                end
                ClusterInfo.setNNode(1);
        %         ClusterInfo.setPrivateKeyFile
                disp(['Currently working on ' curProjectName])
                jobs{k,ii}=batch(@tfmRun,0,{curMD.getFullPath},'profile','nucleus_r2017a'); %'Pool',15);
            elseif wayOfRun==1 %strcmp(bemOrFTTC,'FTTC') && curMD.nFrames_<=10
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
else
    movieSelectorGUI
end
%% saving essential parameters for next MD processing


