% tfmAnalysis.m run registers the movie and ...
% This is specifically designed for movies taken from OMX-SR
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
save([rootAnalysis filesep 'selectedFolders' groupNamesCat '.mat']); % 'rootAnalysis', 'pathDataAll','pathAnalysisAll','pathImgDV','fileImgDV')
%% calibration - This is kind of one-time use. We have to ask users if the calibration files are already made or not
% see if the movie contains 642 or 568 for the bead channel
disp('Checking which among 568 or 640 is the a channel with higher wavelength...')
sampleMD = bfImport([pathImgDV{1} filesep fileImgDV{1}{1}], 'outputDirectory', [pathAnalysisAll{1} filesep fileImgDV{1}{1}]);
numChannelSample = numel(sampleMD.channels_);
laterChan = sampleMD.channels_(numChannelSample);

mockTransformation = input('Do you just want to apply eye transformation to just re-order the channels? (0/1):');
tFormPath=cell(numChannelSample,1);
if mockTransformation
    % Making transformation - starting from an existing one
    [fileMock, pathMock] = uigetfile('*.mat','Pick a real calibration file.');
    userDataObject = load([pathMock filesep fileMock]);
    userData = userDataObject.userData;
    xForm = userData.xForm;
    xForm.tdata.T = eye(3); %curTavg; % Realized that transformation file is wrong. It's better to just use the raw images
    xForm.tdata.Tinv = eye(3); %curTinvAvg;
    for ii=1:numChannelSample
        % Save curXFormAvg into tFormPath1,2,3...
        tFormPath{ii} = [pathMock num2str(sampleMD.channels_(ii).excitationWavelength_) 'to' num2str(laterChan.excitationWavelength_) '_transform_mock.mat'];
        disp(['Saving mock transformation file ' tFormPath{ii} '....'])
        userData.xForm = xForm;
        save(tFormPath{ii},'userData','xForm');
    end
else
    disp('Doing some real calibration ...')
    calibExist = input(['Do you have calibration file against ' num2str(laterChan.excitationWavelength_) '? (0/1)']);
    %%
    if calibExist

        calibPath{1} = '/project/bioinformatics/Danuser_lab/P01adhesion/raw/kdean-shan/OMX_SR/2017_06_27/Bead_Ref_003.dv';
        % calibPathFor642='/project/bioinformatics/Danuser_lab/P01adhesion/raw/kdean-shan/OMX_SR/2017_10_09_beadRef/TetraSpeck_Cal_004.dv';
        calibPathFor642='/project/bioinformatics/Danuser_lab/P01adhesion/raw/kdean-shan/OMX_SR/2018_01_16_referenceBeads/tetraspeck_003.dv';
        [path1,fname1] = fileparts(calibPath{1});

        if laterChan.excitationWavelength_>560 && laterChan.excitationWavelength_<570
            disp('It''s 568.')
            tFormPath2= [path1 filesep '488to562transform_01.mat'];
            tFormPath3= [path1 filesep '562to562transform_01.mat'];
        elseif laterChan.excitationWavelength_>=640
            %% calibration with the special process
            % get tif for each channel
            % import as MD
            disp('It''s 640. Making files with dv file of calibration bead images')
            %% make the transformation from 01 and 04
            % transformCreationGUI
            [path1,fname1] = fileparts(calibPathFor642);
            tFormPath1= [path1 filesep '488to640transform.mat'];
            tFormPath2= [path1 filesep '562to640transform.mat'];
            tFormPath3= [path1 filesep '640to640transform.mat'];
            disp(['Designate transformation file as ' tFormPath1 '.'])
            disp(['Find the ref bead file at ' meanCaliImgPath '.'])
        end
    else
        disp('Let''s make it!')
        if laterChan.excitationWavelength_>560 && laterChan.excitationWavelength_<570
            disp('It''s 568.')
            tFormPath2= [path1 filesep '488to562transform_01.mat'];
            tFormPath3= [path1 filesep '562to562transform_01.mat'];
        elseif laterChan.excitationWavelength_>=640
            orgPath=pwd;
            cd(rootFolder)
            [fileCalibImgDV, pathCalibImgDV] = uigetfile('*.dv','Select dv files of calibration bead images. Use Ctrl for multiselection',...
                                                            'MultiSelect','on');
            cd(orgPath)
            numCalibFiles = numel(fileCalibImgDV);
            caliMD = bfImport([pathCalibImgDV filesep fileCalibImgDV{1}],true);
            numChanCali = numel(caliMD.channels_);
            xFormAll=cell(numChanCali-1,numCalibFiles);
            for ii=1:numCalibFiles
                disp(['Making the best focused, average image from the stack from ' num2str(ii) 'th file...'])
                caliMD = bfImport([pathCalibImgDV filesep fileCalibImgDV{ii}],true);
                numChanCali = numel(caliMD.channels_);
                pathInImages=cell(numCalibFiles-1,1);
                pathBaseImage = cell(1,1);
                % There are three channels: green(488), red(562) and far red (640)
                for jj=1:numChanCali
                    curChan = caliMD.channels_(jj);

                    curBeadStack = curChan.loadStack(1);
                    maxIntenPerFrame = reshape(max(max(curBeadStack)),[],1);
                    [~,maxIntenFrame]=max(maxIntenPerFrame);
                    minFocusedFrame=max(1,maxIntenFrame-2);
                    maxFocusedFrame=max(caliMD.zSize_,maxIntenFrame+2);

                    curStack = curChan.loadStack(1,'z',minFocusedFrame:maxFocusedFrame);
                    curImg = mean(curStack,3);
                    if jj==numChanCali
                        nameWL=640;
                        meanCaliImgPath = [pathCalibImgDV filesep fileCalibImgDV{ii} '_' num2str(nameWL) '.tif'];
                        pathBaseImage = meanCaliImgPath;
                        baseImg = curImg;
                    else
                        nameWL=curChan.excitationWavelength_;
                        meanCaliImgPath = [pathCalibImgDV filesep fileCalibImgDV{ii} '_' num2str(nameWL) '.tif'];
                        pathInImages{jj}=meanCaliImgPath;
                    end
                    imwrite(uint16(curImg),meanCaliImgPath,'Compression','none')
                end   
                numOtherChans = numChanCali-1;

                for jj=1:numOtherChans
                    curInImg = imread(pathInImages{jj});
                    %Let the user know whats going on
                    waitHan = msgbox('After you click "Ok", you will be shown both of the images you want to align. You must click on several points (Called control points) in both images, so that the two images can be aligned based on these points. The more points you click, the better the resulting transform will be! The minimum # of pairs of points for projective transforms is 4 and for polynomial it is 10. Try to spread the points out evenly over the image area, as this will also improve the resulting transformation. Simply close the control-point selection window when you are finished to continue generating the transform.');
                    uiwait(waitHan);            
    
                    %Call this but scale the images first because it displays to
                    %the whole range            
                    [cpIn,cpBase]= cpselect(mat2gray(curInImg),mat2gray(baseImg),'Wait',true);
    
                    if ~isempty(cpIn) && ~isempty(cpBase)
                        msgbox('Please wait, calculating initial transform...');
                        userData.initXform = cp2tform(cpIn,cpBase,'projective');
        
                        userData.initXformFig = fsFigure(.75);
                        if ~isempty(userData.initXform)
                            subplot(1,2,1)
                        end
                        image(cat(3,mat2gray(baseImg),mat2gray(curInImg),zeros(size(baseImg))));
                        hold on,axis image,axis off
                        title('Overlay, before any transformation. Red: base, Green: Input')    
                        if ~isempty(userData.initXform)
                            subplot(1,2,2)        
                            xIn = imtransform(curInImg,userData.initXform,'XData',[1 size(baseImg,2)],...
                                                                'YData',[1 size(baseImg,1)]);
                            image(cat(3,mat2gray(baseImg),mat2gray(xIn),zeros(size(baseImg))));
                            hold on,axis image,axis off
                            title('Overlay, after initial transformation');
                        end
                        drawnow
                        
                        % Transform refinement
                        userData.xForm = findOptimalXform(baseImg,curInImg,0,'projective',userData.initXform.tdata.T);    
                        %Show the pre- and post-final transform alignment, if one was created
                        if isfield(userData, 'initXformFig') && ishandle(userData.initXformFig)
                            delete(userData.initXformFig) 
                        end
                        xForm = userData.xForm;
                        xFormAll{jj,ii}=xForm;
                        userData.initXformFig = fsFigure(.75);
                        xIn = imtransform(curInImg,userData.xForm,'XData',[1 size(baseImg,2)],...
                                                        'YData',[1 size(baseImg,1)]);
                        image(cat(3,mat2gray(baseImg),mat2gray(xIn),zeros(size(baseImg))));
                        hold on,axis image,axis off
                        title('Overlay, after refined transformation');
                        savePath = [pathCalibImgDV filesep num2str(caliMD.channels_(jj).excitationWavelength_) 'to' num2str(laterChan.excitationWavelength_) '_transform_' num2str(ii) '.mat'];
                        disp(['Saving transformation file ' savePath '....'])
                        save(savePath,'-struct','userData','xForm');
                        
                    end                
                end
                % Making transformation ...
    %             disp(['Designate transformation file as ' tFormPath1 '.'])
    %             disp(['Find the ref bead file at ' meanCaliImgPath '.'])
    %             h=transformCreationGUI; % base img: 640, input img: 488, this time I used an existing transform in 2017_02_10
    %             uiwait(h); 
            end
            %% Now we have all the xForms, xFormAll. Now it's a time to average it
            % per channel
            for jj=1:numOtherChans
                try
                    curXFormAll = xFormAll(jj,:);
                    curXFormAll = curXFormAll(~cellfun(@isempty,curXFormAll));
                    curXFormAvg = curXFormAll{1}; %initializing
                    curTs = cellfun(@(x) x.tdata.T,curXFormAll,'unif',false);
                    curTsMat = cat(3,curTs{:});
                    curTinvs = cellfun(@(x) x.tdata.Tinv,curXFormAll,'unif',false);
                    curTinvsMat = cat(3,curTinvs{:});
                    curTavg = mean(curTsMat,3);
                    curTinvAvg = mean(curTinvsMat,3);

                    curXFormAvg.tdata.T = curTavg; % Realized that transformation file is wrong. It's better to just use the raw images
                    curXFormAvg.tdata.Tinv = curTinvAvg;
%                     xFormAvg{jj}=curXFormAvg;
                    % Save curXFormAvg into tFormPath1,2,3...
                    tFormPath{jj} = [pathCalibImgDV num2str(caliMD.channels_(jj).excitationWavelength_) 'to' num2str(laterChan.excitationWavelength_) '_transform_avg.mat'];
                    disp(['Saving transformation file ' tFormPath{jj} '....'])
                    xForm=curXFormAvg;
                    userData.xForm = xForm;
                    save(tFormPath{jj},'userData','xForm');
                catch
                    [tempFile,tempPath] = uigetfile('*.mat');
                    tFormPath{jj} = [tempPath filesep tempFile];
                end
            end
            %% For the base channel
            jj=jj+1;
            try
                curXFormAvg.tdata.T = eye(3);
                curXFormAvg.tdata.Tinv = eye(3);
                tFormPath{jj} = [pathCalibImgDV num2str(caliMD.channels_(jj).excitationWavelength_) 'to' num2str(laterChan.excitationWavelength_) '_transform_avg.mat'];
                disp(['Saving transformation file ' tFormPath{jj} '....'])
                xForm=curXFormAvg;
                userData.xForm = xForm;
                save(tFormPath{jj},'userData','xForm');
            catch
                [tempFile,tempPath] = uigetfile('*.mat');
                tFormPath{jj} = [tempPath filesep tempFile];
            end
        end
    end
end
%% This is temporary. We'll apply the transformation to the raw movies (not making the movie data file). It's kind of overriding the existing ones.
% numConditions = numel(pathDataAll);

% for k=numConditions:-1:1
%     curDataPath = pathDataAll{k};
%     curAnalysisPath = pathAnalysisAll{k};
% %     curDir=dir([curDataPath filesep '*.dv']);
% %     nameFolders = {curDir.name}';
% %     idxRef = contains(nameFolders,'Ref');
% %     curRefDir = curDir(idxRef);
% %     curCellDir = curDir(~idxRef);
%     
%     curCellDir = cellfun(@(x) [curDataPath filesep x],fileImgDV{k},'unif',false);
%     
%     cellDir{k}=curCellDir;
%     % You don't need to transform ref image because it's bead!!
%     numCells = numel(curCellDir);
%     for ii=1:numCells
%         curRawPath=curCellDir{ii}; %[curCellDir(ii).folder filesep curCellDir(ii).name];
%         cellMD=bfImport(curRawPath,true);
%         % if cellMD has z-stack with one time frame, we have to compress it to one frame image
%         if cellMD.zSize_>1 && cellMD.nFrames_==1
%             midPixelsAllFrames = cell(cellMD.zSize_,1);
%             for iiChan=1:numel(cellMD.channels_)
%                 curChan=cellMD.channels_(iiChan);
%                 curChanStack = curChan.loadStack(1);
%                 for jj=1:cellMD.zSize_
%                     curImageFrame = curChanStack(:,:,jj);
%                     curImageFrameSorted = sort(curImageFrame(:),'descend');
%                     midPixelsAllFrames{jj} = curImageFrameSorted(100:300);
%                 end
%                 meanMidInten = cellfun(@mean,midPixelsAllFrames);
%                 % Take top five frames
%                 [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
%                 averagingRange = meanMidIntenIDs(1:5);
% 
%                 curImg = mean(curChanStack(:,:,averagingRange),3); 
%                 [path2,file2] = fileparts(curChan.channelPath_);
%                 %Create a new folder with one frame
%                 path3 = [path2 filesep file2 'OneFrame' filesep num2str(curChan.excitationWavelength_)];
%                 mkdir(path3)
%                 meanImgPath = [path3 filesep num2str(curChan.excitationWavelength_) '.tif'];
%                 imwrite(uint16(curImg),meanImgPath,'Compression','none')
%                 curChanNew(iiChan)=Channel(path3);
%                 curChanNew(iiChan).excitationWavelength_ = curChan.excitationWavelength_;
%                 curChanNew(iiChan).emissionWavelength_ = curChan.emissionWavelength_;
%             end
%             curAnalysisFolder3 = [path2 filesep file2 'OneFrame'];
%             cellMD2 = MovieData(curChanNew,curAnalysisFolder3);
%             cellMD2.setPath(curAnalysisFolder3);
%             cellMD2.setFilename([file2 '.mat']);
%             cellMD2.numAperture_=cellMD.numAperture_;
%             cellMD2.camBitdepth_=cellMD.camBitdepth_;
%             cellMD2.timeInterval_ = cellMD.timeInterval_;
%             cellMD2.pixelSize_= cellMD.pixelSize_; % 60x x 1.8x (new objective config.)
%             cellMD2.sanityCheck;
%             cellMD2.save
%             cellMD = cellMD2;
%         end        
%         %Get the indices of any previous tranformation correction processes
%         iProc = cellMD.getProcessIndex('TransformationProcess', 1, false);
% 
%         %If the process doesn't exist, create it with default settings.
%         if isempty(iProc)
%             iProc = numel(cellMD.processes_)+1;
%             cellMD.addProcess(TransformationProcess(cellMD,cellMD.outputDirectory_));
%         end
% 
%         transfProc = cellMD.getProcess(iProc);
%         p = transfProc.funParams_;
%         curNumChan = numel(cellMD.channels_);
%         if curNumChan==2
%             p.ChannelIndex=[1 2]; %Cell channel
%             firstChan=cellMD.channels_(1);
%             %There is a case where experiment was doen with 488 and 640.
%             %Then 488to640transform.mat should be assigned to channel 1
%             if firstChan.excitationWavelength_>480 && firstChan.excitationWavelength_<500
%                 p.TransformFilePaths={tFormPath{1}, tFormPath{3}};
%             else
%                 p.TransformFilePaths={tFormPath{2}, tFormPath{3}};
%             end
%         elseif curNumChan==3
%             p.ChannelIndex=[1 2 3]; %Cell channel
%             p.TransformFilePaths={tFormPath{1}, tFormPath{2}, tFormPath{3}};
%         end            
%         p.TransformMasks=false;
%         transfProc.setPara(p)
%         transfProc.run
%         cellMD.save
%         curCellDir{ii} = [cellMD.movieDataPath_ filesep cellMD.movieDataFileName_];
%     end
% end
%% loop - now setting up MD!
numConditions = numel(pathDataAll);
iCellRaw=1;
iBeadRaw=2;

for k=numConditions:-1:1
    curDataPath = pathImgDV{k};
    curAnalysisPath = pathAnalysisAll{k};
%     curDir=dir([curDataPath filesep '*.dv']);
%     nameFolders = {curDir.name}';
%     idxRef = contains(nameFolders,'Ref');
%     curRefDir = curDir(idxRef);
%     curCellDir = curDir(~idxRef);
    
    curCellDir = cellfun(@(x) [curDataPath filesep x],fileImgDV{k},'unif',false);
%     curRefDir = cellfun(@(x) [x(1:end-6) 'Ref_' x(end-5:end)],curCellDir,'unif',false);
    curRefDirStruct = dir([curDataPath filesep '*ref*.dv']);  %
    curRefDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirStruct,'unif',false);
    
    cellDir{k}=curCellDir;
    refDir{k}=curRefDir;
    %% Check the ref folders
    % averaging 8-12th sections of the stack
    numCells = numel(curRefDir);
    for ii=1:numCells
        curRef = curRefDir{ii}; %[curRefDir(ii).folder filesep curRefDir(ii).name];
        refMD = bfImport(curRef,true); % 'outputDirectory', [pathAnalysisAll{1} filesep fileImgDV{1}{1}]);
        curRefBeadChan = refMD.channels_(end);
        if refMD.zSize_>1
            % find the best focus
            curRefBeadStack = curRefBeadChan.loadStack(1);
            % I will take 100th to 200th brightest pixels to guess the best
            % focus (usually the very brightest point is from one
            % extraordinary bead) - SH 20171010
            % Get 100th to 200th pixels
            midPixelsAllFrames = cell(refMD.zSize_,1);
            for jj=1:refMD.zSize_
                curRefImageFrame = curRefBeadStack(:,:,jj);
                curRefImageFrameSorted = sort(curRefImageFrame(:),'descend');
                midPixelsAllFrames{jj} = curRefImageFrameSorted(100:300);
            end
            meanMidInten = cellfun(@mean,midPixelsAllFrames);
            % Take top five frames
            [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
            averagingRange = meanMidIntenIDs(1:5);
%             maxProf = reshape(max(max(curRefBeadStack)),[],1);
%             [~,maxIntenFrame]=max(maxProf);
%             minFocusedFrame=max(1,maxIntenFrame-2);
%             maxFocusedFrame=max(refMD.zSize_,maxIntenFrame+2);
            
            curRefBeadStack = curRefBeadChan.loadStack(1,'z',averagingRange);
            meanRefImg = mean(curRefBeadStack,3); 
        else
            meanRefImg = curRefBeadChan.loadImage(1);
        end
    %     figure, imshow(meanRefImg,[])
        % store it somewhere
        [path1,fname1] = fileparts(curRef);
        meanRefImgPath = [path1 filesep fname1 '.tif'];
        imwrite(uint16(meanRefImg),meanRefImgPath,'Compression','none')
        curRefTifPath{ii}=meanRefImgPath;
        % 
    end
    refDirTif{k} = curRefTifPath;

    %% Apply this transforms to the cell channel
    % You don't need to transform ref image because it's bead!!
    for ii=1:numCells
        curRawPath=curCellDir{ii}; %[curCellDir(ii).folder filesep curCellDir(ii).name];
        cellMD=bfImport(curRawPath,true);
        % if cellMD has z-stack with one time frame, we have to compress it to one frame image
        if cellMD.zSize_>1 && cellMD.nFrames_==1
            midPixelsAllFrames = cell(cellMD.zSize_,1);
            for iiChan=1:numel(cellMD.channels_)
                curChan=cellMD.channels_(iiChan);
                curChanStack = curChan.loadStack(1);
                for jj=1:cellMD.zSize_
                    curImageFrame = curChanStack(:,:,jj);
                    curImageFrameSorted = sort(curImageFrame(:),'descend');
                    midPixelsAllFrames{jj} = curImageFrameSorted(100:300);
                end
                meanMidInten = cellfun(@mean,midPixelsAllFrames);
                % Take top five frames
                [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
                averagingRange = meanMidIntenIDs(1:5);

                curImg = mean(curChanStack(:,:,averagingRange),3); 
                [path2,file2] = fileparts(curChan.channelPath_);
                %Create a new folder with one frame
                path3 = [path2 filesep file2 'OneFrame' filesep num2str(curChan.excitationWavelength_)];
                mkdir(path3)
                meanImgPath = [path3 filesep num2str(curChan.excitationWavelength_) '.tif'];
                imwrite(uint16(curImg),meanImgPath,'Compression','none')
                curChanNew(iiChan)=Channel(path3);
                curChanNew(iiChan).excitationWavelength_ = curChan.excitationWavelength_;
                curChanNew(iiChan).emissionWavelength_ = curChan.emissionWavelength_;
            end
            curAnalysisFolder3 = [path2 filesep file2 'OneFrame'];
            cellMD2 = MovieData(curChanNew,curAnalysisFolder3);
            cellMD2.setPath(curAnalysisFolder3);
            cellMD2.setFilename([file2 '.mat']);
            cellMD2.numAperture_=cellMD.numAperture_;
            cellMD2.camBitdepth_=cellMD.camBitdepth_;
            cellMD2.timeInterval_ = cellMD.timeInterval_;
            cellMD2.pixelSize_= cellMD.pixelSize_; % 60x x 1.8x (new objective config.)
            cellMD2.sanityCheck;
            cellMD2.save
            cellMD = cellMD2;
        end        
        %Get the indices of any previous tranformation correction processes
        iProc = cellMD.getProcessIndex('TransformationProcess', 1, false);

        %If the process doesn't exist, create it with default settings.
        if isempty(iProc)
            iProc = numel(cellMD.processes_)+1;
            cellMD.addProcess(TransformationProcess(cellMD,cellMD.outputDirectory_));
        end

        transfProc = cellMD.getProcess(iProc);
        p = transfProc.funParams_;
        curNumChan = numel(cellMD.channels_);
        if curNumChan==2
            p.ChannelIndex=[1 2]; %Cell channel
            firstChan=cellMD.channels_(1);
            %There is a case where experiment was doen with 488 and 640.
            %Then 488to640transform.mat should be assigned to channel 1
%             if firstChan.excitationWavelength_>480 && firstChan.excitationWavelength_<500
                p.TransformFilePaths={tFormPath{1}, tFormPath{2}};
%             else
%                 p.TransformFilePaths={tFormPath{2}, tFormPath{3}};
%             end
        elseif curNumChan==3
            p.ChannelIndex=[1 2 3]; %Cell channel
            p.TransformFilePaths={tFormPath{1}, tFormPath{2}, tFormPath{3}};
        end            
        p.TransformMasks=false;
        transfProc.setPara(p)
        transfProc.run
        cellMD.save
        curCellDir{ii} = [cellMD.movieDataPath_ filesep cellMD.movieDataFileName_];
    end
    %% Register channels appropriately to a new MD
%     orgDir=pwd;
    iNewBead=1;
    iNewCell=2;
    PathAnalysis = curAnalysisPath;

    for ii=1:numCells
        [~, curMDfolder]=fileparts(curCellDir{ii}); %[curCellDir(ii).folder filesep curCellDir(ii).name]);
%         cd([curRawPath filesep curMDfolder])
        curMDname = curCellDir{ii}; %[curRawPath filesep curMDfolder filesep curMDfolder '.mat'];
        cellMD=MovieData.load(curMDname);
        iProc = cellMD.getProcessIndex('TransformationProcess', 1, false);
        transfProc = cellMD.getProcess(iProc);
        oldCh = cellMD.channels_;
        % Now making Channel (Bead=1)
        newCh(iNewBead)=Channel(transfProc.outFilePaths_{end});
        newCh(iNewBead).excitationWavelength_ = oldCh(end).excitationWavelength_;
        newCh(iNewBead).emissionWavelength_ = oldCh(end).emissionWavelength_;
        % Cell=2
        newCh(iNewCell)=Channel(transfProc.outFilePaths_{iCellRaw});
        newCh(iNewCell).excitationWavelength_ = oldCh(iCellRaw).excitationWavelength_;
        newCh(iNewCell).emissionWavelength_ = oldCh(iCellRaw).emissionWavelength_;
        if numel(oldCh)>2
            newCh(iNewCell+1)=Channel(transfProc.outFilePaths_{iCellRaw+1});
            newCh(iNewCell+1).excitationWavelength_ = oldCh(iCellRaw+1).excitationWavelength_;
            newCh(iNewCell+1).emissionWavelength_ = oldCh(iCellRaw+1).emissionWavelength_;
        end

        % new movieData!
        % Constructor needs an array of channels and an output directory (for analysis)
        analysisFolder = [PathAnalysis filesep curMDfolder];
        if ~exist(analysisFolder,'dir')
            mkdir(analysisFolder)
        end
        curMD(ii) = MovieData(newCh,analysisFolder);
        % Set the path where to store the MovieData object.
        curMD(ii).setPath(analysisFolder);
        curMD(ii).setFilename('movieData.mat');

        % Set some additional movie properties
        curMD(ii).numAperture_=cellMD.numAperture_;
        curMD(ii).camBitdepth_=cellMD.camBitdepth_;
        curMD(ii).timeInterval_ = cellMD.timeInterval_;
        curMD(ii).pixelSize_= cellMD.pixelSize_; % 60x x 1.8x (new objective config.)
    %     curMD(ii).nFrames_= cellMD.nFrames_; % 
    %     curMD(ii).imSize_= cellMD.imSize_; % 
    %     curMD(ii).zSize_= cellMD.zSize_; % 
    %     curMD(ii).pixelSizeZ_= cellMD.pixelSizeZ_; % 
        curMD(ii).acquisitionDate_= cellMD.acquisitionDate_; % 60x x 1.8x (new objective config.)
        if ii<numCells
            curMD(ii).timeInterval_= 2;
        else
            curMD(ii).timeInterval_= 3;
        end
        curMD(ii).sanityCheck;
        % Save the movie
        curMD(ii).save
        clear newCh
    end
%     cd(orgDir)
    %%
    ML = MovieList(curMD,PathAnalysis);
    ML.setPath(PathAnalysis);
    ML.setFilename('movieList.mat');
    ML.sanityCheck;
    ML.save
    clear curMD
    clear ML
    clear channels
end
%% Load movieLists for each condition
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
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
    curRefDir=refDirTif{k}; %dir([curDataPath filesep '*Ref*.tif']);
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


