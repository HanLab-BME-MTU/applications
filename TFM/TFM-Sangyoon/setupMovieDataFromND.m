% tfmAnalysis.m run registers the movie and ...
% This is specifically designed for movies taken from Metamorph that
% generates ND files.
%% First Get all project folders
clear
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
        pathProject{ii} = rootFolder;
        groupNames{ii} = finalFolder;
    end
end
%% Now you choose specific dv files of interest per folder
numFolder = numel(pathProject);
orgPath=pwd;
for k=1:numFolder
    cd([pathProject{k} filesep groupNames{k}])
    [fileImgDV{k}, pathImgDV{k}] = uigetfile('*.nd','Select nd files of your interest. Use Ctrl for multiselection',...
    'MultiSelect','on');
end
cd(orgPath)
%% Get the analysis folder (root) and make each analysis folder
rootAnalysis = uigetdir(pwd,'Select the root analysis folder');
% populate groupNames to make absolute paths with root analysis folder
pathAnalysis = cellfun(@(x,y) [y filesep x],groupNames, pathProject, 'unif',false);
%% Make new movieLists for each condition (interactively)
pathDataAll = pathProject;
pathAnalysisAll = pathAnalysis;
%% Print information about pathAnalysis and pathProject
groupNamesCat = strjoin(groupNames);
save([rootAnalysis filesep 'selectedFolders' groupNamesCat '.mat']); % 'rootAnalysis', 'pathDataAll','pathAnalysisAll','pathImgDV','fileImgDV')
%% calibration - This is kind of one-time use. We have to ask users if the calibration files are already made or not
% see if the movie contains 642 or 568 for the bead channel
%% loop - reference frame creation
numConditions = numel(pathDataAll);
iCellRaw=1;
iBeadRaw=2;

thresVariance=0.8; applySobel=true;
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
    curRefDirStruct = dir([curDataPath filesep '*ref*.nd']);  %
    if ~isempty(curRefDirStruct)
        curRefDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirStruct,'unif',false);
    else
        curRefDirStruct = dir([curDataPath filesep '*Ref*.nd']);  %
        if ~isempty(curRefDirStruct)
            curRefDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirStruct,'unif',false);
        else
            curRefDirStruct = dir([curDataPath filesep '*REF*.nd']);  %
            if ~isempty(curRefDirStruct)
                curRefDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirStruct,'unif',false);
            else
                disp('There is no ref, Ref or REF nd files!')
            end
        end
    end
    % There is a case where there are more ref.nd taken than cell.nd. In
    % that case, we don't know what to match. We try to limit the number of
    % nds to the same number as cell nds and see what is the best matching
    % names.
    if numel(curRefDir)>numel(curCellDir)
        curRefDir2 = curCellDir;
        pp=1;
        for ii=1:numel(curCellDir)
            curRefDir2(ii) = curRefDir(pp);
            %Compare the current one with the next one
            curRefDirCandStruct = dir([curRefDir2{ii}(1:end-3) '*.nd']);  %

            if numel(curRefDirCandStruct)==1
                pp=pp+1;
                continue
            elseif numel(curRefDirCandStruct)>1
                % have to find the right one
                curCandDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirCandStruct,'unif',false);
                % in this case, choose the shortest name
                [~,shortestNameID] = min(cellfun(@length,curCandDir));
                curRefDir2(ii) = curCandDir(shortestNameID); 
                pp=pp+1;
                for qq=setdiff(1:numel(curRefDirCandStruct),shortestNameID)
                    pp=pp+1;
                    disp(['The nd file, ' curCandDir{qq} ' is not used due to overlap with ' curRefDir2{ii} '.'])
                end
            end
        end
        curRefDir = curRefDir2;
    end
    
    cellDir{k}=curCellDir;
    refDir{k}=curRefDir;
    %% Check the ref folders
    % averaging 8-12th sections of the stack
    numNDs = numel(curRefDir);
    for ii=1:numNDs
        curRef = curRefDir{ii}; %[curRefDir(ii).folder filesep curRefDir(ii).name];
        refMD = bfImport(curRef,true); % one curRef contains several movies!
%         [meanRefImg,meanRefImgPath] = createBestFocusedImageMD(refMD);
    %     figure, imshow(meanRefImg,[])
        if numel(refMD)==1
            %curRefTifPath{ii}=meanRefImgPath;
        elseif numel(refMD)>1
            nCurMovies=numel(refMD);
            for jj=1:nCurMovies
                curRefMD = refMD(jj);
                curNChan = numel(curRefMD.channels_);
                if curNChan==1
                    chanNames = arrayfun(@(x) x.name_,curRefMD.channels_,'unif',false);
                    folderfilePath = curRefMD.channels_(1).channelPath_;
                    % This is a specific formulation.
                    folderfilePath2=[folderfilePath(1:end-3) '_w1' chanNames{1} '_s' num2str(jj) '.TIF'];
                    curRefTifPath{ii}{jj} = folderfilePath2;
                else
                    %Which one is Cy5?
                    chanNames = arrayfun(@(x) x.name_,curRefMD.channels_,'unif',false);
                    whichChan = cellfun(@(x) strcmp(x,'Cy5'),chanNames);
                    folderfilePath = curRefMD.channels_(whichChan).channelPath_;
                    % This is a specific formulation.
                    folderfilePath2=[folderfilePath(1:end-3) '_w' num2str(find(whichChan)) chanNames{whichChan} '_s' num2str(jj) '.TIF'];
                    curRefTifPath{ii}{jj} = folderfilePath2;
                end
            end
        end
        % 
    end
    refDirTif{k} = curRefTifPath;
end

%% now setting up MD!
for k=numConditions:-1:1
    curDataPath = pathImgDV{k};
    curAnalysisPath = pathAnalysisAll{k};
%     curDir=dir([curDataPath filesep '*.dv']);
%     nameFolders = {curDir.name}';
%     idxRef = contains(nameFolders,'Ref');
%     curRefDir = curDir(idxRef);
%     curCellDir = curDir(~idxRef);
    
    curCellDir = cellfun(@(x) [curDataPath filesep x],fileImgDV{k},'unif',false);
    curRefDir = cellfun(@(x) [x(1:end-6) 'Ref_' x(end-5:end)],curCellDir,'unif',false);
%     curRefDirStruct = dir([curDataPath filesep '*ref*.dv']);  %
%     curRefDir = arrayfun(@(x) [x.folder filesep x.name],curRefDirStruct,'unif',false);
    
    cellDir{k}=curCellDir;
    refDir{k}=curRefDir;
    %% Apply this transforms to the cell channel
    % You don't need to transform ref image because it's bead!!
    pp=0;
    for ii=1:numNDs
        curRawPath=curCellDir{ii}; %[curCellDir(ii).folder filesep curCellDir(ii).name];
        cellMD=bfImport(curRawPath,true); %Now with this Nik's ND, there will be multiple movies per ND
        %identifying what's missing
        sampleMD = cellMD(1);
        if pp==1 && isempty(sampleMD.timeInterval_) && sampleMD.nFrames_>1
            timeLapse = input('What is the time interval in sec?: ');
        else
            timeLapse = 1;
        end

        for jj=1:numel(cellMD)
            pp=pp+1;
            curMD(pp)=cellMD(jj);
            curMD(pp).timeInterval_= timeLapse;
            
            sampleChans = curMD(pp).channels_;
            for qq=1:numel(sampleChans)
                curChan(qq) = sampleChans(qq);
                curChan(qq).emissionWavelength_ = name2wavelength(curChan(qq).name_)*1e9; % in nm
            end
            curMD(pp).setChannel(curChan)
            
            curMD(pp).sanityCheck;
            
            % Save corresponding Ref tiff path
            refDirTifIntegrated{pp}=refDirTif{k}{ii}{jj};
        end
    end
    %%
    refDirTifAll{k}=refDirTifIntegrated;
    PathAnalysis = pathAnalysisAll{k};
    ML = MovieList(curMD,PathAnalysis);
    ML.setPath(PathAnalysis);
    ML.setFilename('movieList.mat');
    ML.sanityCheck;
    ML.save
    clear curMD
    clear ML
    clear channels
    clear refDirTifIntegrated
end
save([rootAnalysis filesep 'selectedFolders' groupNamesCat '.mat'], 'cellDir','refDir', 'refDirTifAll', 'rootAnalysis', 'pathDataAll','pathAnalysisAll','pathImgDV','fileImgDV')
disp('Done making movieDatas and movieList!')
