function [validationAccuracy,validationAccuracyCurCell, validationAccuracyCross,validationAccuracyCrossUpdated,...
    C,CCurCell,CCross,CUpdated, CCrossUpdated]=testCrossCellClassification(pathForColocalization,varargin)
% 
% this function classifyNascentAdhesionTracks(pathForColocalization,varargin) reads from Colocalization folder 
% and classify NA tracks based on fluorescence signal amplitude evolution
% and distance to the edge. 
% group 1: NAs that form at the edge but say there as the edge protrude and
% turn over:
% group 2: NAs that form and mature as the edge protrude: 
% group 3: NAs that move along the edge and quickly turn over

%% input reading
ip =inputParser;
ip.addRequired('pathForColocalization',@ischar)
ip.addParamValue('tracksNA',[],@isstruct); % selcted track ids
ip.addParamValue('movieData',[],@(x) isa(x,'MovieData')); % selcted track ids
ip.addParamValue('iChan',2,@isscalar); % This is the master channle index.
ip.addParamValue('iChanSlave',[],@isscalar); % This is the master channle index.
ip.addParamValue('labeledData',[],@iscell); % This is the master channle index.
% ip.addParamValue('outputPath','analysis1',@ischar)
ip.parse(pathForColocalization,varargin{:});
pathForColocalization=ip.Results.pathForColocalization;
tracksNA=ip.Results.tracksNA;
MD=ip.Results.movieData;
iChan=ip.Results.iChan;
iChanSlave=ip.Results.iChanSlave;
sampleFolders=ip.Results.labeledData;
% outputPath=ip.Results.outputPath;
%% Load processed data
disp('Loading raw files ...')
tic
if isempty(MD)
    coloPath = fileparts(pathForColocalization);
    MDPath = fileparts(coloPath);
    MDfilePath = [MDPath filesep 'movieData.mat'];
    MD = load(MDfilePath,'MD');
    MD = MD.MD;
end
nFrames=MD.nFrames_;


if isempty(tracksNA)
    tracksNA = load([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNA = tracksNA.tracksNA;
end
outputPath = [pathForColocalization filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
dataPath = [pathForColocalization filesep 'data'];
if ~exist(dataPath,'dir')
    mkdir(dataPath);
end
% numFrames = size(imgMap,3);
toc
%% Read other test data set
reuseSelectedGroups = 'n';
if isempty(sampleFolders)
    importSelectedGroups=input('Do you want to import existing trained data (1/0)?: ');
    if ~isempty(importSelectedGroups) && importSelectedGroups
        doneLoadingTrainedData = false;
        jj=0;
    %     T = table();
        while ~doneLoadingTrainedData
            jj=jj+1;
            % load an existing classifier
            [FileName,PathName] = uigetfile('*.mat','Select the trained data');
            curImportFilePath = fullfile(PathName,FileName);
            idGroups = load(curImportFilePath);
            idGroup1Selected = idGroups.idGroup1Selected;
            idGroup2Selected = idGroups.idGroup2Selected;
            idGroup3Selected = idGroups.idGroup3Selected;
            idGroup4Selected = idGroups.idGroup4Selected;
            idGroup5Selected = idGroups.idGroup5Selected;
            idGroup6Selected = idGroups.idGroup6Selected;
            idGroup7Selected = idGroups.idGroup7Selected;
            idGroup8Selected = idGroups.idGroup8Selected;
            idGroup9Selected = idGroups.idGroup9Selected;
            idGroupSelected{jj}={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                        idGroup7Selected,idGroup8Selected,idGroup9Selected};
            try
                curImportFilePathTracks = fullfile(PathName,'idsClassified.mat');
                curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
                curTracksNA{jj} = curTracksNAfile.tracksNA;
                mdPath = fileparts(fileparts(fileparts(fileparts(PathName))));
                curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
                curMD{jj} = curMDFile.MD;
            catch
                curImportFilePathTracks = fullfile(PathName,'idsClassified_org.mat');
                curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
                curTracksNA{jj} = curTracksNAfile.tracksNA;
                mdPath = fileparts(fileparts(fileparts(fileparts(PathName))));
                curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
                curMD{jj} = curMDFile.MD;
            end            
    %         T=[T; extractFeatureNA(curTracksNA,idGroupSelected)];
            doneLoadingTrainedData = input('Done with importing existing trained data (1/0)?: ');
        end
        nTrainingSets=jj;
        reuseSelectedGroups=input('Do you want to add some more on top of it (a), discard imported data (d) or solely use this data for classifier training(u)?: ','s');
    end
else
    nTrainingSets = numel(sampleFolders);
    for jj=1:nTrainingSets
%         curImportFilePath = fullfile(PathName,FileName);
        display(['Loading ' sampleFolders{jj} '...'])
        idGroups = load(sampleFolders{jj});
        PathName=fileparts(sampleFolders{jj});
        idGroup1Selected = idGroups.idGroup1Selected;
        idGroup2Selected = idGroups.idGroup2Selected;
        idGroup3Selected = idGroups.idGroup3Selected;
        idGroup4Selected = idGroups.idGroup4Selected;
        idGroup5Selected = idGroups.idGroup5Selected;
        idGroup6Selected = idGroups.idGroup6Selected;
        idGroup7Selected = idGroups.idGroup7Selected;
        idGroup8Selected = idGroups.idGroup8Selected;
        idGroup9Selected = idGroups.idGroup9Selected;
        idGroupSelected{jj}={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        try
            curImportFilePathTracks = fullfile(PathName,'idsClassified.mat');
            curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
            curTracksNA{jj} = curTracksNAfile.tracksNA;
            mdPath = fileparts(fileparts(fileparts(fileparts(PathName))));
            curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
            curMD{jj} = curMDFile.MD;
        catch
            try
                curImportFilePathTracks = fullfile(PathName,'idsClassified_org.mat');
                curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
                curTracksNA{jj} = curTracksNAfile.tracksNA;
            catch
                curImportFilePathTracks = fullfile(PathName,'tracksNA.mat');
                curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
                curTracksNA{jj} = curTracksNAfile.tracksNA;
            end
            mdPath = fileparts(fileparts(fileparts(fileparts(PathName))));
            try
                curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
            catch
                mdPath = fileparts(fileparts(fileparts(PathName)));
                curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
            end
            curMD{jj} = curMDFile.MD;
        end            
    end
end
%% Some initial loading and setting
try
    imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
    imgMap = imgMap.paxImgStack;
catch
    % This case SDC was not used and first img frame was used.
    paxImage=MD.getChannel(iChan).loadImage(1); 
    [h,w] = size(paxImage);

    paxImgStack = zeros(h,w,nFrames);
    for ii=1:nFrames
        paxImage=MD.getChannel(iChan).loadImage(ii); 
        paxImgStack(:,:,ii) = paxImage;
    end
    imgMap = paxImgStack;
end
nFrames=MD.nFrames_;
iFrameInterest = min(150 ,nFrames);

%% labeling target cell movie
% reuseSelectedGroups = 'n';
importSelectedGroups=true; %input('Do you want to import existing trained data for current cell (1/0)?: ');
if importSelectedGroups
    disp(['Please select selectedGroups.mat in ' pathForColocalization '.']);
    % load an existing classifier
%     [FileName,PathName] = uigetfile('*.mat','Select the trained data');
%     curImportFilePath = fullfile(PathName,FileName);
    curImportFilePath=[pathForColocalization filesep 'data/selectedGroups.mat'];
    idGroups = load(curImportFilePath);
    idGroup1SelectedCur = idGroups.idGroup1Selected;
    idGroup2SelectedCur = idGroups.idGroup2Selected;
    idGroup3SelectedCur = idGroups.idGroup3Selected;
    idGroup4SelectedCur = idGroups.idGroup4Selected;
    idGroup5SelectedCur = idGroups.idGroup5Selected;
    idGroup6SelectedCur = idGroups.idGroup6Selected;
    idGroup7SelectedCur = idGroups.idGroup7Selected;
    idGroup8SelectedCur = idGroups.idGroup8Selected;
    idGroup9SelectedCur = idGroups.idGroup9Selected;
    idGroupSelectedCur={idGroup1SelectedCur,idGroup2SelectedCur,idGroup3SelectedCur,idGroup4SelectedCur,idGroup5SelectedCur,idGroup6SelectedCur,....
                                idGroup7SelectedCur,idGroup8SelectedCur,idGroup9SelectedCur};
end
reuseSelectedGroups='u'; %input('Do you want to add some more on top of it (a), discard imported data (d) or solely use this data for classifier training(u)?: ','s');
if isempty(importSelectedGroups) || ~importSelectedGroups || strcmp(reuseSelectedGroups, 'a') || isempty(reuseSelectedGroups) || strcmp(reuseSelectedGroups, 'd')
    if strcmp(reuseSelectedGroups, 'a')
        display('Click additional tracks that belong to each group ...')
        fh2=figure;
        fh2.Color=[0 0 0];
        hold on
        colors = distinguishable_colors(9,'k');
        % switching colors between group 6 and 9
        tempColor = colors(6,:);
        colors(6,:) = colors(9,:);
        colors(9,:) = tempColor;

        for pp=1:9
            htrackG{pp}=plot(pp,1,'o','Color',colors(pp,:));
        end
        legend([htrackG{1} htrackG{2} htrackG{3} htrackG{4} htrackG{5} htrackG{6} htrackG{7} htrackG{8} htrackG{9}],...
            {'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
            'G4:retracting','G5:stable at the edge','G6:noise or very transient',...
            'G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
        legend('boxoff')
        for pp=1:9
            htrackG{pp}.Visible='off';
        end
        fh2.CurrentAxes.Color=[0 0 0];
        fh2.CurrentAxes.Visible='off';

        [idTracksAdditional, iGroupAdditional] = showAdhesionTracks(pathForColocalization,'all','tracksNA',tracksNA,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD);
%         idTracks = [idTracks idTracksAdditional];
%         iGroup = [iGroup iGroupAdditional];
        [idGroup1SelectedAdded,idGroup2SelectedAdded,idGroup3SelectedAdded,idGroup4SelectedAdded,idGroup5SelectedAdded,idGroup6SelectedAdded,...
            idGroup7SelectedAdded,idGroup8SelectedAdded,idGroup9SelectedAdded] = ...
            sortIDTracks(idTracksAdditional,iGroupAdditional);
        idGroup1SelectedAll=[idGroup1SelectedCur idGroup1SelectedAdded]; idGroup2SelectedAll=[idGroup2SelectedCur idGroup2SelectedAdded]; 
        idGroup3SelectedAll=[idGroup3SelectedCur idGroup3SelectedAdded]; idGroup4SelectedAll=[idGroup4SelectedCur idGroup4SelectedAdded]; 
        idGroup5SelectedAll=[idGroup5SelectedCur idGroup5SelectedAdded]; idGroup6SelectedAll=[idGroup6SelectedCur idGroup6SelectedAdded]; 
        idGroup7SelectedAll=[idGroup7SelectedCur idGroup7SelectedAdded]; idGroup8SelectedAll=[idGroup8SelectedCur idGroup8SelectedAdded]; 
        idGroup9SelectedAll=[idGroup9SelectedCur idGroup9SelectedAdded]; 

        idGroupSelectedAll={idGroup1SelectedAll,idGroup2SelectedAll,idGroup3SelectedAll,...
                                    idGroup4SelectedAll,idGroup5SelectedAll,idGroup6SelectedAll,...
                                    idGroup7SelectedAll,idGroup8SelectedAll,idGroup9SelectedAll};
%         TcurCell = extractFeatureNA(tracksNA,idGroupSelectedAll,kk,MD);
    else
        display('Click tracks that belong to each group ...')
    %     newTracksNA=tracksNA(~idxMatureNAs);
    %     idNAs = find(~idxMatureNAs);
        [idTracks, iGroup] = showAdhesionTracks(pathForColocalization,'all','tracksNA',tracksNA,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD);
        [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
            idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
            sortIDTracks(idTracks,iGroup);
        idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        bigEnoughGroups=cellfun(@length,idGroupSelected);
        bigEnoughGroups=find(bigEnoughGroups>=5);
        idGroupSelectedAll = idGroupSelected;
        idGroupSelectedAll(setdiff(1:9,bigEnoughGroups))={[]};
%         [TcurCell,allData]=extractFeatureNA(tracksNA,idGroupFiltered);
    end
    idGroup1SelectedTemp=idGroup1Selected; idGroup2SelectedTemp=idGroup2Selected; idGroup3SelectedTemp=idGroup3Selected; idGroup4SelectedTemp=idGroup4Selected; 
    idGroup5SelectedTemp=idGroup5Selected; idGroup6SelectedTemp=idGroup6Selected; idGroup7SelectedTemp=idGroup7Selected; idGroup8SelectedTemp=idGroup8Selected; 
    idGroup9SelectedTemp=idGroup9Selected; 

    idGroup1Selected=idGroup1SelectedAll; idGroup2Selected=idGroup2SelectedAll; idGroup3Selected=idGroup3SelectedAll; 
    idGroup4Selected=idGroup4SelectedAll; idGroup5Selected=idGroup5SelectedAll; idGroup6Selected=idGroup6SelectedAll; 
    idGroup7Selected=idGroup7SelectedAll; idGroup8Selected=idGroup8SelectedAll; idGroup9Selected=idGroup9SelectedAll; 

    save([pathForColocalization filesep 'data' filesep 'selectedGroups.mat'],'idGroup1Selected',...
        'idGroup2Selected','idGroup3Selected','idGroup4Selected','idGroup5Selected','idGroup6Selected',...
        'idGroup7Selected','idGroup8Selected','idGroup9Selected');

    idGroup1Selected=idGroup1SelectedTemp; idGroup2Selected=idGroup2SelectedTemp; idGroup3Selected=idGroup3SelectedTemp; 
    idGroup4Selected=idGroup4SelectedTemp; idGroup5Selected=idGroup5SelectedTemp; idGroup6Selected=idGroup6SelectedTemp; 
    idGroup7Selected=idGroup7SelectedTemp; idGroup8Selected=idGroup8SelectedTemp; idGroup9Selected=idGroup9SelectedTemp; 
end


%% Making classifier out of imported labeled data
C=cell(1,7); CCurCell=cell(1,7); CCross=cell(1,7); CUpdated=cell(1,7); CCrossUpdated=cell(1,7);
validationAccuracy=zeros(1,7);
validationAccuracyCurCell=zeros(1,7);
validationAccuracyCross=zeros(1,7);
validationAccuracyCrossUpdated=zeros(1,7);
for kk=1:7; % normalization method
    % Using only imported labeled data
    display('Making classifier out of imported labeled data...')
    pathEachNormalization=[pathForColocalization filesep 'normMethod' num2str(kk)];
    pathForNorm = [pathEachNormalization filesep 'classificationFromImportedLabeledData'];
    mkClrDir(pathForNorm)
    T = table();
    for jj=1:nTrainingSets
        T=[T; extractFeatureNA(curTracksNA{jj},idGroupSelected{jj},kk,curMD{jj})];
    %         T=extractFeatureNA(curTracksNA{jj},idGroupSelected{jj},ii,curMD);
    end
    %     bigEnoughGroups=cellfun(@length,idGroupSelected{jj});
    %     bigEnoughGroups=find(bigEnoughGroups>=5);
    %     idGroupFiltered = idGroupSelected{jj};
    %     idGroupFiltered(setdiff(1:9,bigEnoughGroups))={[]};

    [trainedClassifierSVM, validationAccuracySVM, CSVM, orderSVM] = trainClassifierNA(T);

    % I will use SVM no matter what, because it will be compatible with
    trainedClassifier=trainedClassifierSVM;
    validationAccuracy(kk)=validationAccuracySVM;
    C{kk}=CSVM;
    order=orderSVM;
    mkClrDir([pathForNorm filesep 'data'])
    classifierInfo = fopen([pathForNorm filesep 'data' filesep 'trainedClassifier is from SVM.txt'],'w');
    fprintf(classifierInfo, 'This is from quadratic SVM. \n');
    fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracy(kk)) '. \n']);
    disp(['Validation accuracy is ' num2str(validationAccuracy(kk)) '.'])
    fclose(classifierInfo);
    save([pathForNorm filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifier')

    % normalize confusion matrix
    for ii=1:size(C{kk},1)
        C{kk}(ii,:) = C{kk}(ii,:)/sum(C{kk}(ii,:));
    end
    response = T.Group;
    % Get the unique resonses
    totalGroups = unique(response);

    figure; confAxis=axes; imagesc(C{kk}); title('Confusion Matrix')
    set(confAxis,'xtick',1:size(C{kk},1))
    set(confAxis,'xticklabel',order')
    set(confAxis,'XTickLabelRotation',45)
    set(confAxis,'ytick',1:size(C{kk},2))
    set(confAxis,'yticklabel',order')
    xlabel('Prediction outcome')
    ylabel('Actual labels')

    c = colorbar;
    c.Label.String = 'normalized prediction';
    mkClrDir([pathForNorm filesep 'eps'])
    mkClrDir([pathForNorm filesep 'figs'])
    print('-depsc2', '-r300', [pathForNorm filesep 'eps' filesep 'confusionMatrix.eps']);
    savefig([pathForNorm filesep 'figs' filesep 'confusionMatrix.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNorm filesep 'eps' filesep 'confusionMatrix.tif'])
    close

    T = sortrows(T,size(T,2));
    features =table2array(T(:,1:end-1));
    species = table2array(T(:,end));
    nGroups = length(totalGroups);
    % normalize features
    for i = 1 : size(features,2)
        features(:,i) = (features(:,i) - min(features(:,i)))./(max(features(:,i)) - min(features(:,i)));
    end
    figure; imagesc(features');hold on
    c = colorbar;
    c.Label.String = 'feature value';
    print('-depsc2', '-r300', [pathForNorm filesep 'eps' filesep 'featureSpace.eps']);
    savefig([pathForNorm filesep 'figs' filesep 'featureSpace.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNorm filesep 'eps' filesep 'featureSpace.tif'])
    close

    D = pdist(features);
    D1 =  squareform(D);
    figure; imagesc(D1);
    title('similarityAmongTrainedData')
    c = colorbar;
    c.Label.String = 'p-dist';
    for ii=1:nGroups
    %         x0 = find(strcmp(species,['Group' num2str(ii)]),1);
    %         w = sum(strcmp(species,['Group' num2str(ii)]));
        x0 = find(strcmp(species,totalGroups{ii}),1);
        w = sum(strcmp(species,totalGroups{ii}));
        rectangle('Position',[x0-0.5 x0-0.5 w w],'EdgeColor','w','LineWidth',0.5)
    end
    print('-depsc2', '-r300', [pathForNorm filesep 'eps' filesep 'similarityAmongTrainedData.eps']);
    savefig([pathForNorm filesep 'figs' filesep 'similarityAmongTrainedData.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNorm filesep 'eps' filesep 'similarityAmongTrainedData.tif'])
    close

    Dfeats = pdist(features');
    Dfeats1 =  squareform(Dfeats);
    figure; imagesc(Dfeats1); title('similarityAmongFeatures')
    c = colorbar;
    c.Label.String = 'p-dist';

    print('-depsc2', '-r300', [pathForNorm filesep 'eps' filesep 'similarityAmongFeatures.eps']);
    savefig([pathForNorm filesep 'figs' filesep 'similarityAmongFeatures.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNorm filesep 'eps' filesep 'similarityAmongFeatures.tif'])
    close

    disp('The order is :')
    disp(order)

    %% Build classifier here
    display('Making classifier out of own labeled data...')
    pathForNormOwn = [pathEachNormalization filesep 'classificationFromOwnLabels'];
    mkClrDir(pathForNormOwn)
    TcurCell = extractFeatureNA(tracksNA,idGroupSelectedCur,kk,MD);
    [trainedClassifierCurCell, validationAccuracyCurCell(kk), CCurCell{kk}, orderCurCell] = trainClassifierNA(TcurCell);

    % I will use SVM no matter what, because it will be compatible with
    mkClrDir([pathForNormOwn filesep 'data'])
    classifierInfo = fopen([pathForNormOwn filesep 'data' filesep 'trainedClassifier is from SVM.txt'],'w');
    fprintf(classifierInfo, 'This is from quadratic SVM. \n');
    fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracyCurCell(kk)) '. \n']);
    disp(['Validation accuracy is ' num2str(validationAccuracyCurCell(kk)) '.'])
    fclose(classifierInfo);
    save([pathForNormOwn filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifierCurCell')

    % normalize confusion matrix
    for ii=1:size(CCurCell{kk},1)
        CCurCell{kk}(ii,:) = CCurCell{kk}(ii,:)/sum(CCurCell{kk}(ii,:));
    end
%     responseCurCell = TcurCell.Group;
    % Get the unique resonses
    totalGroupsCurCell = trainedClassifierCurCell.ClassNames;

    figure; confAxis=axes; imagesc(CCurCell{kk}); title('Confusion Matrix')
    set(confAxis,'xtick',1:size(CCurCell{kk},1))
    set(confAxis,'xticklabel',orderCurCell')
    set(confAxis,'XTickLabelRotation',45)
    set(confAxis,'ytick',1:size(CCurCell{kk},2))
    set(confAxis,'yticklabel',orderCurCell')
    xlabel('Prediction outcome')
    ylabel('Actual labels')
    c = colorbar;
    c.Label.String = 'normalized prediction';
    mkClrDir([pathForNormOwn filesep 'eps'])
    mkClrDir([pathForNormOwn filesep 'figs'])
    print('-depsc2', '-r300', [pathForNormOwn filesep 'eps' filesep 'confusionMatrix.eps']);
    savefig([pathForNormOwn filesep 'figs' filesep 'confusionMatrix.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormOwn filesep 'eps' filesep 'confusionMatrix.tif'])
    close

    TcurCell = sortrows(TcurCell,size(TcurCell,2));
    featuresCurCell =table2array(TcurCell(:,1:end-1));
    species = table2array(TcurCell(:,end));
    nGroupsCurCell = length(totalGroupsCurCell);
    % normalize features
    for i = 1 : size(featuresCurCell,2)
        featuresCurCell(:,i) = (featuresCurCell(:,i) - min(featuresCurCell(:,i)))./(max(featuresCurCell(:,i)) - min(featuresCurCell(:,i)));
    end
    figure; imagesc(featuresCurCell');hold on
    c = colorbar;
    c.Label.String = 'feature value';
    print('-depsc2', '-r300', [pathForNormOwn filesep 'eps' filesep 'featureSpace.eps']);
    savefig([pathForNormOwn filesep 'figs' filesep 'featureSpace.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormOwn filesep 'eps' filesep 'featureSpace.tif'])
    close

    D = pdist(featuresCurCell);
    D1 =  squareform(D);
    figure; imagesc(D1);
    title('similarityAmongTrainedData')
    c = colorbar;
    c.Label.String = 'p-dist';
    for ii=1:nGroupsCurCell
        x0 = find(strcmp(species,totalGroupsCurCell{ii}),1);
        w = sum(strcmp(species,totalGroupsCurCell{ii}));
        rectangle('Position',[x0-0.5 x0-0.5 w w],'EdgeColor','w','LineWidth',0.5)
    end
    print('-depsc2', '-r300', [pathForNormOwn filesep 'eps' filesep 'similarityAmongTrainedData.eps']);
    savefig([pathForNormOwn filesep 'figs' filesep 'similarityAmongTrainedData.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormOwn filesep 'eps' filesep 'similarityAmongTrainedData.tif'])
    close

    Dfeats = pdist(featuresCurCell');
    Dfeats1 =  squareform(Dfeats);
    figure; imagesc(Dfeats1); title('similarityAmongFeatures')
    c = colorbar;
    c.Label.String = 'p-dist';

    print('-depsc2', '-r300', [pathForNormOwn filesep 'eps' filesep 'similarityAmongFeatures.eps']);
    savefig([pathForNormOwn filesep 'figs' filesep 'similarityAmongFeatures.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormOwn filesep 'eps' filesep 'similarityAmongFeatures.tif'])
    close

    disp('The order is :')
    disp(order)
    %% Apply previous classifier to the current labeled data
    idxTinClassifier=ismember(TcurCell.Group,trainedClassifier.ClassNames);
    TcurCell= TcurCell(idxTinClassifier,:);
    [validationAccuracyCross(kk), CCross{kk}, orderCross] = validateClassifier(trainedClassifier,TcurCell);

    % Plotting confusion matrix
    % normalize confusion matrix
    for ii=1:size(CCross{kk},1)
        CCross{kk}(ii,:) = CCross{kk}(ii,:)/sum(CCross{kk}(ii,:));
    end
    curCCross=CCross{kk};
    colToRemove=find(all(isnan(curCCross), 1));
    rowToRemove=find(all(isnan(curCCross), 2));
    rowToRemove=unique([rowToRemove colToRemove]);
    curCCross(rowToRemove,:)=[];
    curCCross(:,rowToRemove)=[];
    orderCross(rowToRemove)=[];
    
    figure; confAxis=axes; imagesc(curCCross); title('Confusion Matrix')
    set(confAxis,'xtick',1:size(curCCross,1))
    set(confAxis,'xticklabel',orderCross')
    set(confAxis,'XTickLabelRotation',45)
    set(confAxis,'ytick',1:size(curCCross,2))
    set(confAxis,'yticklabel',orderCross')
    xlabel('Prediction outcome')
    ylabel('Actual labels')
    c = colorbar;
    c.Label.String = 'normalized prediction';

    pathForNormCrossCell = [pathEachNormalization filesep 'classificationFromLabelsToThisCell'];
    mkClrDir([pathForNormCrossCell filesep 'eps'])
    mkClrDir([pathForNormCrossCell filesep 'figs'])
    print('-depsc2', '-r300', [pathForNormCrossCell filesep 'eps' filesep 'confusionMatrix.eps']);
    savefig([pathForNormCrossCell filesep 'figs' filesep 'confusionMatrix.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormCrossCell filesep 'eps' filesep 'confusionMatrix.tif'])
    close

    mkClrDir([pathForNormCrossCell filesep 'data'])
    classifierInfo = fopen([pathForNormCrossCell filesep 'data' filesep 'Cross-validation accuracy.txt'],'w');
    fprintf(classifierInfo, ['Cross Validation accuracy to the cell''s labeled data is ' num2str(validationAccuracyCross(kk)) '. \n']);
    disp(['Validation accuracy is ' num2str(validationAccuracyCross(kk)) '.'])
    fclose(classifierInfo);
    
    %% Applying trainedClassifier to this cell
    [~,allData] = extractFeatureNA(tracksNA,[],kk,MD);
    allDataClass = predict(trainedClassifier,allData);

    figure; imshow(imgMap(:,:,iFrameInterest),[]), hold on
    [htrackLine, htrackCircles] = drawClassifiedTracks(allDataClass,tracksNA,iFrameInterest,[],false);
    classNames={'G1: turn-over','G2: maturing','G3: moving along protruding edge',...
        'G4: retracting','G5: stable at the edge','G6: noise or very transient','G7: adhesions at stalling edge','G8: strong stable adhesion', 'G9: weak stable adhesion inside'};
    existingClasses=~cellfun(@isempty,htrackCircles);
    % cellfun(@(x) x{1},htrackCircles(existingClasses),'UniformOutput',false)
    % this didn't work. I have to create line object array using for loop
    markerArray=[];
    for qq=find(existingClasses)
        markerArray=[markerArray htrackCircles{qq}{1}];
    end
    legend(markerArray,classNames(existingClasses),'TextColor','w','Location','best','FontSize',8,'FontWeight','bold','Box','off')
    % legend([htrackCircles{1}{1} htrackCircles{2}{1}],classNames{existingClasses},'TextColor','w','Location','best')
    print('-depsc2', '-r300', [pathForNorm filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
    savefig([pathForNorm filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNorm filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.tif'])
    close
    save([pathForNorm filesep 'data' filesep 'allDataClass.mat'],'allDataClass','-v7.3') % predicted classes
    %% Applying trainedClassifierCurCell to this cell
    allDataClassCurCell = predict(trainedClassifierCurCell,allData);

    figure; imshow(imgMap(:,:,iFrameInterest),[]), hold on
    [htrackLine, htrackCircles] = drawClassifiedTracks(allDataClassCurCell,tracksNA,iFrameInterest,[],false);
    existingClasses=~cellfun(@isempty,htrackCircles);
    % cellfun(@(x) x{1},htrackCircles(existingClasses),'UniformOutput',false)
    % this didn't work. I have to create line object array using for loop
    markerArray=[];
    for qq=find(existingClasses)
        markerArray=[markerArray htrackCircles{qq}{1}];
    end
    legend(markerArray,classNames(existingClasses),'TextColor','w','Location','best','FontSize',8,'FontWeight','bold','Box','off')
    % legend([htrackCircles{1}{1} htrackCircles{2}{1}],classNames{existingClasses},'TextColor','w','Location','best')
    print('-depsc2', '-r300', [pathForNormOwn filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
    savefig([pathForNormOwn filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormOwn filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.tif'])
    close
    save([pathForNormOwn filesep 'data' filesep 'allDataClass.mat'],'allDataClassCurCell','-v7.3') % predicted classes
    %% difference
    indSames=strcmp(allDataClassCurCell,allDataClass);
    disp(['Prediction match is: ' num2str(sum(indSames)/length(indSames)*100) '%.'])
    classifierInfo = fopen([pathForNormOwn filesep 'data' filesep 'predictionMatch.txt'],'w');
    fprintf(classifierInfo, ['Prediction match is: ' num2str(sum(indSames)/length(indSames)*100) '%.']);
    fclose(classifierInfo);
    %% Adding 5 labeled data each per group
    % idxTinClassifier=ismember(TcurCell.Group,trainedClassifier.ClassNames);
    sampleSize=5;
    uniqGroupNames=unique(TcurCell.Group);
    T_ccSample=table();
    for ii=1:numel(uniqGroupNames)
        idxThisGroupFirst=find(strcmp(TcurCell.Group,uniqGroupNames{ii}),1);
        if ismember(ii,[2,4,5,9])
            idxThisGroupNext=find(strcmp(TcurCell.Group,uniqGroupNames{ii+1}),1);
            T_ccSample=[T_ccSample;TcurCell(idxThisGroupFirst:idxThisGroupNext-1,:)];
        else
            T_ccSample=[T_ccSample;TcurCell(idxThisGroupFirst:idxThisGroupFirst+sampleSize-1,:)];
        end
    end
    Tupdated=[T; T_ccSample];
    [trainedClassifierUpdated, validationAccuracyUpdated, CUpdated{kk}, orderUpdated] = trainClassifierNA(Tupdated);

    [validationAccuracyCrossUpdated(kk), CCrossUpdated{kk},orderCrossUpdated] = validateClassifier(trainedClassifierUpdated,TcurCell);

    pathForNormUpdated = [pathEachNormalization filesep 'classificationWithPartialData'];
    mkClrDir([pathForNormUpdated filesep 'data'])
    classifierInfo = fopen([pathForNormUpdated filesep 'data' filesep 'trainedClassifier''s validation accuracy' num2str(validationAccuracyUpdated) '.txt'],'w');
    fprintf(classifierInfo, 'This is from quadratic SVM. \n');
    fprintf(classifierInfo, ['Validation accuracy of its own is ' num2str(validationAccuracyCurCell(kk)) '. \n']);
    fprintf(classifierInfo, ['Validation accuracy to the cell''s labeled data is ' num2str(validationAccuracyCrossUpdated(kk)) '. \n']);
    disp(['Validation accuracy is ' num2str(validationAccuracyCurCell(kk)) '.'])
    disp(['Cross-Validation accuracy is ' num2str(validationAccuracyCrossUpdated(kk)) '.'])
    fclose(classifierInfo);
    save([pathForNormUpdated filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifierUpdated')

    for ii=1:size(CCrossUpdated{kk},1)
        CCrossUpdated{kk}(ii,:) = CCrossUpdated{kk}(ii,:)/sum(CCrossUpdated{kk}(ii,:));
    end
    curCCrossUpdated=CCrossUpdated{kk};
    colToRemove=find(all(isnan(curCCrossUpdated), 1));
    rowToRemove=find(all(isnan(curCCrossUpdated), 2));
    rowToRemove=unique([rowToRemove colToRemove]);
    curCCrossUpdated(rowToRemove,:)=[];
    curCCrossUpdated(:,rowToRemove)=[];
    orderCrossUpdated(rowToRemove)=[];
    figure; confAxis=axes; imagesc(curCCrossUpdated); title('Confusion Matrix by updated classifer on this cell labeled data')

    set(confAxis,'xtick',1:size(curCCrossUpdated,1))
    set(confAxis,'xticklabel',orderCrossUpdated')
    set(confAxis,'XTickLabelRotation',45)
    set(confAxis,'ytick',1:size(curCCrossUpdated,2))
    set(confAxis,'yticklabel',orderCrossUpdated')
    xlabel('Prediction outcome')
    ylabel('Actual labels')

    c = colorbar;
    c.Label.String = 'normalized prediction';
    mkClrDir([pathForNormUpdated filesep 'eps'])
    mkClrDir([pathForNormUpdated filesep 'figs'])
    print('-depsc2', '-r300', [pathForNormUpdated filesep 'eps' filesep 'confusionMatrixCrossUpdated.eps']);
    savefig([pathForNormUpdated filesep 'figs' filesep 'confusionMatrixCrossUpdated.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormUpdated filesep 'eps' filesep 'confusionMatrixCrossUpdated.tif'])
    close

    for ii=1:size(CUpdated{kk},1)
        CUpdated{kk}(ii,:) = CUpdated{kk}(ii,:)/sum(CUpdated{kk}(ii,:));
    end
    figure; confAxis=axes; imagesc(CUpdated{kk}); title('Confusion Matrix by own updated labels')
    set(confAxis,'xtick',1:size(CUpdated{kk},1))
    set(confAxis,'xticklabel',orderUpdated')
    set(confAxis,'XTickLabelRotation',45)
    set(confAxis,'ytick',1:size(CUpdated{kk},2))
    set(confAxis,'yticklabel',orderUpdated')
    xlabel('Prediction outcome')
    ylabel('Actual labels')

    c = colorbar;
    c.Label.String = 'normalized prediction';
    print('-depsc2', '-r300', [pathForNormUpdated filesep 'eps' filesep 'confusionMatrixAmongUpdatedLabels.eps']);
    savefig([pathForNormUpdated filesep 'figs' filesep 'confusionMatrixAmongUpdatedLabels.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormUpdated filesep 'eps' filesep 'confusionMatrixAmongUpdatedLabels.tif'])
    close
    %% Applying for prediction again
    allDataClassUpdated = predict(trainedClassifierUpdated,allData);

    figure; imshow(imgMap(:,:,iFrameInterest),[]), hold on
    [htrackLine, htrackCircles] = drawClassifiedTracks(allDataClassUpdated,tracksNA,iFrameInterest,[],false);
    existingClasses=~cellfun(@isempty,htrackCircles);
    % cellfun(@(x) x{1},htrackCircles(existingClasses),'UniformOutput',false)
    % this didn't work. I have to create line object array using for loop
    markerArray=[];
    for qq=find(existingClasses)
        markerArray=[markerArray htrackCircles{qq}{1}];
    end
    legend(markerArray,classNames(existingClasses),'TextColor','w','Location','best','FontSize',8,'FontWeight','bold','Box','off')
    % legend([htrackCircles{1}{1} htrackCircles{2}{1}],classNames{existingClasses},'TextColor','w','Location','best')
    print('-depsc2', '-r300', [pathForNormUpdated filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
    savefig([pathForNormUpdated filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified.fig'])
    print('-dtiff', '-loose', '-r300', [pathForNormUpdated filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.tif'])
    close
    save([pathForNormUpdated filesep 'data' filesep 'allDataClass.mat'],'allDataClassCurCell','-v7.3') % predicted classes
    %% difference
    indSames=strcmp(allDataClassCurCell,allDataClassUpdated);
    disp(['Prediction match is: ' num2str(sum(indSames)/length(indSames)*100) '%.'])
    classifierInfo = fopen([pathForNormUpdated filesep 'data' filesep 'predictionMatch.txt'],'w');
    fprintf(classifierInfo, ['Prediction match is: ' num2str(sum(indSames)/length(indSames)*100) '%.']);
    fclose(classifierInfo);
end
end
function  [validationAccuracy,C,order] = validateClassifier(trainedClassifier,datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs',...
    'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs',...
    'maxIntensityNAs', 'timeToMaxInten', 'edgeVariation'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Perform cross-validation
% figure; imagesc(predictors);
% predictedLabels = predict(trainedClassifier,predictors);
% [predictedLabels,NegLoss,PBScore] = trainedClassifier.predict(predictors);
predictedLabels = trainedClassifier.predict(predictors);
% if the number of response is different (e.g. no Group6), merge into
% Group9.
results = nan(1,numel(predictedLabels));
for i = 1 : numel(predictedLabels)
    results(i) = strcmp(predictedLabels{i},response{i});
end
validationAccuracy=sum(results)/length(results);
% confusion matrix
ascendingOrder = unique(sort([response; predictedLabels]));
[C,order] = confusionmat(response,predictedLabels,'order',ascendingOrder);
end
