function []=testFeatureNormalizationForClassification(pathForColocalization,varargin)
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
ip.parse(pathForColocalization,varargin{:});
pathForColocalization=ip.Results.pathForColocalization;
tracksNA=ip.Results.tracksNA;
MD=ip.Results.movieData;
iChan=ip.Results.iChan;
iChanSlave=ip.Results.iChanSlave;
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
% try
%     tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
%     tMap = tMap.tMap;
% catch
%     TFMpackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
%     forceProc =TFMpackage.processes_{4};
% %     forceField = forceProc.loadChannelOutput;
%     tMap = load(forceProc.outFilePaths_{2});
%     tMap = tMap.tMap;
% end
outputPath = [pathForColocalization filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
dataPath = [pathForColocalization filesep 'data'];
if ~exist(dataPath,'dir')
    mkdir(dataPath);
end
numFrames = size(imgMap,3);
startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
% movieData to find out pixel size
pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron
toc
%% Read other test data set
reuseSelectedGroups = 'n';
importSelectedGroups=input('Do you want to import existing trained data (1/0)?: ');
if ~isempty(importSelectedGroups) && importSelectedGroups
    doneLoadingTrainedData = false;
    T = table();
    while ~doneLoadingTrainedData
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
        idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        try
            curImportFilePathTracks = fullfile(PathName,'idsClassified.mat');
            curTracksNA = load(curImportFilePathTracks,'tracksNA');
            curTracksNA = curTracksNA.tracksNA;
        catch
            curImportFilePathTracks = fullfile(PathName,'idsClassified_org.mat');
            curTracksNA = load(curImportFilePathTracks,'tracksNA');
            curTracksNA = curTracksNA.tracksNA;
        end            
%         T=[T; extractFeatureNA(curTracksNA,idGroupSelected)];
        doneLoadingTrainedData = input('Done with importing existing trained data (1/0)?: ');
    end
    reuseSelectedGroups=input('Do you want to add some more on top of it (a), discard imported data (d) or solely use this data for classifier training(u)?: ','s');
end
%% New training data labeling - for third step - select 'u' in previous question
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
        [idTracksAdditional, iGroupAdditional] = showAdhesionTracks(pathForColocalization,'all','tracksNA',tracksNA,'trainedData',T,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD);
%         idTracks = [idTracks idTracksAdditional];
%         iGroup = [iGroup iGroupAdditional];
        [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
            idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
            sortIDTracks(idTracksAdditional,iGroupAdditional);
        idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        [curT,allData,meas] = extractFeatureNA(tracksNA,idGroupSelected);
        T=[T; curT];
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
        idGroupFiltered = idGroupSelected;
        idGroupFiltered(setdiff(1:9,bigEnoughGroups))={[]};
        [T,allData]=extractFeatureNA(tracksNA,idGroupFiltered);
    end
    save([pathForColocalization filesep 'data' filesep 'selectedGroups.mat'],'idGroup1Selected',...
        'idGroup2Selected','idGroup3Selected','idGroup4Selected','idGroup5Selected','idGroup6Selected',...
        'idGroup7Selected','idGroup8Selected','idGroup9Selected');
end

%% feature extraction
if importSelectedGroups && strcmp(reuseSelectedGroups, 'u')
    bigEnoughGroups=cellfun(@length,idGroupSelected);
    bigEnoughGroups=find(bigEnoughGroups>=5);
    idGroupFiltered = idGroupSelected;
    idGroupFiltered(setdiff(1:9,bigEnoughGroups))={[]};
    [T,allData]=extractFeatureNA(tracksNA,idGroupFiltered);
else
    T=extractFeatureNA(curTracksNA,idGroupSelected);
end
[trainedClassifierSVM, validationAccuracySVM, CSVM, orderSVM] = trainClassifierNA(T);
% I will use SVM no matter what, because it will be compatible with
trainedClassifier=trainedClassifierSVM;
validationAccuracy=validationAccuracySVM;
C=CSVM;
order=orderSVM;
classifierInfo = fopen([pathForColocalization filesep 'data' filesep 'trainedClassifier is from SVM.txt'],'w');
fprintf(classifierInfo, 'This is from quadratic SVM. \n');
fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracy) '. \n']);
disp(['Validation accuracy is ' num2str(validationAccuracy) '.'])
fclose(classifierInfo);
save([pathForColocalization filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifier')

% normalize confusion matrix
for ii=1:size(C,1)
    C(ii,:) = C(ii,:)/sum(C(ii,:));
end
response = T.Group;
% Get the unique resonses
totalGroups = unique(response);

figure; confAxis=axes; imagesc(C); title('Confusion Matrix')
set(confAxis,'xticklabel',totalGroups')
set(confAxis,'yticklabel',totalGroups')
c = colorbar;
c.Label.String = 'normalized prediction';
print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'confusionMatrix.eps']);
savefig([pathForColocalization filesep 'figs' filesep 'confusionMatrix.fig'])
print('-dtiff', '-loose', '-r300', [pathForColocalization filesep 'eps' filesep 'confusionMatrix.tif'])

T = sortrows(T,13);
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
print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'featureSpace.eps']);
savefig([pathForColocalization filesep 'figs' filesep 'featureSpace.fig'])
print('-dtiff', '-loose', '-r300', [pathForColocalization filesep 'eps' filesep 'featureSpace.tif'])

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
print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongTrainedData.eps']);
savefig([pathForColocalization filesep 'figs' filesep 'similarityAmongTrainedData.fig'])
print('-dtiff', '-loose', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongTrainedData.tif'])

Dfeats = pdist(features');
Dfeats1 =  squareform(Dfeats);
figure; imagesc(Dfeats1); title('similarityAmongFeatures')
c = colorbar;
c.Label.String = 'p-dist';

print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongFeatures.eps']);
savefig([pathForColocalization filesep 'figs' filesep 'similarityAmongFeatures.fig'])
print('-dtiff', '-loose', '-r300', [pathForColocalization filesep 'eps' filesep 'similarityAmongFeatures.tif'])

disp('The order is :')
disp(order)

    
    
allDataClass = predict(trainedClassifier,allData);

figure, imshow(imgMap(:,:,end),[])
hold on
colors = distinguishable_colors(9,'k');
idGroup1 = strcmp(allDataClass,'Group1');
idGroup2 = strcmp(allDataClass,'Group2');
idGroup3 = strcmp(allDataClass,'Group3');
idGroup4 = strcmp(allDataClass,'Group4');
idGroup5 = strcmp(allDataClass,'Group5');
idGroup6 = strcmp(allDataClass,'Group6');
idGroup7 = strcmp(allDataClass,'Group7');
idGroup8 = strcmp(allDataClass,'Group8');
idGroup9 = strcmp(allDataClass,'Group9');
htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1));
htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup2),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup2));
htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup3),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup3));
htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(4,:)),tracksNA(idGroup4),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(4,:)),tracksNA(idGroup4));
htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(5,:)),tracksNA(idGroup5),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(5,:)),tracksNA(idGroup5));
htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(9,:)),tracksNA(idGroup6),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(9,:)),tracksNA(idGroup6));
htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(7,:)),tracksNA(idGroup7),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(7,:)),tracksNA(idGroup7));
htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(8,:)),tracksNA(idGroup8),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(8,:)),tracksNA(idGroup8));
htrackG9=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(6,:)),tracksNA(idGroup9),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(6,:)),tracksNA(idGroup9));
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1} htrackG9{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
% legend('boxoff')
print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
savefig([pathForColocalization filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified.fig'])
print('-dtiff', '-loose', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.tif'])
save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','tracksNA','-v7.3')
end
function  [validationAccuracy,C,order] = validateClassifier(trainedClassifier,datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Perform cross-validation
% figure; imagesc(predictors);
% predictedLabels = predict(trainedClassifier,predictors);
% [predictedLabels,NegLoss,PBScore] = trainedClassifier.predict(predictors);
predictedLabels = trainedClassifier.predict(predictors);
results = nan(1,numel(predictedLabels));
for i = 1 : numel(predictedLabels)
    results(i) = strcmp(predictedLabels{i},response{i});
end
validationAccuracy=sum(results)/length(results);
% confusion matrix
[C,order] = confusionmat(response,predictedLabels);
end
