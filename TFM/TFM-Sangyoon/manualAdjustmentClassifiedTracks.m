function [] = manualAdjustmentClassifiedTracks( MD )
%function [] = manualLabelingAfterInitialClassification( MD ) opens an
% interactive window with classified adhesion tracks so that the user can
% select several adhesions to change their assignment to a right class. 

% Show some instruction window to tell a user what to do

%% Let the user know whats going on
waitHan = msgbox({'You will see the cell window and classified adhesion tracks with color';
'label. After closing this window, please select colored (classified)';
'adhesions to correct their assignment to add the adjusted class';
'assignment to the idsClassified.mat, which will be updated after this process.'});

uiwait(waitHan);            

%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
% p = classProc.funParams_;

%% Load the selectedGroups.mat
iChan = find(classProc.checkChannelOutput);
idsClassifiedStruct = load(classProc.outFilePaths_{4,iChan});

%% Load tracksNA
adhAnalProc = faPackage.getProcess(7);
tracksNA=adhAnalProc.loadChannelOutput(iChan,'output','tracksNA');

%% Load imgStack, forceStack and anyother stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
% Launch pickAdhesion window with labeled adhesions with a right color and
% unlabed ones with white color. Get the right classes per newly selected
% adhesions
[IDs, iGroups, iPrevGroups,tracksNA]=pickAdhesionTracksInteractive(tracksNA, imgStack,...
    'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2, 'idSelected',idsClassifiedStruct);


%% Combine the indices into a new selectedGroups (take out ids from
% idGroupSelectedStruct if changed in iGroups, and put into a new
% idGroupSelectedStruct).
% % check if the included iGroups were changing from iPrevGroups
% iChangedGroups = iGroups ~= iPrevGroups;
idGroupSelectedStruct = classProc.loadChannelOutput(iChan, 'output', 'selectedGroups');
numGroups = 9;
idGroupSelectedCell=cell(1,numGroups);
% go over each label and assign and de-assign the group
for kk=1:numel(IDs)
    curID = IDs(kk);
    curPrevGroup = iPrevGroups(kk);
    curGroup = iGroups(kk);
    % Deassign if it was belong to another group
    groupPrevName = ['idGroup' num2str(curPrevGroup) 'Selected'];
    if ismember(curID,idGroupSelectedStruct.(groupPrevName))
        idGroupSelectedStruct.(groupPrevName)(ismember(idGroupSelectedStruct.(groupPrevName),curID))=[];
    end
    % Assign to the current group
    groupCurName = ['idGroup' num2str(curGroup) 'Selected'];
    if ~ismember(curID,idGroupSelectedStruct.(groupCurName))
        idGroupSelectedStruct.(groupCurName)=[idGroupSelectedStruct.(groupCurName) curID];
    end
end
for curGroup=1:numGroups
    groupCurName = ['idGroup' num2str(curGroup) 'Selected'];
    idGroupSelectedCell{curGroup}=idGroupSelectedStruct.(groupCurName);
end
%% Perform classification for entire population of tracksNA
% including previous training data
p=classProc.funParams_;
sampleFolders=p.labelData;

nTrainingSets = numel(sampleFolders);
if ~isempty(sampleFolders)
    for jj=1:nTrainingSets
%         curImportFilePath = fullfile(PathName,FileName);
        disp(['Loading ' sampleFolders{jj} '...'])
        [idGroupSelectedAll{jj},MDAll{jj},curTracksNA{jj}]=loadSampleFolder(sampleFolders{jj});
    end
end

% Current one
nTrainingSets=nTrainingSets+1;
curTracksNA{nTrainingSets}=tracksNA;
idGroupSelectedAll{nTrainingSets}=idGroupSelectedCell;
MDAll{nTrainingSets}=MD;


T = table();kk=2;
for jj=1:nTrainingSets
    T=[T; extractFeatureNA(curTracksNA{jj},idGroupSelectedAll{jj},kk,MDAll{jj})];
%         T=extractFeatureNA(curTracksNA{jj},idGroupSelected{jj},ii,curMD);
end

disp(['Training SVM classifier (' num2str(size(T,1)) ' labeled data)...'])
tic
[trainedClassifier, validationAccuracySVM, CSVM, orderSVM] = trainClassifierNA(T);
toc
validationAccuracy=validationAccuracySVM;
C=CSVM;
order=orderSVM;
classifierInfo = fopen([p.OutputDirectory filesep 'data' filesep 'trainedClassifier is from SVM.txt'],'w');
fprintf(classifierInfo, 'This is from quadratic SVM. \n');
fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracy) '. \n']);

%% normalize confusion matrix
for ii=1:size(C,1)
    C(ii,:) = C(ii,:)/sum(C(ii,:));
end
response = T.Group;
% Get the unique resonses
totalGroups = unique(response);

confFig=figure; confAxis=axes; imagesc(C); title('Confusion Matrix')
set(confAxis,'xtick',1:size(C,1))
set(confAxis,'xticklabel',order')
set(confAxis,'XTickLabelRotation',45)
set(confAxis,'ytick',1:size(C,2))
set(confAxis,'yticklabel',order')
xlabel('Prediction outcome')
ylabel('Actual labels')

c = colorbar;
c.Label.String = 'normalized prediction';
print(confFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'confusionMatrix.eps']);
savefig(confFig,[p.OutputDirectory filesep 'figs' filesep 'confusionMatrix.fig'])
print(confFig,'-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'confusionMatrix.tif'])
%% Update idClassified.mat
[~,allData] = extractFeatureNA(tracksNA,[],2,MD);
allDataClass = predict(trainedClassifier,allData);
idGroup1 = strcmp(allDataClass,'Group1');
idGroup2 = strcmp(allDataClass,'Group2');
idGroup3 = strcmp(allDataClass,'Group3');
idGroup4 = strcmp(allDataClass,'Group4');
idGroup5 = strcmp(allDataClass,'Group5');
idGroup6 = strcmp(allDataClass,'Group6');
idGroup7 = strcmp(allDataClass,'Group7');
idGroup8 = strcmp(allDataClass,'Group8');
idGroup9 = strcmp(allDataClass,'Group9');

% minority redemption
indexG1 = idGroupSelectedStruct.idGroup1Selected; 
indexG2 = idGroupSelectedStruct.idGroup2Selected; 
indexG3 = idGroupSelectedStruct.idGroup3Selected; 

if sum(idGroup1)<numel(indexG1)
    idGroup1(indexG1)=true;
end
if sum(idGroup2)<numel(indexG2)
    idGroup2(indexG2)=true;
end
if sum(idGroup3)<numel(indexG3)
    idGroup3(indexG3)=true;
end

save(classProc.outFilePaths_{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3')

disp('Classification Adjustment Done!')
end

