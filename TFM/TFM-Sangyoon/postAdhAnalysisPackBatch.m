%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
if ~ischar(pathSFolders) && pathSFolders==0
    analysisFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~analysisFolderSelectionDone
        ii=ii+1;
%         curPathProject = uigetdir(rootFolder,'Select each analysis folder that contains movieList.mat (Click Cancel when no more)');
        [curMLFile,curPathProject] = uigetfile(rootFolder,'Select the movie list file one per each attempt (Click Cancel when no more)');
        if ~ischar(curPathProject) && curPathProject==0
            analysisFolderSelectionDone=true;
        else
            [~,finalFolder] = fileparts(curPathProject);
            pathAnalysisAll{ii} = curPathProject;
            groupNames{ii} = finalFolder;
            MLFileNamesAll{ii} = curMLFile;
            MLNames{ii} = 'movieList.mat';
        end
    end
    if analysisFolderSelectionDone && ii==1
        MLSelectionDone = false;
        ii=0;
        while ~MLSelectionDone
            ii=ii+1;
            [nameML, curPathML] = uigetfile('*.mat','Select each movieList (Click Cancel when no more)');
            if ~ischar(curPathML) || isempty(curPathML)
                MLSelectionDone=true;
            else
                curML = load(fullfile(curPathML,nameML),'ML'); curML=curML.ML;
                pathAnalysisAll{ii} = curML.getPath;
                try
                    groupNames{ii} = nameML(10:end-4); %Excluding first 'movieList' and last '.mat'
                catch % when movieList is just movieList.mat, use the name of the containing folder
                    [~,finalFolder] = fileparts(pathAnalysisAll{ii});
                    groupNames{ii} = finalFolder;
                end
                MLNames{ii} = nameML;
            end
        end
    end
    specificName = strjoin(groupNames);
    rootAnalysis = pathAnalysisAll{1};
    save([rootAnalysis filesep 'selectedFolders' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
else
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
    specificName=fileSFolders(16:end);
    for k=1:numel(pathAnalysisAll)
        MLFileNamesAll{k} = 'movieList.mat';
    end
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);

for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLFileNamesAll{k}]);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummary_Adhesion' specificName '/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary_Adhesion' specificName '/Data'];
mkdir(dataPath)
%% Collecting general adhesion-reated features
numClasses=9;
N=zeros(numConditions,1);
NAdensityGroup = cell(numConditions,1);
FAareaGroup = cell(numConditions,1);
FAlengthGroup = cell(numConditions,1);
FAdensityGroup = cell(numConditions,1);
NAstructGroup= cell(numConditions,1);
FAstructGroup= cell(numConditions,1);
initRiseGroup = cell(numConditions,1);
peakGroup = cell(numConditions,1);
endTimeGroup = cell(numConditions,1);

numPureFAsGroup = cell(numConditions,1);
cellAreaGroup = cell(numConditions,1);
NADensityGroup = cell(numConditions,1);
FADensityGroup = cell(numConditions,1);
lifeTimeAllGroup = cell(numConditions,1);
lifeTimeFailingNAsGroup = cell(numConditions,1);
lifeTimeMaturingNAsGroup = cell(numConditions,1);
maturingRatioGroup = cell(numConditions,1);
assemRateGroup = cell(numConditions,1);
disassemRateGroup = cell(numConditions,1);
nucleatingNARatioGroup = cell(numConditions,1);
disassemNARatioGroup = cell(numConditions,1);

numG1Group = cell(numConditions,1);
numG2Group = cell(numConditions,1);
numG3Group = cell(numConditions,1);
numG4Group = cell(numConditions,1);
numG5Group = cell(numConditions,1);
numG6Group = cell(numConditions,1);
numG7Group = cell(numConditions,1);
numG8Group = cell(numConditions,1);
numG9Group = cell(numConditions,1);

meanRelativePopG1Group = cell(numConditions,1);
meanRelativePopG2Group = cell(numConditions,1);
meanRelativePopG3Group = cell(numConditions,1);
meanRelativePopG4Group = cell(numConditions,1);
meanRelativePopG5Group = cell(numConditions,1);
meanRelativePopG6Group = cell(numConditions,1);
meanRelativePopG7Group = cell(numConditions,1);
meanRelativePopG8Group = cell(numConditions,1);
meanRelativePopG9Group = cell(numConditions,1);

meanAdhDensityG1Group = cell(numConditions,1);
meanAdhDensityG2Group = cell(numConditions,1);
meanAdhDensityG3Group = cell(numConditions,1);
meanAdhDensityG4Group = cell(numConditions,1);
meanAdhDensityG5Group = cell(numConditions,1);
meanAdhDensityG6Group = cell(numConditions,1);
meanAdhDensityG7Group = cell(numConditions,1);
meanAdhDensityG8Group = cell(numConditions,1);
meanAdhDensityG9Group = cell(numConditions,1);

iAdhChan = input('Your adhesion channel of interest? (default: 2): ');
if isempty(iAdhChan); iAdhChan=2; end

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);
    curNAdensity = zeros(N(ii),1);
    curFAarea = cell(N(ii),1);
    curFAlength = cell(N(ii),1);
    curFAdensity = zeros(N(ii),1);
    curML = MLAll(ii);

    tracksAdhs=cell(N(ii),1);
    indFail=cell(N(ii),1);
    indMature=cell(N(ii),1);
    lifeTimeNAfailing=cell(N(ii),1);
    lifeTimeNAmaturing=cell(N(ii),1);
    maturingRatio=zeros(N(ii),1);
    NADensity=cell(N(ii),1);
    FADensity=cell(N(ii),1);
    numNAsInBand=cell(N(ii),1);
    focalAdhInfo=cell(N(ii),1);    
    assemRateCell=cell(N(ii),1);    
    disassemRateCell=cell(N(ii),1);    
    nucleationRatio=cell(N(ii),1);    
    disassemblingNARatio=cell(N(ii),1);   
    
    meanMedianFAarea=zeros(N(ii),1);
    meanMedianFAlength=zeros(N(ii),1);
    meanNumPureFAs=zeros(N(ii),1);
    meanCellArea=zeros(N(ii),1);
    meanNADensity=zeros(N(ii),1);
    meanFADensity=zeros(N(ii),1);
    meanLifeTimeAll=zeros(N(ii),1);
    meanLifeTimeFailingNAs=zeros(N(ii),1);
    meanLifeTimeMaturingNAs=zeros(N(ii),1);
    meanMaturingRatio=zeros(N(ii),1);
    meanAssemRate=zeros(N(ii),1);
    meanDisassemRate=zeros(N(ii),1);
    meanNucleatingNARatio=zeros(N(ii),1);
    meanDisassemNARatio=zeros(N(ii),1);
    
    numG1=zeros(N(ii),1);
    numG2=zeros(N(ii),1);
    numG3=zeros(N(ii),1);
    numG4=zeros(N(ii),1);
    numG5=zeros(N(ii),1);
    numG6=zeros(N(ii),1);
    numG7=zeros(N(ii),1);
    numG8=zeros(N(ii),1);
    numG9=zeros(N(ii),1);
    
    meanRelativePopG1=zeros(N(ii),1);
    meanRelativePopG2=zeros(N(ii),1);
    meanRelativePopG3=zeros(N(ii),1);
    meanRelativePopG4=zeros(N(ii),1);
    meanRelativePopG5=zeros(N(ii),1);
    meanRelativePopG6=zeros(N(ii),1);
    meanRelativePopG7=zeros(N(ii),1);
    meanRelativePopG8=zeros(N(ii),1);
    meanRelativePopG9=zeros(N(ii),1);
    
    meanAdhDensityG1=zeros(N(ii),1);
    meanAdhDensityG2=zeros(N(ii),1);
    meanAdhDensityG3=zeros(N(ii),1);
    meanAdhDensityG4=zeros(N(ii),1);
    meanAdhDensityG5=zeros(N(ii),1);
    meanAdhDensityG6=zeros(N(ii),1);
    meanAdhDensityG7=zeros(N(ii),1);
    meanAdhDensityG8=zeros(N(ii),1);
    meanAdhDensityG9=zeros(N(ii),1);
    
    curMovies = curML.movies_;
    for k=1:N(ii)
        disp(['ii: ' num2str(ii) '  k:' num2str(k)])
        % get the movie
        curMovie=curMovies{k};
        % get the general adhesion-reated features
        iCurPack = curMovie.getPackageIndex('FocalAdhesionPackage');
        curFAPackage = curMovie.getPackage(iCurPack);
        % get the general adhesion-reated features
        iCurAnalyProc = 7;
        curAnalProc = curFAPackage.getProcess(iCurAnalyProc);
        
%     outFilename = [chanDirName '_Chan' num2str(i) '_focalAdhInfo'];
%     outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_NAFADensity'];
%     outFilePaths{3,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_allAnalysisFA'];
%     outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_assemRates'];
%     outFilePaths{5,i} = [p.OutputDirectory filesep outFilename '.mat'];
        
        faInfoStruct=load(curAnalProc.outFilePaths_{2,iAdhChan});
        curFocalAdhInfo = faInfoStruct.focalAdhInfo;
        curMedianFAareaAllFrames = arrayfun(@(x) x.medianFAarea,curFocalAdhInfo);
        meanMedianFAarea(k) = mean(curMedianFAareaAllFrames);
        curMedianFAlengthAllFrames = arrayfun(@(x) x.medianLength,curFocalAdhInfo);
        meanMedianFAlength(k) = mean(curMedianFAlengthAllFrames);
        curNumPureFAs = arrayfun(@(x) x.numberPureFA,curFocalAdhInfo);
        meanNumPureFAs(k) = mean(curNumPureFAs);
        curCellAreas = arrayfun(@(x) x.cellArea,curFocalAdhInfo); % in um2
        meanCellArea(k) = mean(curCellAreas);
        
        nafaStruct=load(curAnalProc.outFilePaths_{3,iAdhChan});
        curNADensity = nafaStruct.NADensity;
        meanNADensity(k) = mean(curNADensity);
        curFADensity = nafaStruct.FADensity;
        meanFADensity(k) = meanNumPureFAs(k)/meanCellArea(k); % in num/um2

        maturingStruct=load(curAnalProc.outFilePaths_{4,iAdhChan});
        meanLifeTimeAll(k)=mean(maturingStruct.lifeTimeAll);
        meanLifeTimeFailingNAs(k)=mean(maturingStruct.lifeTimeNAfailing);
        meanLifeTimeMaturingNAs(k)=mean(maturingStruct.lifeTimeNAmaturing);
        meanMaturingRatio(k)=mean(maturingStruct.maturingRatio);
        
        assemRateStruct=load(curAnalProc.outFilePaths_{5,iAdhChan});
        meanAssemRate(k)=nanmean(assemRateStruct.assemRateCell);
        meanDisassemRate(k)=nanmean(assemRateStruct.disassemRateCell);
        meanNucleatingNARatio(k) = mean(assemRateStruct.nucleationRatio);
        meanDisassemNARatio(k) = mean(assemRateStruct.disassemblingNARatio);
        
        % get the classification-reated features
        iClassiProc = 8;
        classProc = curFAPackage.getProcess(iClassiProc);
%     outFilename = [chanDirName '_Chan' num2str(i) '_idsClassified'];
%     outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '.mat'];
        idsStruct=load(classProc.outFilePaths_{4,iAdhChan});
        numG1(k) = sum(idsStruct.idGroup1);
        numG2(k) = sum(idsStruct.idGroup2);
        numG3(k) = sum(idsStruct.idGroup3);
        numG7(k) = sum(idsStruct.idGroup7);
        numG9(k) = sum(idsStruct.idGroup9);
        numG5(k) = sum(idsStruct.idGroup5);
        numG8(k) = sum(idsStruct.idGroup8);
        numG4(k) = sum(idsStruct.idGroup4);
        
        classiOutFolder = fileparts(classProc.outFilePaths_{4,iAdhChan});
        classiDataFolder = [classiOutFolder filesep 'data'];
        % Mean relative population for each group
        meanRelativePopG1(k) = sum(idsStruct.idGroup1)/length(idsStruct.idGroup1);
        meanRelativePopG2(k) = sum(idsStruct.idGroup2)/length(idsStruct.idGroup1);
        meanRelativePopG3(k) = sum(idsStruct.idGroup3)/length(idsStruct.idGroup1);
        meanRelativePopG4(k) = sum(idsStruct.idGroup4)/length(idsStruct.idGroup1);
        meanRelativePopG5(k) = sum(idsStruct.idGroup5)/length(idsStruct.idGroup1);
        meanRelativePopG6(k) = sum(idsStruct.idGroup6)/length(idsStruct.idGroup1);
        meanRelativePopG7(k) = sum(idsStruct.idGroup7)/length(idsStruct.idGroup1);
        meanRelativePopG8(k) = sum(idsStruct.idGroup8)/length(idsStruct.idGroup1);
        meanRelativePopG9(k) = sum(idsStruct.idGroup9)/length(idsStruct.idGroup1);
        % Mean adh density per each group
        meanAdhDensityG1(k) = sum(idsStruct.idGroup1)/meanCellArea(k);
        meanAdhDensityG2(k) = sum(idsStruct.idGroup2)/meanCellArea(k);
        meanAdhDensityG3(k) = sum(idsStruct.idGroup3)/meanCellArea(k);
        meanAdhDensityG4(k) = sum(idsStruct.idGroup4)/meanCellArea(k);
        meanAdhDensityG5(k) = sum(idsStruct.idGroup5)/meanCellArea(k);
        meanAdhDensityG6(k) = sum(idsStruct.idGroup6)/meanCellArea(k);
        meanAdhDensityG7(k) = sum(idsStruct.idGroup7)/meanCellArea(k);
        meanAdhDensityG8(k) = sum(idsStruct.idGroup8)/meanCellArea(k);
        meanAdhDensityG9(k) = sum(idsStruct.idGroup9)/meanCellArea(k);
        % Other feature related properties are calculated in the step 11
        
    end
    
    
    initRiseAgainstForceEachClass=cell(numClasses,1);
    peakTimeAgainstForceEachClass=cell(numClasses,1);
    endTimeAgainstForceEachClass=cell(numClasses,1);
    iForceSlave = 1;
    for pp=1:numClasses
        for k=1:N(ii)
            curMovie=curMovies{k};
            iCurPack = curMovie.getPackageIndex('FocalAdhesionPackage');
            curFAPackage = curMovie.getPackage(iCurPack);
            % Initial rise
            iInitRiseProc = 11;
            initRiseProc = curFAPackage.getProcess(iInitRiseProc);
            % I decided to get the adjusted time lags
            initOutFolder = fileparts(initRiseProc.outFilePaths_{2,iForceSlave});
            initDataPath = [initOutFolder filesep 'data'];
            
            nameTitle=['initialLag Class' num2str(pp)];
            initRiseStruct = load([initDataPath filesep nameTitle],'initialLagTogetherAdjusted','nameList2');   
            curInitTimeLag=initRiseStruct.initialLagTogetherAdjusted;
            initRiseAgainstForceEachClass{pp}{k} = curInitTimeLag;
            
            nameTitle=['peakLag Class' num2str(pp)];
            peakStruct = load([initDataPath filesep nameTitle],'peakLagTogetherAdjusted','nameList2');   
            curPeakTimeLag=peakStruct.peakLagTogetherAdjusted;
            peakTimeAgainstForceEachClass{pp}{k} = curPeakTimeLag;
            
            nameTitle=['endingLag Class' num2str(pp)];
            endTimeStruct = load([initDataPath filesep nameTitle],'endingLagTogetherAdjusted','nameList2');   
            curEndTimeLag=endTimeStruct.endingLagTogetherAdjusted;
            endTimeAgainstForceEachClass{pp}{k} = curEndTimeLag;
        end
    end
    FAareaGroup{ii,1}=meanMedianFAarea;
    FAlengthGroup{ii,1}=meanMedianFAlength;
    numPureFAsGroup{ii,1}=meanNumPureFAs;
    cellAreaGroup{ii,1}=meanCellArea;
    NADensityGroup{ii,1}=meanNADensity;
    FADensityGroup{ii,1}=meanFADensity;
    lifeTimeAllGroup{ii,1}=meanLifeTimeAll;
    lifeTimeFailingNAsGroup{ii,1}=meanLifeTimeFailingNAs;
    lifeTimeMaturingNAsGroup{ii,1}=meanLifeTimeMaturingNAs;
    maturingRatioGroup{ii,1}=meanMaturingRatio;
    assemRateGroup{ii,1}=meanAssemRate;
    disassemRateGroup{ii,1}=meanDisassemRate;
    nucleatingNARatioGroup{ii,1}=meanNucleatingNARatio;
    disassemNARatioGroup{ii,1}=meanDisassemNARatio;
    
    numG1Group{ii,1}=numG1;
    numG2Group{ii,1}=numG2;
    numG3Group{ii,1}=numG3;
    numG4Group{ii,1}=numG4;
    numG5Group{ii,1}=numG5;
    numG6Group{ii,1}=numG6;
    numG7Group{ii,1}=numG7;
    numG8Group{ii,1}=numG8;
    numG9Group{ii,1}=numG9;

    meanRelativePopG1Group{ii,1}=meanRelativePopG1;
    meanRelativePopG2Group{ii,1}=meanRelativePopG2;
    meanRelativePopG3Group{ii,1}=meanRelativePopG3;
    meanRelativePopG4Group{ii,1}=meanRelativePopG4;
    meanRelativePopG5Group{ii,1}=meanRelativePopG5;
    meanRelativePopG6Group{ii,1}=meanRelativePopG6;
    meanRelativePopG7Group{ii,1}=meanRelativePopG7;
    meanRelativePopG8Group{ii,1}=meanRelativePopG8;
    meanRelativePopG9Group{ii,1}=meanRelativePopG9;

    meanAdhDensityG1Group{ii,1}=meanAdhDensityG1;
    meanAdhDensityG2Group{ii,1}=meanAdhDensityG2;
    meanAdhDensityG3Group{ii,1}=meanAdhDensityG3;
    meanAdhDensityG4Group{ii,1}=meanAdhDensityG4;
    meanAdhDensityG5Group{ii,1}=meanAdhDensityG5;
    meanAdhDensityG6Group{ii,1}=meanAdhDensityG6;
    meanAdhDensityG7Group{ii,1}=meanAdhDensityG7;
    meanAdhDensityG8Group{ii,1}=meanAdhDensityG8;
    meanAdhDensityG9Group{ii,1}=meanAdhDensityG9;
    
    initRiseGroup{ii,1}=initRiseAgainstForceEachClass;
    peakGroup{ii,1}=peakTimeAgainstForceEachClass;
    endTimeGroup{ii,1}=endTimeAgainstForceEachClass;
end
disp('Done')
%% setting up group name
<<<<<<< HEAD
for ii=1:numConditions
    [pathFolder, finalFolder]=fileparts(pathAnalysisAll{ii});
    if isempty(finalFolder)
        [~, finalFolder]=fileparts(pathFolder);
    end
    groupNames{ii} = finalFolder;
end
=======
% for ii=1:numConditions
%     [~, finalFolder]=fileparts(pathAnalysisAll{ii});
%     groupNames{ii} = finalFolder;
% end
>>>>>>> 1f4ab8ef40f405b766c09b86da19b0e5da324b43
%% Plotting each - initRiseGroup - all classes
nameList=groupNames'; %{'pLVX' 'P29S'};
for curGroup=1:9
    initRiseGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),initRiseGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(initRiseGroupEach,nameList,1,false,true);
    ylabel(['Initial rise in group ' num2str(curGroup) ' (sec)'])
    title(['Initial rise in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],'-dtiff')

    tableInitRiseGroupEach=table(initRiseGroupEach,'RowNames',nameList);
    writetable(tableInitRiseGroupEach,[dataPath filesep 'initialRizeG' num2str(curGroup) '.csv'])
end
%% num of each class
numGroupAll={numG1Group,numG2Group,numG3Group,numG4Group,numG5Group,numG6Group,numG7Group,numG8Group,numG9Group};
for curGroup=1:9
    curNumGroup=numGroupAll{curGroup};
    h1=figure; 
    boxPlotCellArray(curNumGroup,nameList,1,false,true);
    ylabel(['Number of adhesions in group ' num2str(curGroup) ' (#)'])
    title(['Number of adhesions in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'numGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-dtiff')
    tableAdhNumGroupEach=table(curNumGroup,'RowNames',nameList);
    writetable(tableAdhNumGroupEach,[dataPath filesep 'numGroup' num2str(curGroup) '.csv'])
end
%% meanRelative population
meanRelPopGroupAll={meanRelativePopG1Group,meanRelativePopG2Group,meanRelativePopG3Group,...
    meanRelativePopG4Group,meanRelativePopG5Group,meanRelativePopG6Group,...
    meanRelativePopG7Group,meanRelativePopG8Group,meanRelativePopG9Group};
for curGroup=1:9
    curmeanPopGroup=meanRelPopGroupAll{curGroup};
    h1=figure; 
    boxPlotCellArray(curmeanPopGroup,nameList,1,false,true);
    ylabel(['Mean relative adhesion population in group ' num2str(curGroup) ' (1)'])
    title(['Mean relative adhesion population in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],'-dtiff')

    tableMeanPopGroupEach=table(curmeanPopGroup,'RowNames',nameList);
    writetable(tableMeanPopGroupEach,[dataPath filesep 'meanRelPopGroup' num2str(curGroup) '.csv'])
end
%% RAtio between g1 and g2
ratioG1G2=cellfun(@(x,y) y./x,meanRelativePopG1Group,meanRelativePopG2Group,'unif',false);
h1=figure; 
boxPlotCellArray(ratioG1G2,nameList,1,false,true);
ylabel(['Mean relative adhesion population of g2 compared to g1'])
title('Mean relative adhesion population of g2 compared to g1')
hgexport(h1,[figPath filesep 'meanRelPopG2overG1'],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'meanRelPopG2overG1'],'-v7.3')
print(h1,[figPath filesep 'meanRelPopG2overG1'],'-dtiff')

tableRelPopG2overG1=table(ratioG1G2,'RowNames',nameList);
writetable(tableRelPopG2overG1,[dataPath filesep 'meanRelPopG2overG1' num2str(curGroup) '.csv'])

%% mean adhesion density - significantly higher NA density in WT compared to R8
meanAdhDensityGroupAll={meanAdhDensityG1Group,meanAdhDensityG2Group,meanAdhDensityG3Group,...
    meanAdhDensityG4Group,meanAdhDensityG5Group,meanAdhDensityG6Group,...
    meanAdhDensityG7Group,meanAdhDensityG8Group,meanAdhDensityG9Group};
for curGroup=1:9
    curAdhDenGroup=meanAdhDensityGroupAll{curGroup};
    h1=figure; 
    boxPlotCellArray(curAdhDenGroup,nameList,1,false,true);
    ylabel(['Mean adhesion density in group ' num2str(curGroup) ' (num/um^2)'])
    title(['Mean adhesion density in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-dtiff')

    tableAdhDenGroupEach=table(curAdhDenGroup,'RowNames',nameList);
    writetable(tableAdhDenGroupEach,[dataPath filesep 'meanAdhDenGroup' num2str(curGroup) '.csv'])
end
%% cell area - 
    h1=figure; 
    boxPlotCellArray(cellAreaGroup,nameList,1,false,true);
    ylabel(['Cell spread area (um^2)'])
    title(['Cell spread area'])
    hgexport(h1,[figPath filesep 'cellArea'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'cellArea'],'-v7.3')
    print(h1,[figPath filesep 'cellArea'],'-dtiff')

    tableCellArea=table(cellAreaGroup,'RowNames',nameList);
    writetable(tableCellArea,[dataPath filesep 'cellArea.csv'])
%% FA area - 
    h1=figure; 
    boxPlotCellArray(FAareaGroup,nameList,1,false,true);
    ylabel(['Focal adhesion area (um^2)'])
    title(['Focal adhesion area'])
    hgexport(h1,[figPath filesep 'faArea'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faArea'],'-v7.3')
    print(h1,[figPath filesep 'faArea'],'-dtiff')

    tableFAArea=table(FAareaGroup,'RowNames',nameList);
    writetable(tableFAArea,[dataPath filesep 'faArea.csv'])
%% NA density
    h1=figure; 
    boxPlotCellArray(NADensityGroup,nameList,1,false,true);
    ylabel(['NA density (#/um^2)'])
    title(['NA density'])
    hgexport(h1,[figPath filesep 'naDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'naDensity'],'-v7.3')
    print(h1,[figPath filesep 'naDensity'],'-dtiff')

    tableNADensity=table(NADensityGroup,'RowNames',nameList);
    writetable(tableNADensity,[dataPath filesep 'naDensity.csv'])
%% FA density
    h1=figure; 
    FADensityGroup=cellfun(@(x,y) x./y,numPureFAsGroup,cellAreaGroup,'unif',false);
    boxPlotCellArray(FADensityGroup,nameList,1,false,true);
    ylabel(['FA density (#/um^2)'])
    title(['FA density'])
    hgexport(h1,[figPath filesep 'faDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faDensity'],'-v7.3')
    print(h1,[figPath filesep 'faDensity'],'-dtiff')

    tableFADensity=table(FADensityGroup,'RowNames',nameList);
    writetable(tableFADensity,[dataPath filesep 'faDensity.csv'])
%% numPureFAsGroup
    h1=figure; 
    boxPlotCellArray(numPureFAsGroup,nameList,1,false,true);
    ylabel(['The number of pure FAs in a cell (#/cell)'])
    title(['The number of pure FAs in a cell'])
    hgexport(h1,[figPath filesep 'numFAs'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numFAs'],'-v7.3')
    print(h1,[figPath filesep 'numFAs'],'-dtiff')

    tableNumPureFAs=table(numPureFAsGroup,'RowNames',nameList);
    writetable(tableNumPureFAs,[dataPath filesep 'numFAs.csv'])
%% maturingRatioGroup
    h1=figure; 
    boxPlotCellArray(maturingRatioGroup,nameList,1,false,true);
    ylabel(['Maturing ratio (1)'])
    title(['Maturing ratio of NAs to FAs'])
    hgexport(h1,[figPath filesep 'maturingRatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatio'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatio'],'-dtiff')

    tableMaturingRatio=table(maturingRatioGroup,'RowNames',nameList);
    writetable(tableMaturingRatio,[dataPath filesep 'maturingRatio.csv'])
%% nucleatingNARatioGroup
    h1=figure; 
    boxPlotCellArray(nucleatingNARatioGroup,nameList,1,false,true);
    ylabel(['nucleating NA Ratio (1)'])
    title(['newly nucleating NA Ratio'])
    hgexport(h1,[figPath filesep 'nucleatingNARatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'nucleatingNARatio'],'-v7.3')
    print(h1,[figPath filesep 'nucleatingNARatio'],'-dtiff')

    tableNucleatingNARatio=table(nucleatingNARatioGroup,'RowNames',nameList);
    writetable(tableNucleatingNARatio,[dataPath filesep 'nucleatingNARatio.csv'])
%% assemRateGroup
    h1=figure; 
    boxPlotCellArray(assemRateGroup,nameList,1,false,true);
    ylabel(['assemRateGroup (1)'])
    title(['assemRateGroup'])
    hgexport(h1,[figPath filesep 'assemRate'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'assemRate'],'-v7.3')
    print(h1,[figPath filesep 'assemRate'],'-dtiff')

    tableAssemRateGroup=table(assemRateGroup,'RowNames',nameList);
    writetable(tableAssemRateGroup,[dataPath filesep 'assemRate.csv'])
%% Plotting each - peakGroup - all classes - usually not interesting: there is no lag.
for curGroup=1:9
    peakLagGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),peakGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(peakLagGroupEach,nameList,1,false,true);
    ylabel(['Peak lag in group ' num2str(curGroup) ' (sec)'])
    title(['Peak lag in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'peakLagG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'peakLagG' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'peakLagG' num2str(curGroup)],'-dtiff')

    tablePeakLagGroupEach=table(peakLagGroupEach,'RowNames',nameList);
    writetable(tablePeakLagGroupEach,[dataPath filesep 'peakLagG' num2str(curGroup) '.csv'])
end
%% Plotting each - endTimeGroup - all classes - usually not interesting: there is no lag.
for curGroup=1:9
    endLagGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),endTimeGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(endLagGroupEach,nameList,1,false,true);
    ylabel(['Ending Time lag in group ' num2str(curGroup) ' (sec)'])
    title(['Ending Time lag in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'endLagG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'endLagG' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'endLagG' num2str(curGroup)],'-dtiff')

    tableEndLagGroupEach=table(endLagGroupEach,'RowNames',nameList);
    writetable(tableEndLagGroupEach,[dataPath filesep 'endLagG' num2str(curGroup) '.csv'])
end
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
