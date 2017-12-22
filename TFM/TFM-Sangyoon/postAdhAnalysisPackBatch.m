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
            [~,finalFolder] = fileparts(curPathProject);
            pathAnalysisAll{ii} = curPathProject;
            groupNames{ii} = finalFolder;
        end
    end
    specificName = strjoin(groupNames);
else
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
    specificName=fileSFolders(16:end);
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);

for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummary' specificName '/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary' specificName '/Data'];
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
        curCellAreas = arrayfun(@(x) x.cellArea,curFocalAdhInfo);
        meanCellArea(k) = mean(curCellAreas);
        
        nafaStruct=load(curAnalProc.outFilePaths_{3,iAdhChan});
        curNADensity = nafaStruct.NADensity;
        meanNADensity(k) = mean(curNADensity);
        curFADensity = nafaStruct.FADensity;
        meanFADensity(k) = mean(curFADensity);

        maturingStruct=load(curAnalProc.outFilePaths_{4,iAdhChan});
        meanLifeTimeAll(k)=mean(maturingStruct.lifeTimeAll);
        meanLifeTimeFailingNAs(k)=mean(maturingStruct.lifeTimeNAfailing);
        meanLifeTimeMaturingNAs(k)=mean(maturingStruct.lifeTimeNAmaturing);
        meanMaturingRatio(k)=mean(maturingStruct.maturingRatio);
        
        assemRateStruct=load(curAnalProc.outFilePaths_{5,iAdhChan});
        meanAssemRate(k)=mean(assemRateStruct.assemRateCell);
        meanDisassemRate(k)=mean(assemRateStruct.disassemRateCell);
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
            
            nameTitle=['initialLag_Class' num2str(pp)];
            initRiseStruct = load([initDataPath filesep nameTitle],'initialLagTogetherAdjusted','nameList2');   
            curInitTimeLag=initRiseStruct.initialLagTogetherAdjusted;
            initRiseAgainstForceEachClass{pp}{k} = curInitTimeLag;
            
            nameTitle=['peakLag_Class' num2str(pp)];
            peakStruct = load([initDataPath filesep nameTitle],'peakLagTogetherAdjusted','nameList2');   
            curPeakTimeLag=peakStruct.peakLagTogetherAdjusted;
            peakTimeAgainstForceEachClass{pp}{k} = curPeakTimeLag;
            
            nameTitle=['endingLag_Class' num2str(pp)];
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

    meanAdhDensityG1Group{ii,1}=meanRelativePopG1;
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
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
%% Plotting each - SE
nameList=groupNames'; %{'pLVX' 'P29S'};
SEGroupCell = cellfun(@(x) cell2mat(x),SEGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEGroupCell,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy')
hgexport(h1,strcat(figPath,'/strainEnergyTotal'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyTotal'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyTotal.tif'),'-dtiff')
%% Strain energy density
SEDenGroupCell = cellfun(@(x) cell2mat(x),SEDenGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCell,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density')
hgexport(h1,strcat(figPath,'/SEDensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensity'),'-v7.3')
print(h1,strcat(figPath,'/SEDensity.tif'),'-dtiff')
%% Total force
totForceCell = cellfun(@(x) cell2mat(x),totalForceGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(totForceCell,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in force blobs')
hgexport(h1,strcat(figPath,'/totForce'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForce'),'-v7.3')
print(h1,strcat(figPath,'/totForce.tif'),'-dtiff')
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
