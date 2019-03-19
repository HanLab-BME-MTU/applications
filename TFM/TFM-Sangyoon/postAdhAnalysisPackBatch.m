%% open necessary MLs
MLdirect=false;
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
groupNames=[];
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
            [curPathProject2,finalFolder] = fileparts(curPathProject);
            pathAnalysisAll{ii} = curPathProject;
            if isempty(finalFolder)
                [~,finalFolder] = fileparts(curPathProject2);
            end
            groupNames{ii} = finalFolder;
            MLFileNamesAll{ii} = curMLFile;
            MLNames{ii} = curMLFile;
            MLdirect=true;
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
    MLFileNamesAll = selectedFolders.MLNames;%'movieList.mat';
%     for k=1:numel(pathAnalysisAll)
%         MLFileNamesAll{k} = selectedFolders.MLNames{k};%'movieList.mat';
%     end
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);

for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLFileNamesAll{k}]);
end
%% Output
% rootAnalysis = fileparts(pathAnalysisAll{1});
rootAnalysis = pathAnalysisAll{1};
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
halfBccGroup = cell(numConditions,1);

earlyAmpSlopeGroup = cell(numConditions,1);
earlyForceSlopeGroup = cell(numConditions,1);

numPureFAsGroup = cell(numConditions,1);
cellAreaGroup = cell(numConditions,1);
NADensityGroup = cell(numConditions,1);
FADensityGroup = cell(numConditions,1);
lifeTimeAllGroup = cell(numConditions,1);
lifeTimeFailingNAsGroup = cell(numConditions,1);
lifeTimeMaturingNAsGroup = cell(numConditions,1);
maturingRatioGroup = cell(numConditions,1);
maturingRatioNAtoFAGroup = cell(numConditions,1);
maturingRatioFCtoFAGroup = cell(numConditions,1);
lifeTimeFAsGroup = cell(numConditions,1);
stableNAFCratioGroup = cell(numConditions,1);
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

intensitiesInNAsGroup = cell(numConditions,1);
intensitiesInFCsGroup = cell(numConditions,1);
intensitiesInFAsGroup = cell(numConditions,1);

mainBccPeakValuesGroupGroup = cell(numConditions,1);
mainTimeToPeakGroupGroup = cell(numConditions,1);
sideBccPeakValuesGrouppGroup = cell(numConditions,1);
sideTimeToPeakGroupGroup = cell(numConditions,1);

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
    
    allFAareaGroup=cell(N(ii),1);
    allFAlengthGroup=cell(N(ii),1);
    meanNumPureFAs=zeros(N(ii),1);
    meanCellArea=zeros(N(ii),1);
    naDensityAllFrames=cell(N(ii),1);
    faDensityAllFrames=cell(N(ii),1);
    meanLifeTimeAll=zeros(N(ii),1);
    meanLifeTimeFailingNAs=zeros(N(ii),1);
    meanLifeTimeMaturingNAs=zeros(N(ii),1);
    meanMaturingRatio=zeros(N(ii),1);
    meanMaturingRatioNAtoFA=NaN(N(ii),1);
    meanMaturingRatioFCtoFA=NaN(N(ii),1);
    lifeTimeFAsAll=cell(N(ii),1);
    meanStableNAFCratio=NaN(N(ii),1);
    meanAssemRate=cell(N(ii),1);
    meanDisassemRate=cell(N(ii),1);
    meanNucleatingNARatio=cell(N(ii),1);
    meanDisassemNARatio=cell(N(ii),1);
    
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
    
    meanAdhDensityG1=cell(N(ii),1);
    meanAdhDensityG2=cell(N(ii),1);
    meanAdhDensityG3=cell(N(ii),1);
    meanAdhDensityG4=cell(N(ii),1);
    meanAdhDensityG5=cell(N(ii),1);
    meanAdhDensityG6=cell(N(ii),1);
    meanAdhDensityG7=cell(N(ii),1);
    meanAdhDensityG8=cell(N(ii),1);
    meanAdhDensityG9=cell(N(ii),1);

    intensitiesInNAs=cell(N(ii),1);   
    intensitiesInFCs=cell(N(ii),1);   
    intensitiesInFAs=cell(N(ii),1);   
    
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
        allFAareaGroup{k} = cell2mat(arrayfun(@(x) x.area',curFocalAdhInfo,'unif',false));
%         allFAareaGroup{k} = mean(curAllFAareaAllFrames);
        allFAlengthGroup{k} = cell2mat(arrayfun(@(x) x.length',curFocalAdhInfo,'unif',false));
%         curMedianFAlengthAllFrames = arrayfun(@(x) x.medianLength,curFocalAdhInfo);
%         allFAlengthGroup(k) = mean(curMedianFAlengthAllFrames);
        curNumPureFAs = arrayfun(@(x) x.numberPureFA,curFocalAdhInfo);
        meanNumPureFAs(k) = mean(curNumPureFAs);
        curCellAreas = arrayfun(@(x) x.cellArea,curFocalAdhInfo); % in um2
        meanCellArea(k) = mean(curCellAreas);
        
        nafaStruct=load(curAnalProc.outFilePaths_{3,iAdhChan});
        naDensityAllFrames{k} = nafaStruct.NADensity;
%         meanNADensity(k) = mean(curNADensity);
        faDensityAllFrames{k} = nafaStruct.FADensity;
%         meanFADensity(k) = meanNumPureFAs(k)/meanCellArea(k); % in num/um2

        maturingStruct=load(curAnalProc.outFilePaths_{4,iAdhChan});
        meanLifeTimeAll(k)=mean(maturingStruct.lifeTimeAll);
        meanLifeTimeFailingNAs(k)=mean(maturingStruct.lifeTimeNAfailing);
        meanLifeTimeMaturingNAs(k)=mean(maturingStruct.lifeTimeNAmaturing);
        meanMaturingRatio(k)=mean(maturingStruct.maturingRatio);
        if isfield(maturingStruct,'maturingRatioNAtoFA') % This is the newer fields
            meanMaturingRatioNAtoFA(k)=mean(maturingStruct.maturingRatioNAtoFA);
            meanMaturingRatioFCtoFA(k)=mean(maturingStruct.maturingRatioFCtoFA);
            meanStableNAFCratio(k)=mean(maturingStruct.stableNAFCratio);
            lifeTimeFAsAll{k}=maturingStruct.lifeTimeFAs;
        else
            meanMaturingRatioNAtoFA(k)=NaN;
            meanMaturingRatioFCtoFA(k)=NaN;
            meanStableNAFCratio(k)=NaN;
            lifeTimeFAsAll{k}=NaN;
        end
        
        assemRateStruct=load(curAnalProc.outFilePaths_{5,iAdhChan});
        meanAssemRate{k}=(assemRateStruct.assemRateCell);
        meanDisassemRate{k}=(assemRateStruct.disassemRateCell);
        meanNucleatingNARatio{k} = (assemRateStruct.nucleationRatio);
        meanDisassemNARatio{k} = (assemRateStruct.disassemblingNARatio);
        
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
        
        curPresenceStruct = load(classProc.outFilePaths_{5,iAdhChan},'presenceAll');
        curPresence = curPresenceStruct.presenceAll;
        
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
        meanAdhDensityG1{k} = sum(idsStruct.idGroup1 & curPresence)./curCellAreas'; meanAdhDensityG1{k}(1)=[];
        meanAdhDensityG2{k} = sum(idsStruct.idGroup2 & curPresence)./curCellAreas'; meanAdhDensityG2{k}(1)=[];
        meanAdhDensityG3{k} = sum(idsStruct.idGroup3 & curPresence)./curCellAreas'; meanAdhDensityG3{k}(1)=[];
        meanAdhDensityG4{k} = sum(idsStruct.idGroup4 & curPresence)./curCellAreas'; meanAdhDensityG4{k}(1)=[];
        meanAdhDensityG5{k} = sum(idsStruct.idGroup5 & curPresence)./curCellAreas'; meanAdhDensityG5{k}(1)=[];
        meanAdhDensityG6{k} = sum(idsStruct.idGroup6 & curPresence)./curCellAreas'; meanAdhDensityG6{k}(1)=[];
        meanAdhDensityG7{k} = sum(idsStruct.idGroup7 & curPresence)./curCellAreas'; meanAdhDensityG7{k}(1)=[];
        meanAdhDensityG8{k} = sum(idsStruct.idGroup8 & curPresence)./curCellAreas'; meanAdhDensityG8{k}(1)=[];
        meanAdhDensityG9{k} = sum(idsStruct.idGroup9 & curPresence)./curCellAreas'; meanAdhDensityG9{k}(1)=[];
%         % Individual adh density per each group per frame -- I will do
%         this later SH 11/14/2018
%         s = struct2table(tracksNA);
%         startingFrameExtra = s.startingFrameExtra;
%         endingFrameExtra = s.endingFrameExtra;
%         validTracks = validState & startingFrameExtra <= iFrame & endingFrameExtra >= iFrame;                    
%         numG1PerFrame = 1
%         adhDensityG1{k}=1
        
        % About process 9
        iTheOtherProc = 9;
        tOtherProc = curFAPackage.getProcess(iTheOtherProc);
        if ~isempty(tOtherProc)
            iOther = input('Another channel other than the main adhesion channel? (default: 3): ');
            if isempty(iOther); iOther=3; end

            intenGroupStruct=load(tOtherProc.outFilePaths_{2,iOther});
            intenGroup=intenGroupStruct.intensityGroup; 
            intensitiesInNAs{k} = intenGroup{1};
            intensitiesInFCs{k} = intenGroup{2};
            intensitiesInFAs{k} = intenGroup{3};
        else
            intensitiesInNAs{k} = [];
            intensitiesInFCs{k} = [];
            intensitiesInFAs{k} = [];
        end        
        % Other feature related properties are calculated in the step 11
        
    end
    
    
    initRiseAgainstForceEachClass=cell(numClasses,1);
    peakTimeAgainstForceEachClass=cell(numClasses,1);
    endTimeAgainstForceEachClass=cell(numClasses,1);
    halfBccTogetherEachClass=cell(numClasses,1);
    
    earlyAmpSlopeEachClass=cell(numClasses,1);
    earlyForceSlopeEachClass=cell(numClasses,1);
    
    mainBccPeakValuesGroup=cell(2,1);
    mainTimeToPeakGroup=cell(2,1);
    sideBccPeakValuesGroup=cell(2,1);
    sideTimeToPeakGroup=cell(2,1);
    
    iForceSlave = 1;
    iInitRiseProc = 11;
    for pp=1:numClasses
        for k=1:N(ii)
            curMovie=curMovies{k};
            iCurPack = curMovie.getPackageIndex('FocalAdhesionPackage');
            curFAPackage = curMovie.getPackage(iCurPack);
            % Initial rise
            initRiseProc = curFAPackage.getProcess(iInitRiseProc);
            % I decided to get the adjusted time lags
            if ~isempty(initRiseProc)
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

                nameTitle=['halfBccTogetherAdjusted' num2str(pp)];
                halfBccTogetherStruct = load([initDataPath filesep nameTitle],'halfBccTogetherAdjusted','nameList2');    
                halfBccTogetherEachClass{pp}{k} = halfBccTogetherStruct.halfBccTogetherAdjusted;

                earlyAmpSlopeStr=load([initDataPath filesep 'earlyAmpSlopeAllGroups.mat'],'earlyAmpSlope');
                earlyAmpSlopeEachClass{pp}{k} = earlyAmpSlopeStr.earlyAmpSlope{pp};

                earlyForceSlopeStr=load([initDataPath filesep 'earlyForceSlopeAllGroups.mat'],'earlyForceSlope');
                earlyForceSlopeEachClass{pp}{k} = earlyForceSlopeStr.earlyForceSlope{pp};


                if ismember(pp,[1 2])
                    try
                        mainBccPeakValuesGroupStruct = load([initDataPath filesep 'mainBccPeakValuesGroup.mat'],'mainBccPeakValuesGroup');
                        mainBccPeakValuesGroup{pp}{k} = mainBccPeakValuesGroupStruct.mainBccPeakValuesGroup{pp};
                        mainTimeToPeakGroupStruct = load([initDataPath filesep 'mainTimeToPeakGroup.mat'],'mainTimeToPeakGroup');
                        mainTimeToPeakGroup{pp}{k} = mainTimeToPeakGroupStruct.mainTimeToPeakGroup{pp};
                        sideBccPeakValuesGroupStruct = load([initDataPath filesep 'sideBccPeakValuesGroup.mat'],'sideBccPeakValuesGroup');
                        sideBccPeakValuesGroup{pp}{k} = sideBccPeakValuesGroupStruct.sideBccPeakValuesGroup{pp};
                        sideTimeToPeakGroupStruct = load([initDataPath filesep 'sideTimeToPeakGroup.mat'],'sideTimeToPeakGroup');
                        sideTimeToPeakGroup{pp}{k} = sideTimeToPeakGroupStruct.sideTimeToPeakGroup{pp};
                    catch
                        try
                            mainBccPeakValuesGroupStruct = load([initDataPath filesep 'mainBccPeakValues-G' num2str(pp) '.mat'],'mainBccPeakValues');
                            mainBccPeakValuesGroup{pp}{k} = mainBccPeakValuesGroupStruct.mainBccPeakValues;
                            mainTimeToPeakGroupStruct = load([initDataPath filesep 'mainTimeToPeak-G' num2str(pp) '.mat']);%,'mainTimeToPeakGroup');
                            mainTimeToPeakGroup{pp}{k} = mainTimeToPeakGroupStruct.mainTimeToPeak;
                            sideBccPeakValuesGroupStruct = load([initDataPath filesep 'sideBccPeakValues-G' num2str(pp) '.mat'],'sideBccPeakValues');
                            sideBccPeakValuesGroup{pp}{k} = sideBccPeakValuesGroupStruct.sideBccPeakValues;
                            sideTimeToPeakGroupStruct = load([initDataPath filesep 'sideTimeToPeak-G' num2str(pp) '.mat'],'sideTimeToPeak');
                            sideTimeToPeakGroup{pp}{k} = sideTimeToPeakGroupStruct.sideTimeToPeak;
                        catch
                            mainBccPeakValuesGroup{pp}{k} = [];
                            mainTimeToPeakGroup{pp}{k} = [];
                            sideBccPeakValuesGroup{pp}{k} = [];
                            sideTimeToPeakGroup{pp}{k} = [];
                        end
                    end
                end
            end
        end
    end
    FAareaGroup{ii,1}=allFAareaGroup;
    FAlengthGroup{ii,1}=allFAlengthGroup;
    numPureFAsGroup{ii,1}=meanNumPureFAs;
    cellAreaGroup{ii,1}=meanCellArea;
    NADensityGroup{ii,1}=naDensityAllFrames;
    FADensityGroup{ii,1}=faDensityAllFrames;
    lifeTimeAllGroup{ii,1}=meanLifeTimeAll;
    lifeTimeFailingNAsGroup{ii,1}=meanLifeTimeFailingNAs;
    lifeTimeMaturingNAsGroup{ii,1}=meanLifeTimeMaturingNAs;
    maturingRatioGroup{ii,1}=meanMaturingRatio;
    maturingRatioNAtoFAGroup{ii,1}=meanMaturingRatioNAtoFA;
    maturingRatioFCtoFAGroup{ii,1}=meanMaturingRatioFCtoFA;
    stableNAFCratioGroup{ii,1}=meanStableNAFCratio;
    lifeTimeFAsGroup{ii,1} = lifeTimeFAsAll;
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

    intensitiesInNAsGroup{ii,1}=intensitiesInNAs;
    intensitiesInFCsGroup{ii,1}=intensitiesInFCs;
    intensitiesInFAsGroup{ii,1}=intensitiesInFAs;

    if ~isempty(initRiseProc)
        initRiseGroup{ii,1}=initRiseAgainstForceEachClass;
        peakGroup{ii,1}=peakTimeAgainstForceEachClass;
        endTimeGroup{ii,1}=endTimeAgainstForceEachClass;

        halfBccGroup{ii,1}=halfBccTogetherEachClass;

        earlyAmpSlopeGroup{ii,1}=earlyAmpSlopeEachClass;
        earlyForceSlopeGroup{ii,1}=earlyForceSlopeEachClass;

        mainBccPeakValuesGroupGroup{ii,1}=mainBccPeakValuesGroup;
        mainTimeToPeakGroupGroup{ii,1}=mainTimeToPeakGroup;
        sideBccPeakValuesGrouppGroup{ii,1}=sideBccPeakValuesGroup;
        sideTimeToPeakGroupGroup{ii,1}=sideTimeToPeakGroup;
    end
end
disp('Done')
%% setting up group name
if MLdirect && isempty(groupNames)
%     if strcmp(MLFileNamesAll{1}(end-7:end-4),'List')
    groupNames=MLFileNamesAll;
elseif isempty(groupNames)
    for ii=1:numConditions
        [pathFolder, finalFolder]=fileparts(pathAnalysisAll{ii});
        if isempty(finalFolder)
            [~, finalFolder]=fileparts(pathFolder);
        end
        groupNames{ii} = finalFolder;
    end
end
% for ii=1:numConditions
%     [~, finalFolder]=fileparts(pathAnalysisAll{ii});
%     groupNames{ii} = finalFolder;
% end
%% Plotting each - initRiseGroup - all classes
nameList=groupNames'; %{'pLVX' 'P29S'};
if ~isempty(initRiseProc)
    for curGroup=1:9
        initRiseGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),initRiseGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(initRiseGroupEach,nameList,1,false,true);
        ylabel(['Initial rise in group ' num2str(curGroup) ' (sec)'])
        title(['Initial rise in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'initialRizeG' num2str(curGroup)],'-dtiff')

        tableHalfBccGroup=table(initRiseGroupEach,'RowNames',nameList);
        writetable(tableHalfBccGroup,[dataPath filesep 'initialRizeG' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% halfBccTime
if ~isempty(initRiseProc)
    for curGroup=1:9
        halfBccGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),halfBccGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(halfBccGroupEach,nameList,1,false,true);
        ylabel(['halfBccGroup ' num2str(curGroup) ' (sec)'])
        title(['halfBccGroup in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'halfBccGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'halfBccGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'halfBccGroup' num2str(curGroup)],'-dtiff')

        tableHalfBccGroup=table(halfBccGroupEach,'RowNames',nameList);
        writetable(tableHalfBccGroup,[dataPath filesep 'halfBccGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% earlyAmpSlope
if ~isempty(initRiseProc)
    for curGroup=1:9
        earlyAmpSlopeGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),earlyAmpSlopeGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(earlyAmpSlopeGroupEach,nameList,1,false,true);
        ylabel(['earlyAmpSlope ' num2str(curGroup) ' (sec)'])
        title(['earlyAmpSlope in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],'-dtiff')

        tableEarlyAmpSlope=table(earlyAmpSlopeGroupEach,'RowNames',nameList);
        writetable(tableEarlyAmpSlope,[dataPath filesep 'earlyAmpSlope' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% earlyAmpSlope - G1G2 only
if ~isempty(initRiseProc)
    curGroup=1;
    earlyAmpSlopeG1 = cellfun(@(x) cell2mat(x{curGroup}'),earlyAmpSlopeGroup,'unif',false);
    curGroup=2;
    earlyAmpSlopeG2 = cellfun(@(x) cell2mat(x{curGroup}'),earlyAmpSlopeGroup,'unif',false);
    earlyAmpSlopeG1G2=cell(numel(earlyAmpSlopeG1)*2,1);
    nameListG1G2=cell(numel(earlyAmpSlopeG1)*2,1);
    for ii=1:numel(earlyAmpSlopeG1)
        earlyAmpSlopeG1G2(2*(ii-1)+1)=earlyAmpSlopeG1(ii);
        earlyAmpSlopeG1G2(2*(ii))=earlyAmpSlopeG2(ii);
    end
    if MLdirect
        for ii=1:numel(earlyAmpSlopeG1)
            nameListG1G2{2*(ii-1)+1}=[nameList{ii}(10:end-4) '-G1'];
            nameListG1G2{2*(ii)}=[nameList{ii}(10:end-4) '-G2'];
        end
    else
        for ii=1:numel(earlyAmpSlopeG1)
            nameListG1G2{2*(ii-1)+1}=[nameList{ii} '-G1'];
            nameListG1G2{2*(ii)}=[nameList{ii} '-G2'];
        end
    end
    h1=figure; 
    boxPlotCellArray(earlyAmpSlopeG1G2,nameListG1G2,1,false,true);
    ylabel('earlyAmpSlope (a.u./sec)')
    title('earlyAmpSlope in G1 vs. G2')
    hgexport(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],'-dtiff')
end            
%% earlyForceSlope
if ~isempty(initRiseProc)
    for curGroup=1:9
        earlyForceSlopeGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),earlyForceSlopeGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(earlyForceSlopeGroupEach,nameList,1,false,true);
        ylabel(['earlyForceSlope ' num2str(curGroup) ' (sec)'])
        title(['earlyForceSlope in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],'-dtiff')

        tableEarlyForceSlope=table(earlyForceSlopeGroupEach,'RowNames',nameList);
        writetable(tableEarlyForceSlope,[dataPath filesep 'earlyForceSlope' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% earlyForceSlope - G1G2 only
if ~isempty(initRiseProc)
    curGroup=1;
    earlyForceSlopeG1 = cellfun(@(x) cell2mat(x{curGroup}'),earlyForceSlopeGroup,'unif',false);
    curGroup=2;
    earlyForceSlopeG2 = cellfun(@(x) cell2mat(x{curGroup}'),earlyForceSlopeGroup,'unif',false);
    earlyForceSlopeG1G2=cell(numel(earlyAmpSlopeG1)*2,1);
    for ii=1:numel(earlyForceSlopeG1)
        earlyForceSlopeG1G2(2*(ii-1)+1)=earlyForceSlopeG1(ii);
        earlyForceSlopeG1G2(2*(ii))=earlyForceSlopeG2(ii);
    end

    h1=figure; 
    % barPlotCellArray(earlyForceSlopeG1G2,nameListG1G2);
    boxPlotCellArray(earlyForceSlopeG1G2,nameListG1G2,1,false,true);
    ylabel('earlyForceSlope (a.u./sec)')
    title('earlyForceSlope in G1 vs. G2')
    hgexport(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],'-dtiff')
end            
%% Plotting each - mainTimeToPeakGroupGroup - G1 and G2
if ~isempty(initRiseProc)
    for curGroup=1:2
        mainTimeToPeakGroupGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),mainTimeToPeakGroupGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(mainTimeToPeakGroupGroupEach,nameList,1,false,true);
        ylabel(['mainTimeToPeakGroupGroupEach ' num2str(curGroup) ' (sec)'])
        title(['mainTimeToPeakGroupGroupEach ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],'-dtiff')
    end
    %% Plotting each - mainBccPeakValuesGroupGroup - G1 and G2
    for curGroup=1:2
        mainBccPeakValuesGroupGrouppEach = cellfun(@(x) cell2mat(x{curGroup}'),mainBccPeakValuesGroupGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(mainBccPeakValuesGroupGrouppEach,nameList,1,false,true);
        ylabel(['mainBccPeakValuesGroupGrouppEach ' num2str(curGroup) ' (sec)'])
        title(['mainBccPeakValuesGroupGrouppEach ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],'-dtiff')
    end
    %% Plotting each - sideTimeToPeakGroupGroup - G1 and G2
    for curGroup=1:2
        sideTimeToPeakGroupGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),sideTimeToPeakGroupGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(sideTimeToPeakGroupGroupEach,nameList,1,false,true);
        ylabel(['sideTimeToPeakGroupGroup ' num2str(curGroup) ' (sec)'])
        title(['sideTimeToPeakGroupGroup ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],'-dtiff')
    end
    %% Plotting each - sideBccPeakValuesGrouppGroup - G1 and G2
    for curGroup=1:2
        sideBccPeakValuesGrouppGrouppEach = cellfun(@(x) cell2mat(x{curGroup}'),sideBccPeakValuesGrouppGroup,'unif',false);
        h1=figure; 
        boxPlotCellArray(sideBccPeakValuesGrouppGrouppEach,nameList,1,false,true);
        ylabel(['sideBccPeakValuesGrouppGroup ' num2str(curGroup) ' (sec)'])
        title(['sideBccPeakValuesGrouppGroup ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],'-dtiff')
    end
end
%% num of each class
numGroupOnlyTwo={numG1Group,numG2Group}; %,numG3Group,numG4Group,numG5Group,numG6Group,numG7Group,numG8Group,numG9Group};
for curGroup=1:2
    curNumGroup=numGroupOnlyTwo{curGroup};
    h1=figure; 
    boxPlotCellArray(curNumGroup,nameList,1,false,true);
    ylabel(['Number of adhesions in group ' num2str(curGroup) ' (#)'])
    title(['Number of adhesions in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'numGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-dtiff')
    tableAdhNumGroupEach=table(curNumGroup,'RowNames',nameList);
    writetable(tableAdhNumGroupEach,[dataPath filesep 'numGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
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
    writetable(tableMeanPopGroupEach,[dataPath filesep 'meanRelPopGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
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
writetable(tableRelPopG2overG1,[dataPath filesep 'meanRelPopG2overG1' num2str(curGroup) '.csv'],'WriteRowNames',true)

%% mean adhesion density - significantly higher NA density in WT compared to R8
meanAdhDensityGroupAll={meanAdhDensityG1Group,meanAdhDensityG2Group,meanAdhDensityG3Group,...
    meanAdhDensityG4Group,meanAdhDensityG5Group,meanAdhDensityG6Group,...
    meanAdhDensityG7Group,meanAdhDensityG8Group,meanAdhDensityG9Group};
for curGroup=1:9
    curAdhDenGroup=meanAdhDensityGroupAll{curGroup};
    h1=figure; 
    curAdhDenGroupCell = cellfun(@(x) cell2mat(x'),curAdhDenGroup,'unif',false);

    boxPlotCellArray(curAdhDenGroupCell,nameList,1,false,true);
    ylabel(['Mean adhesion density in group ' num2str(curGroup) ' (num/um^2)'])
    title(['Mean adhesion density in group ' num2str(curGroup)])
    hgexport(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-dtiff')

    tableAdhDenGroupEach=table(curAdhDenGroup,'RowNames',nameList);
    writetable(tableAdhDenGroupEach,[dataPath filesep 'meanAdhDenGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
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
    writetable(tableCellArea,[dataPath filesep 'cellArea.csv'],'WriteRowNames',true)
%% FA area - 
    FAareaGroupCell = cellfun(@(x) cell2mat(x),FAareaGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(FAareaGroupCell,nameList,1,false,true);
    ylabel(['Focal adhesion area (um^2)'])
    title(['Focal adhesion area'])
    hgexport(h1,[figPath filesep 'faArea'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faArea'],'-v7.3')
    print(h1,[figPath filesep 'faArea'],'-dtiff')

    tableFAArea=table(FAareaGroupCell,'RowNames',nameList);
    writetable(tableFAArea,[dataPath filesep 'faArea.csv'],'WriteRowNames',true)
%% FA length - 
    FAlengthGroupCell = cellfun(@(x) cell2mat(x),FAlengthGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(FAlengthGroupCell,nameList,1,false,true);
    ylabel(['Focal adhesion length (um)'])
    title(['Focal adhesion length'])
    hgexport(h1,[figPath filesep 'faLength'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faLength'],'-v7.3')
    print(h1,[figPath filesep 'faLength'],'-dtiff')

    tableFALength=table(FAlengthGroupCell,'RowNames',nameList);
    writetable(tableFALength,[dataPath filesep 'faLength.csv'],'WriteRowNames',true)
%% NA density
    NADensityGroupCell = cellfun(@(x) cell2mat(x),NADensityGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(NADensityGroupCell,nameList,1,false,true);
    ylabel(['NA density (#/um^2)'])
    title(['NA density'])
    hgexport(h1,[figPath filesep 'naDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'naDensity'],'-v7.3')
    print(h1,[figPath filesep 'naDensity'],'-dtiff')

    tableNADensity=table(NADensityGroup,'RowNames',nameList);
    writetable(tableNADensity,[dataPath filesep 'naDensity.csv'],'WriteRowNames',true)
%% FA density
    FADensityGroupCell = cellfun(@(x) cell2mat(x),FADensityGroup,'unif',false);
    h1=figure; 
%     FADensityGroup=cellfun(@(x,y) x./y,numPureFAsGroup,cellAreaGroup,'unif',false);
    boxPlotCellArray(FADensityGroupCell,nameList,1,false,true);
    ylabel(['FA density (#/um^2)'])
    title(['FA density'])
    hgexport(h1,[figPath filesep 'faDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faDensity'],'-v7.3')
    print(h1,[figPath filesep 'faDensity'],'-dtiff')

    tableFADensity=table(FADensityGroup,'RowNames',nameList);
    writetable(tableFADensity,[dataPath filesep 'faDensity.csv'],'WriteRowNames',true)
%% numPureFAsGroup
    h1=figure; 
    boxPlotCellArray(numPureFAsGroup,nameList,1,false,true);
    ylabel(['The number of pure FAs in a cell (#/cell)'])
    title(['The number of pure FAs in a cell'])
    hgexport(h1,[figPath filesep 'numFAs'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numFAs'],'-v7.3')
    print(h1,[figPath filesep 'numFAs'],'-dtiff')

    tableNumPureFAs=table(numPureFAsGroup,'RowNames',nameList);
    writetable(tableNumPureFAs,[dataPath filesep 'numFAs.csv'],'WriteRowNames',true)
%% maturingRatioGroup
    h1=figure; 
    boxPlotCellArray(maturingRatioGroup,nameList,1,false,true);
    ylabel(['Maturing ratio (1)'])
    title(['Maturing ratio of NAs to FCs'])
    hgexport(h1,[figPath filesep 'maturingRatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatio'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatio'],'-dtiff')

    tableMaturingRatio=table(maturingRatioGroup,'RowNames',nameList);
    writetable(tableMaturingRatio,[dataPath filesep 'maturingRatio.csv'],'WriteRowNames',true)
%% maturingRatioNAtoFAGroup
    h1=figure; 
    boxPlotCellArray(maturingRatioNAtoFAGroup,nameList,1,false,true);
    ylabel(['Maturing ratio from NA to FA (1)'])
    title(['Maturing ratio of NAs to FAs'])
    hgexport(h1,[figPath filesep 'maturingRatioNAtoFA'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatioNAtoFA'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatioNAtoFA'],'-dtiff')

    tableMaturingRatioNAtoFAGroup=table(maturingRatioNAtoFAGroup,'RowNames',nameList);
    writetable(tableMaturingRatioNAtoFAGroup,[dataPath filesep 'maturingRatioNAtoFA.csv'],'WriteRowNames',true)
%% maturingRatioFCtoFAGroup
    h1=figure; 
    boxPlotCellArray(maturingRatioFCtoFAGroup,nameList,1,false,true);
    ylabel('Maturing ratio from FCs to FAs (1)')
    title(['Maturing ratio of FCs to FAs'])
    hgexport(h1,[figPath filesep 'maturingRatioFCtoFA'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatioFCtoFA'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatioFCtoFA'],'-dtiff')

    tableMaturingRatioFCtoFA=table(maturingRatioFCtoFAGroup,'RowNames',nameList);
    writetable(tableMaturingRatioFCtoFA,[dataPath filesep 'maturingRatioFCtoFA.csv'],'WriteRowNames',true)
%% stableNAFCratioGroup
    h1=figure; 
    boxPlotCellArray(stableNAFCratioGroup,nameList,1,false,true);
    ylabel('The ratio of stable NAs and FCs (1)')
    title(['The ratio of stable NAs and FCs'])
    hgexport(h1,[figPath filesep 'stableNAFCratio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'stableNAFCratio'],'-v7.3')
    print(h1,[figPath filesep 'stableNAFCratio'],'-dtiff')

    tableStableNAFCratio=table(stableNAFCratioGroup,'RowNames',nameList);
    writetable(tableStableNAFCratio,[dataPath filesep 'stableNAFCratio.csv'],'WriteRowNames',true)
%% lifeTimeFAsGroup
    lifeTimeFAsGroupCell = cellfun(@(x) cell2mat(x),lifeTimeFAsGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(lifeTimeFAsGroupCell,nameList,curMovie.timeInterval_/60,false,true);
%     barPlotCellArray(lifeTimeFAsGroupCell,nameList,curMovie.timeInterval_/60);
    ylabel('Life time of FAs (min)')
    title(['Life time of FAs'])
    hgexport(h1,[figPath filesep 'lifeTimeFAs'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeFAs'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeFAs'],'-dtiff')

    tableStableNAFCratio=table(stableNAFCratioGroup,'RowNames',nameList);
    writetable(tableStableNAFCratio,[dataPath filesep 'lifeTimeFAs.csv'],'WriteRowNames',true)

%% intensity of the other channel at each adhesion type
    h1=figure; 
    intensitiesInNAsGroupCell = cellfun(@(x) cell2mat(x),intensitiesInNAsGroup,'unif',false);
    intensitiesInFCsGroupCell = cellfun(@(x) cell2mat(x),intensitiesInFCsGroup,'unif',false);
    intensitiesInFAsGroupCell = cellfun(@(x) cell2mat(x),intensitiesInFAsGroup,'unif',false);
    intensityGroupAll=cell(numel(intensitiesInNAsGroup)*3,1);
    intensityGroup={intensitiesInNAsGroupCell, intensitiesInFCsGroupCell, intensitiesInFAsGroupCell};

    nameListAdh={'NA','FC','FA'};
    nameListAdhComb=cell(numel(intensitiesInNAsGroup)*3,1);
    
    for ii=1:numel(intensitiesInNAsGroup)
        p=ii-1;
        for jj=1:3
            nameListAdhComb{3*p+jj,1} = [nameList{ii} '-' nameListAdh{jj}];
            intensityGroupAll{3*p+jj,1} = intensityGroup{jj}{ii};
        end
        
    end
    boxPlotCellArray(intensityGroupAll,nameListAdhComb,1,false,true);
    ylabel(['Fluorescence Intensity (a.u.)'])
    title(['F.I. of the other channel'])
    hgexport(h1,[figPath filesep 'intenTheOtherChannel'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'intenTheOtherChannel'],'-v7.3')
    print(h1,[figPath filesep 'intenTheOtherChannel'],'-dtiff')

    tableIntensityTheOther=table(intensityGroupAll,'RowNames',nameListAdhComb);
    writetable(tableIntensityTheOther,[dataPath filesep 'intenTheOtherChannel.csv'],'WriteRowNames',true)    
%% nucleatingNARatioGroup
    nucleatingNARatioGroupCell = cellfun(@(x) cell2mat(x),nucleatingNARatioGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(nucleatingNARatioGroupCell,nameList,1,false,true);
    ylabel(['nucleating NA Ratio (1)'])
    title(['newly nucleating NA Ratio'])
    hgexport(h1,[figPath filesep 'nucleatingNARatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'nucleatingNARatio'],'-v7.3')
    print(h1,[figPath filesep 'nucleatingNARatio'],'-dtiff')

    tableNucleatingNARatio=table(nucleatingNARatioGroup,'RowNames',nameList);
    writetable(tableNucleatingNARatio,[dataPath filesep 'nucleatingNARatio.csv'],'WriteRowNames',true)
%% assemRateGroup
    assemRateGroupCell = cellfun(@(x) cell2mat(x),assemRateGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(assemRateGroupCell,nameList,1,false,true);
    ylabel(['assemRateGroup (1)'])
    title(['assemRateGroup'])
    hgexport(h1,[figPath filesep 'assemRate'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'assemRate'],'-v7.3')
    print(h1,[figPath filesep 'assemRate'],'-dtiff')

    tableAssemRateGroup=table(assemRateGroup,'RowNames',nameList);
    writetable(tableAssemRateGroup,[dataPath filesep 'assemRate.csv'],'WriteRowNames',true)
%% disassemRateGroup
    disassemRateGroupCell = cellfun(@(x) cell2mat(x),disassemRateGroup,'unif',false);
    h1=figure; 
    boxPlotCellArray(disassemRateGroupCell,nameList,1,false,true);
    ylabel(['disassemRate (1)'])
    title(['disassembly rate'])
    hgexport(h1,[figPath filesep 'disassemRate'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'disassemRate'],'-v7.3')
    print(h1,[figPath filesep 'disassemRate'],'-dtiff')

    tableDisassemRateGroup=table(disassemRateGroupCell,'RowNames',nameList);
    writetable(tableDisassemRateGroup,[dataPath filesep 'disassemRate.csv'],'WriteRowNames',true)
%% Plotting each - peakGroup - all classes - usually not interesting: there is no lag.
if ~isempty(initRiseProc)
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
        writetable(tablePeakLagGroupEach,[dataPath filesep 'peakLagG' num2str(curGroup) '.csv'],'WriteRowNames',true)
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
        writetable(tableEndLagGroupEach,[dataPath filesep 'endLagG' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
