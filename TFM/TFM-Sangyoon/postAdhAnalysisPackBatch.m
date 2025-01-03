%% open necessary MLs
[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,~,MLdirect]=chooseSelectedFolders;
% Asking user
disp('The current names are: ')
nameList = groupNames';
disp(nameList)
namesOK = input('Do you want to rename your condition names? (y/n)','s');
if strcmp(namesOK, 'y')
    for ii=1:numel(nameList)
        curName = input(['For ' nameList{ii} ': '], 's');
        if ~isempty(curName)
            nameList{ii} = curName;
        end
    end
    specificName = strjoin(nameList);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
% rootAnalysis = pathAnalysisAll{1};
summaryPath = [rootAnalysis '/AnalysisSummary_AdhesionDynamic' specificName];
ii=0;
while exist(summaryPath, 'dir')
    ii=ii+1;
    summaryPath = [rootAnalysis '/AnalysisSummary_AdhesionDynamic' specificName num2str(ii)];
end
figPath = [summaryPath '/Figs'];
mkdir(figPath)
dataPath = [summaryPath '/Data'];
mkdir(dataPath)
save([rootAnalysis filesep 'selectedFolders' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')

%% Loading
numConditions=numel(pathAnalysisAll);

for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end
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
fractionForceTransmittingGroup = cell(numConditions,1);
numEachGroupForceTransmittingGroup = cell(numConditions,1);
assemRateEachGroup = cell(numConditions,1);
disassemRateEachGroup = cell(numConditions,1);

assemRate2EachGroup = cell(numConditions,1);
disassemRate2EachGroup = cell(numConditions,1);

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
assemRate2Group = cell(numConditions,1);
disassemRate2Group = cell(numConditions,1);
nucleatingNARatioGroup = cell(numConditions,1);
disassemNARatioGroup = cell(numConditions,1);

ampEmergingArrayGroup = cell(numConditions,1);
ampFailingArrayGroup = cell(numConditions,1);
ampMaturingArrayGroup = cell(numConditions,1);

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

amplitudeInNAsGroup = cell(numConditions,1);
amplitudeInFCsGroup = cell(numConditions,1);
amplitudeInFAsGroup = cell(numConditions,1);

mainBccPeakValuesGroupGroup = cell(numConditions,1);
mainTimeToPeakGroupGroup = cell(numConditions,1);
sideBccPeakValuesGrouppGroup = cell(numConditions,1);
sideTimeToPeakGroupGroup = cell(numConditions,1);

% iAdhChan = input('Your adhesion channel of interest? (default: 2): ');
% if isempty(iAdhChan); iAdhChan=2; end

forceMeasured=true;

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
    meanLifeTimeAll=cell(N(ii),1);
    meanLifeTimeFailingNAs=cell(N(ii),1);
    meanLifeTimeMaturingNAs=cell(N(ii),1);
    meanMaturingRatio=zeros(N(ii),1);
    meanMaturingRatioNAtoFA=NaN(N(ii),1);
    meanMaturingRatioFCtoFA=NaN(N(ii),1);
    lifeTimeFAsAll=cell(N(ii),1);
    meanStableNAFCratio=NaN(N(ii),1);
    meanAssemRate=cell(N(ii),1);
    meanDisassemRate=cell(N(ii),1);
    meanNucleatingNARatio=cell(N(ii),1);
    meanDisassemNARatio=cell(N(ii),1);
    
    ampMaturingArray=cell(N(ii),1);
    ampFailingArray=cell(N(ii),1);
    ampEmergingArray=cell(N(ii),1);
    
    numG1=NaN(N(ii),1);
    numG2=NaN(N(ii),1);
    numG3=NaN(N(ii),1);
    numG4=NaN(N(ii),1);
    numG5=NaN(N(ii),1);
    numG6=NaN(N(ii),1);
    numG7=NaN(N(ii),1);
    numG8=NaN(N(ii),1);
    numG9=NaN(N(ii),1);
    
    meanRelativePopG1=NaN(N(ii),1);
    meanRelativePopG2=NaN(N(ii),1);
    meanRelativePopG3=NaN(N(ii),1);
    meanRelativePopG4=NaN(N(ii),1);
    meanRelativePopG5=NaN(N(ii),1);
    meanRelativePopG6=NaN(N(ii),1);
    meanRelativePopG7=NaN(N(ii),1);
    meanRelativePopG8=NaN(N(ii),1);
    meanRelativePopG9=NaN(N(ii),1);
    
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
    
    amplitudeInNAs=cell(N(ii),1);   
    amplitudeInFCs=cell(N(ii),1);   
    amplitudeInFAs=cell(N(ii),1);   
    
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
        iAdhChan = find(curAnalProc.checkChannelOutput);

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
        if size(nafaStruct.NADensity,2)==1
            naDensityAllFrames{k} = nafaStruct.NADensity;
    %         meanNADensity(k) = mean(curNADensity);
            faDensityAllFrames{k} = nafaStruct.FADensity;
    %         meanFADensity(k) = meanNumPureFAs(k)/meanCellArea(k); % in num/um2
        else
            naDensityAllFrames{k} = nafaStruct.NADensity';
            faDensityAllFrames{k} = nafaStruct.FADensity';
        end

        maturingStruct=load(curAnalProc.outFilePaths_{4,iAdhChan});
        meanLifeTimeAll{k}=(maturingStruct.lifeTimeAll);
        meanLifeTimeFailingNAs{k}=(maturingStruct.lifeTimeNAfailing);
        meanLifeTimeMaturingNAs{k}=(maturingStruct.lifeTimeNAmaturing);
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
        meanAssemRate2{k}=(assemRateStruct.assemRateCell);
        meanDisassemRate2{k}=(assemRateStruct.disassemRateCell);
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
        numG6(k) = sum(idsStruct.idGroup6);
        
        try
            curPresenceStruct = load(classProc.outFilePaths_{5,iAdhChan},'presenceAll');
            curPresence = curPresenceStruct.presenceAll;
        catch % This is for backward compatibility - SH 20190403
            tracksNA = curAnalProc.loadChannelOutput(iAdhChan,'output','tracksNA'); %load(classProc.outFilePaths_{5,iAdhChan},'tableTracksNA');
            tableTracksNA = struct2table(tracksNA);
            curPresence = tableTracksNA.presence; presenceAll=curPresence;
            try
                movefile(classProc.outFilePaths_{5,iAdhChan},[classProc.outFilePaths_{5,iAdhChan}(1:end-4) 'Original.mat'])
            catch
                disp(['Overwriting ' classProc.outFilePaths_{5,iAdhChan} ' ...'])
            end
            save(classProc.outFilePaths_{5,iAdhChan},'presenceAll')
            clear tracksNA tableTracksNA
        end
        
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
            theOtherFunParams = tOtherProc.funParams_;
            if ~isempty(tOtherProc)
                iOther = theOtherFunParams.iChanSlave; %input('Another channel other than the main adhesion channel? (default: 3): ');
                if isempty(iOther); iOther=3; end

                intenGroupStruct=load(tOtherProc.outFilePaths_{2,iOther});
                intenGroup=intenGroupStruct.intensityGroup; 
                intensitiesInNAs{k} = intenGroup{1};
                intensitiesInFCs{k} = intenGroup{2};
                intensitiesInFAs{k} = intenGroup{3};

                try
                    amp2Group=intenGroupStruct.amplitudeGroup; 
                    amplitudeInNAs{k} = amp2Group{1};
                    amplitudeInFCs{k} = amp2Group{2};
                    amplitudeInFAs{k} = amp2Group{3};
                catch
                    amplitudeInNAs{k} = [];
                    amplitudeInFCs{k} = [];
                    amplitudeInFAs{k} = [];
                end                    
            else
                intensitiesInNAs{k} = [];
                intensitiesInFCs{k} = [];
                intensitiesInFAs{k} = [];

                amplitudeInNAs{k} = [];
                amplitudeInFCs{k} = [];
                amplitudeInFAs{k} = [];
            end        
            % Other feature related properties are calculated in the step 11
        end        
    end
    
    
    initRiseAgainstForceEachClass=cell(numClasses,1);
    peakTimeAgainstForceEachClass=cell(numClasses,1);
    endTimeAgainstForceEachClass=cell(numClasses,1);
    halfBccTogetherEachClass=cell(numClasses,1);
    
    earlyAmpSlopeEachClass=cell(numClasses,1);
    earlyForceSlopeEachClass=cell(numClasses,1);
    assemRateEachClass=cell(numClasses,1);
    disassemRateEachClass=cell(numClasses,1);
    
    fractionForceTransmittingEachClass=cell(numClasses,1);
    numEachGroupForceTransmitting=cell(numClasses,1);
    
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
                try
                    initOutFolder = fileparts(initRiseProc.outFilePaths_{2,iForceSlave});
                catch
                    try
                        initOutFolder = fileparts(initRiseProc.outFilePaths_{1,iForceSlave});
                    catch
                        % get whatever that's not empty
                        iFilledOutput = find(~cellfun(@isempty,initRiseProc.outFilePaths_));
                        if ~isempty(iFilledOutput)
                            initOutFolder = fileparts(initRiseProc.outFilePaths_{iFilledOutput(1)});
                        else
                            disp('No init rise process results! Skipping ...')
                            continue
                        end
                    end
                end
                initDataPath = [initOutFolder filesep 'data'];

                nameTitle=['initialLag Class' num2str(pp)]; curClassPath=[initDataPath filesep nameTitle '.mat'];
                if exist(curClassPath,'file')
                    initRiseStruct = load(curClassPath,'initialLagTogetherAdjusted','nameList2');   
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

                    nameTitle=['halfBccTogetherAdjusted' num2str(pp)];
                    halfBccTogetherStruct = load([initDataPath filesep nameTitle],'halfBccTogetherAdjusted','nameList2');    
                    halfBccTogetherEachClass{pp}{k} = halfBccTogetherStruct.halfBccTogetherAdjusted;
                    
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
                    
                    %Force related
                    try
                        earlyForceSlopeStr=load([initDataPath filesep 'earlyForceSlopeAllGroups.mat'],'earlyForceSlope');
                        earlyForceSlopeEachClass{pp}{k} = earlyForceSlopeStr.earlyForceSlope{pp};

                        fractionForceTransmittingStr=load([initDataPath filesep 'forceTransmitting.mat'],'fractionForceTransmitting','numEachGroup','numEachGroupForceTransmitting','groupLabel');
                        fractionForceTransmittingEachClass{pp}(k) = fractionForceTransmittingStr.fractionForceTransmitting(pp);               
                        numEachGroupForceTransmitting{pp}(k) = fractionForceTransmittingStr.numEachGroup(pp);               
                    catch
                        disp('No force measured')
                        forceMeasured=false;
                    end
                end
                
                % AssemRate and DisassemRate related
                earlyAmpSlopPath=[initDataPath filesep 'earlyAmpSlopeAllGroups.mat'];
                if exist(earlyAmpSlopPath,'file')
                    earlyAmpSlopeStr=load([initDataPath filesep 'earlyAmpSlopeAllGroups.mat'],'earlyAmpSlope');
                    earlyAmpSlopeEachClass{pp}{k} = earlyAmpSlopeStr.earlyAmpSlope{pp};
                end
                if exist([initDataPath filesep 'assemRate.mat'],'file')
                    assemRateStr=load([initDataPath filesep 'assemRate.mat'],'assemRate');
                    assemRateEachClass{pp}{k} = assemRateStr.assemRate{pp};
                    disassemRateStr=load([initDataPath filesep 'disassemRate.mat'],'disassemRate');
                    disassemRateEachClass{pp}{k} = disassemRateStr.disassemRate{pp};                    
                else
                    assemRateEachClass{pp}{k} = [];
                    disassemRateEachClass{pp}{k} = [];                    
                end
                if exist([initDataPath filesep 'assemRate2.mat'],'file')
                    if pp==1
                        assemRate2Str=load([initDataPath filesep 'assemRate2.mat'],'assemRate2');
                        disassemRate2Str=load([initDataPath filesep 'disassemRate2.mat'],'disassemRate2');
                        meanAssemRate2{k} = cell2mat(assemRate2Str.assemRate2');
                        meanDisassemRate2{k} = cell2mat(disassemRate2Str.disassemRate2');
                    end
                    assemRate2EachClass{pp}{k} = assemRate2Str.assemRate2{pp};
                    disassemRate2EachClass{pp}{k} = disassemRate2Str.disassemRate2{pp};                    
                else
                    assemRate2EachClass{pp}{k} = [];
                    disassemRate2EachClass{pp}{k} = [];                    
                end
                % Amp profile related
                if pp==1 % this needs only once
                    ampArrayPath=[initDataPath filesep 'ampArrays.mat'];
                    if exist(ampArrayPath,'file')
                        ampArrayStr=load(ampArrayPath);
                        ampEmergingArray{k}=ampArrayStr.ampEmergingArray;
                        ampFailingArray{k}=ampArrayStr.ampFailingArray;
                        ampMaturingArray{k}=ampArrayStr.ampMaturingArray;
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
    if exist([initDataPath filesep 'assemRate2.mat'],'file')
        assemRate2Group{ii,1}=meanAssemRate2;
        disassemRate2Group{ii,1}=meanDisassemRate2;
    end
    nucleatingNARatioGroup{ii,1}=meanNucleatingNARatio;
    disassemNARatioGroup{ii,1}=meanDisassemNARatio;
    
    ampEmergingArrayGroup{ii,1} = ampEmergingArray;
    ampFailingArrayGroup{ii,1} = ampFailingArray;
    ampMaturingArrayGroup{ii,1} = ampMaturingArray;
    
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
    
    amplitudeInNAsGroup{ii,1}=amplitudeInNAs;
    amplitudeInFCsGroup{ii,1}=amplitudeInFCs;
    amplitudeInFAsGroup{ii,1}=amplitudeInFAs;

    if ~isempty(initRiseProc)
        initRiseGroup{ii,1}=initRiseAgainstForceEachClass;
        peakGroup{ii,1}=peakTimeAgainstForceEachClass;
        endTimeGroup{ii,1}=endTimeAgainstForceEachClass;

        halfBccGroup{ii,1}=halfBccTogetherEachClass;

        earlyAmpSlopeGroup{ii,1}=earlyAmpSlopeEachClass;
        earlyForceSlopeGroup{ii,1}=earlyForceSlopeEachClass;
        assemRateEachGroup{ii,1}=assemRateEachClass;
        disassemRateEachGroup{ii,1}=disassemRateEachClass;
        if exist([initDataPath filesep 'assemRate2.mat'],'file')
            assemRate2EachGroup{ii,1}=assemRate2EachClass;
            disassemRate2EachGroup{ii,1}=disassemRate2EachClass;
        end
        fractionForceTransmittingGroup{ii,1} = fractionForceTransmittingEachClass;
        numEachGroupForceTransmittingGroup{ii,1} = numEachGroupForceTransmitting;

        mainBccPeakValuesGroupGroup{ii,1}=mainBccPeakValuesGroup;
        mainTimeToPeakGroupGroup{ii,1}=mainTimeToPeakGroup;
        sideBccPeakValuesGrouppGroup{ii,1}=sideBccPeakValuesGroup;
        sideTimeToPeakGroupGroup{ii,1}=sideTimeToPeakGroup;
    end
end
disp('Done')
%% Plotting each - initRiseGroup against force - all classes
if ~isempty(initRiseProc) && exist('initRiseStruct','var')
    numSlaves = numel(initRiseStruct.nameList2);
    for ii=1:numSlaves
        for curGroup=1:9
            try
                initRiseGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) y{ii},x{curGroup}(~cellfun(@isempty,x{curGroup}))','unif',false)),initRiseGroup,'unif',false);
                h1=figure; 
                boxPlotCellArray(initRiseGroupEach,nameList,1,false,true);
                ylabel(['Initial rise in group ' num2str(curGroup) ' (sec): ' initRiseStruct.nameList2{ii}])
                title(['Initial rise in group ' num2str(curGroup) ': ' initRiseStruct.nameList2{ii}])
                hgexport(h1,[figPath filesep 'initialRiseG' num2str(curGroup) initRiseStruct.nameList2{ii}],hgexport('factorystyle'),'Format','eps')
                hgsave(h1,[figPath filesep 'initialRiseG' num2str(curGroup) initRiseStruct.nameList2{ii}],'-v7.3')
                print(h1,[figPath filesep 'initialRiseG' num2str(curGroup) initRiseStruct.nameList2{ii}],'-dtiff')

                tableInitialRiseGroup=table(initRiseGroupEach,'RowNames',nameList);
                writetable(tableInitialRiseGroup,[dataPath filesep 'initialRiseG' num2str(curGroup) initRiseStruct.nameList2{ii} '.csv'],'WriteRowNames',true)
            catch
                disp(['Nothing in ' num2str(curGroup) 'th group'])
            end
        end
    end
end
%% halfBccTime
if ~isempty(initRiseProc) && exist('halfBccGroup','var') && ~isempty(halfBccGroup{1}{1})
    for curGroup=1:9
        try        
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
        catch
            disp(['Nothing in ' num2str(curGroup) 'th group'])
        end
    end
end
%% earlyAmpSlope
if ~isempty(initRiseProc) && exist('earlyAmpSlopeGroup','var') && ~isempty(earlyAmpSlopeGroup{1}{1})
    for curGroup=1:9
        try
            earlyAmpSlopeGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),earlyAmpSlopeGroup,'unif',false);
            h1=figure; 
            boxPlotCellArray(earlyAmpSlopeGroupEach,nameList,1,false,true);
            ylabel(['earlyAmpSlope ' num2str(curGroup) ' (a.u./min)'])
            title(['earlyAmpSlope in group ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'earlyAmpSlope' num2str(curGroup)],'-dtiff')

            tableEarlyAmpSlope=table(earlyAmpSlopeGroupEach,'RowNames',nameList);
            writetable(tableEarlyAmpSlope,[dataPath filesep 'earlyAmpSlope' num2str(curGroup) '.csv'],'WriteRowNames',true)
        catch
            disp(['Nothing in ' num2str(curGroup) 'th group'])
        end
    end
end
%% earlyAmpSlope - G1G2 only
if ~isempty(initRiseProc) && exist('earlyAmpSlopeGroup','var') && ~isempty(earlyAmpSlopeGroup{1}{1})
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
    ylabel('earlyAmpSlope (a.u./min)')
    title('earlyAmpSlope in G1 vs. G2')
    hgexport(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],'-v7.3')
    print(h1,[figPath filesep 'earlyAmpSlopeG1G2' num2str(curGroup)],'-dtiff')
end            
%% earlyForceSlope
if ~isempty(initRiseProc) && exist('earlyForceSlopeGroup','var') && ~isempty(earlyForceSlopeGroup{1}{1})
    for curGroup=1:9
        earlyForceSlopeGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),earlyForceSlopeGroup,'unif',false);
        h1=figure; 
        plotSuccess=boxPlotCellArray(earlyForceSlopeGroupEach,nameList,1,false,true);
        if plotSuccess
            ylabel(['earlyForceSlope ' num2str(curGroup) ' (sec)'])
            title(['earlyForceSlope in group ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'earlyForceSlope' num2str(curGroup)],'-dtiff')

            tableEarlyForceSlope=table(earlyForceSlopeGroupEach,'RowNames',nameList);
            writetable(tableEarlyForceSlope,[dataPath filesep 'earlyForceSlope' num2str(curGroup) '.csv'],'WriteRowNames',true)
        end
    end
end
%% earlyForceSlope - G1G2 only
if ~isempty(initRiseProc) && exist('earlyForceSlopeGroup','var') && ~isempty(earlyForceSlopeGroup{1}{1})
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
    plotSuccess=boxPlotCellArray(earlyForceSlopeG1G2,nameListG1G2,1,false,true);
    if plotSuccess
        ylabel('earlyForceSlope (a.u./sec)')
        title('earlyForceSlope in G1 vs. G2')
        hgexport(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'earlyForceSlopeG1G2' num2str(curGroup)],'-dtiff')
    end
end  
%% fraction of force transmitting
if ~isempty(initRiseProc) && ~isempty(curFAPackage.getProcess(10))
    for ii=1:numConditions
        curFractionForceTransmitting = fractionForceTransmittingGroup{ii,1};
        curNumEachGroupForceTransmitting = numEachGroupForceTransmittingGroup{ii,1};
        try
            groupLabel = arrayfun(@(x,y) ['G' num2str(x) '(M=' num2str(numel(y{1})) ',N=' num2str(sum(y{1})) ')'],1:9,curNumEachGroupForceTransmitting,'unif',false);
        catch
            groupLabel = arrayfun(@(x,y) ['G' num2str(x) '(M=' num2str(numel(y{1})) ',N=' num2str(sum(y{1})) ')'],1:9,curNumEachGroupForceTransmitting','unif',false);
        end
        
        h1=figure; 
        plotSuccess=boxPlotCellArray(curFractionForceTransmitting,groupLabel',1,false,true);
        if plotSuccess
            ylabel('Fraction of force transmitting NAs (1)')
            title(['Fraction of force transmitting NAs. ' groupNames{ii}])
            hgexport(h1,[figPath filesep 'fractionForceTransmitting_' groupNames{ii}],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'fractionForceTransmitting_' groupNames{ii}],'-v7.3')
            print(h1,[figPath filesep 'fractionForceTransmitting_' groupNames{ii}],'-dtiff')
        end
    end
end
%% Plotting each - mainTimeToPeakGroupGroup - G1 and G2
if ~isempty(initRiseProc) && exist('mainTimeToPeakGroupGroup','var') && ~isempty(mainTimeToPeakGroupGroup{1}{1})
    for curGroup=1:2
        mainTimeToPeakGroupGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),mainTimeToPeakGroupGroup,'unif',false);
        h1=figure; 
        plotSuccess=boxPlotCellArray(mainTimeToPeakGroupGroupEach,nameList,1,false,true);
        if plotSuccess
            ylabel(['mainTimeToPeakGroupGroupEach ' num2str(curGroup) ' (sec)'])
            title(['mainTimeToPeakGroupGroupEach ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'mainTimeToPeakGroup' num2str(curGroup)],'-dtiff')
        end
    end
    %% Plotting each - mainBccPeakValuesGroupGroup - G1 and G2
    for curGroup=1:2
        mainBccPeakValuesGroupGrouppEach = cellfun(@(x) cell2mat(x{curGroup}'),mainBccPeakValuesGroupGroup,'unif',false);
        h1=figure; 
        plotSuccess=boxPlotCellArray(mainBccPeakValuesGroupGrouppEach,nameList,1,false,true);
        if plotSuccess
            ylabel(['mainBccPeakValuesGroupGrouppEach ' num2str(curGroup) ' (sec)'])
            title(['mainBccPeakValuesGroupGrouppEach ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'mainBccPeakValuesGroupGrouppEach' num2str(curGroup)],'-dtiff')
        end
    end
    %% Plotting each - sideTimeToPeakGroupGroup - G1 and G2
    for curGroup=1:2
        sideTimeToPeakGroupGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),sideTimeToPeakGroupGroup,'unif',false);
        h1=figure; 
        plotSuccess=boxPlotCellArray(sideTimeToPeakGroupGroupEach,nameList,1,false,true);
        if plotSuccess
            ylabel(['sideTimeToPeakGroupGroup ' num2str(curGroup) ' (sec)'])
            title(['sideTimeToPeakGroupGroup ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'sideTimeToPeakGroupGroup' num2str(curGroup)],'-dtiff')
        end
    end
    %% Plotting each - sideBccPeakValuesGrouppGroup - G1 and G2
    for curGroup=1:2
        sideBccPeakValuesGrouppGrouppEach = cellfun(@(x) cell2mat(x{curGroup}'),sideBccPeakValuesGrouppGroup,'unif',false);
        h1=figure; 
        plotSuccess=boxPlotCellArray(sideBccPeakValuesGrouppGrouppEach,nameList,1,false,true);
        if plotSuccess
            ylabel(['sideBccPeakValuesGrouppGroup ' num2str(curGroup) ' (sec)'])
            title(['sideBccPeakValuesGrouppGroup ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'sideBccPeakValuesGrouppGroup' num2str(curGroup)],'-dtiff')
        end
    end
end
%% num of each class
numGroupOnlyTwo={numG1Group,numG2Group}; %,numG3Group,numG4Group,numG5Group,numG6Group,numG7Group,numG8Group,numG9Group};
for curGroup=1:2
    curNumGroup=numGroupOnlyTwo{curGroup};
    h1=figure; 
    plotSuccess=boxPlotCellArray(curNumGroup,nameList,1,false,true);
    if plotSuccess
        ylabel(['Number of adhesions in group ' num2str(curGroup) ' (#)'])
        title(['Number of adhesions in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'numGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'numGroup' num2str(curGroup)],'-dtiff')
        tableAdhNumGroupEach=table(curNumGroup,'RowNames',nameList);
        writetable(tableAdhNumGroupEach,[dataPath filesep 'numGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end
%% num of adhesion for all classes - new nameList
nameListAllTogether=cell(1,numConditions*9);
for ii=1:numConditions
    for jj=1:9
        nameListAllTogether{1,(jj-1)*numConditions+ii}= [nameList{ii} 'G' num2str(jj)];
    end
end
%% num of adhesion for all classes
numGroupAll={numG1Group,numG2Group,numG3Group,numG4Group,numG5Group,numG6Group,numG7Group,numG8Group,numG9Group};
numGroupAllLinear = cell(1,numConditions*9);
for ii=1:numConditions
    for jj=1:9
        numGroupAllLinear{1,(jj-1)*2+ii}=numGroupAll{jj}{ii,1};
    end
end
h1=figure; 
% barPlotCellArray(numGroupAllLinear,nameListAllTogether,1);
plotSuccess=boxPlotCellArray(numGroupAllLinear,nameListAllTogether,1,false, true);
if plotSuccess
    ylabel(['Number of adhesions in each group (#)'])
    title(['Number of adhesions in each group '])
    hgexport(h1,[figPath filesep 'numGroupAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numGroupAll'],'-v7.3')
    print(h1,[figPath filesep 'numGroupAll'],'-dtiff')
    tablenumGroupAll=table(numGroupAllLinear','RowNames',nameListAllTogether');
    writetable(tablenumGroupAll,[dataPath filesep 'numGroupAll.csv'],'WriteRowNames',true)
end

%% meanRelative population
meanRelPopGroupAll={meanRelativePopG1Group,meanRelativePopG2Group,meanRelativePopG3Group,...
    meanRelativePopG4Group,meanRelativePopG5Group,meanRelativePopG6Group,...
    meanRelativePopG7Group,meanRelativePopG8Group,meanRelativePopG9Group};
for curGroup=1:9
    curmeanPopGroup=meanRelPopGroupAll{curGroup};
    h1=figure; 
    plotSuccess=boxPlotCellArray(curmeanPopGroup,nameList,1,false,true);
    if plotSuccess
        ylabel(['Mean relative adhesion population in group ' num2str(curGroup) ' (1)'])
        title(['Mean relative adhesion population in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'meanRelPopGroup' num2str(curGroup)],'-dtiff')

        tableMeanPopGroupEach=table(curmeanPopGroup,'RowNames',nameList);
        writetable(tableMeanPopGroupEach,[dataPath filesep 'meanRelPopGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end

%% relative population for all classes
meanRelPopGroupAllLinear = cell(1,numConditions*9);
for ii=1:numConditions
    for jj=1:9
        meanRelPopGroupAllLinear{1,(jj-1)*2+ii}=meanRelPopGroupAll{jj}{ii,1};
    end
end
h1=figure; 
% barPlotCellArray(numGroupAllLinear,nameListAllTogether,1);
plotSuccess=boxPlotCellArray(meanRelPopGroupAllLinear,nameListAllTogether,1,false, true);
if plotSuccess
    ylabel(['Relative adhesion population in each group (#/#)'])
    title(['Relative adhesion population in each group '])
    hgexport(h1,[figPath filesep 'meanRelPopAllGroup'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanRelPopAllGroup'],'-v7.3')
    print(h1,[figPath filesep 'meanRelPopAllGroup'],'-dtiff')
    tableAll=table(meanRelPopGroupAllLinear','RowNames',nameListAllTogether');
    writetable(tableAll,[dataPath filesep 'meanRelPopAllGroup.csv'],'WriteRowNames',true)
end
%% RAtio between g1 and g2
ratioG1G2=cellfun(@(x,y) y./x,meanRelativePopG1Group,meanRelativePopG2Group,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(ratioG1G2,nameList,1,false,true);
if plotSuccess
    ylabel(['Mean relative adhesion population of g2 compared to g1'])
    title('Mean relative adhesion population of g2 compared to g1')
    hgexport(h1,[figPath filesep 'meanRelPopG2overG1'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanRelPopG2overG1'],'-v7.3')
    print(h1,[figPath filesep 'meanRelPopG2overG1'],'-dtiff')

    tableRelPopG2overG1=table(ratioG1G2,'RowNames',nameList);
    writetable(tableRelPopG2overG1,[dataPath filesep 'meanRelPopG2overG1' num2str(curGroup) '.csv'],'WriteRowNames',true)
end
%% mean adhesion density - significantly higher NA density in WT compared to R8
meanAdhDensityGroupAll={meanAdhDensityG1Group,meanAdhDensityG2Group,meanAdhDensityG3Group,...
    meanAdhDensityG4Group,meanAdhDensityG5Group,meanAdhDensityG6Group,...
    meanAdhDensityG7Group,meanAdhDensityG8Group,meanAdhDensityG9Group};
for curGroup=1:9
    curAdhDenGroup=meanAdhDensityGroupAll{curGroup};
    h1=figure; 
    curAdhDenGroupCell = cellfun(@(x) cell2mat(x'),curAdhDenGroup,'unif',false);

    plotSuccess=boxPlotCellArray(curAdhDenGroupCell,nameList,1,false,true);
    if plotSuccess
        ylabel(['Mean adhesion density in group ' num2str(curGroup) ' (num/um^2)'])
        title(['Mean adhesion density in group ' num2str(curGroup)])
        hgexport(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-v7.3')
        print(h1,[figPath filesep 'meanAdhDenGroup' num2str(curGroup)],'-dtiff')

        tableAdhDenGroupEach=table(curAdhDenGroup,'RowNames',nameList);
        writetable(tableAdhDenGroupEach,[dataPath filesep 'meanAdhDenGroup' num2str(curGroup) '.csv'],'WriteRowNames',true)
    end
end

%% mean adhesion density for all classes
meanAdhDensityGroupAllLinear = cell(1,numConditions*9);
for ii=1:numConditions
    for jj=1:9
        meanAdhDensityGroupAllLinear{1,(jj-1)*2+ii}=cell2mat(meanAdhDensityGroupAll{jj}{ii,1}');
    end
end
h1=figure; 
% barPlotCellArray(numGroupAllLinear,nameListAllTogether,1);
plotSuccess=boxPlotCellArray(meanAdhDensityGroupAllLinear,nameListAllTogether,1,false, true);
if plotSuccess
    ylabel(['Mean adhesion density in each group (#/#)'])
    title(['Mean adhesion density in each group '])
    hgexport(h1,[figPath filesep 'meanAdhDenAllGroup'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'meanAdhDenAllGroup'],'-v7.3')
    print(h1,[figPath filesep 'meanAdhDenAllGroup'],'-dtiff')
    tableAll=table(meanAdhDensityGroupAllLinear','RowNames',nameListAllTogether');
    writetable(tableAll,[dataPath filesep 'meanAdhDenAllGroup.csv'],'WriteRowNames',true)
end
%% cell area - 
convertL=curMovie.pixelSize_/1000; %um/pix
convertArea = convertL^2;

h1=figure; 
plotSuccess=boxPlotCellArray(cellAreaGroup,nameList,convertArea,false,true);
if plotSuccess
    ylabel(['Cell spread area (um^2)'])
    title(['Cell spread area'])
    hgexport(h1,[figPath filesep 'cellArea'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'cellArea'],'-v7.3')
    print(h1,[figPath filesep 'cellArea'],'-dtiff')

    tableCellArea=table(cellAreaGroup,'RowNames',nameList);
    writetable(tableCellArea,[dataPath filesep 'cellArea.csv'],'WriteRowNames',true)
end
%% FA area - 
FAareaGroupCell = cellfun(@(x) cell2mat(x),FAareaGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(FAareaGroupCell,nameList,convertArea,false,true);
if plotSuccess
    ylabel(['Focal adhesion area (\mum^2)'])
    title(['Focal adhesion area'])
    hgexport(h1,[figPath filesep 'faArea'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faArea'],'-v7.3')
    print(h1,[figPath filesep 'faArea'],'-dtiff')

    tableFAArea=table(FAareaGroupCell,'RowNames',nameList);
    writetable(tableFAArea,[dataPath filesep 'faArea.csv'],'WriteRowNames',true)
end
%% Percentage of FA area
percLT=5;
perc=percLT/100;
FAareaGroupCellSmall=cell(numel(FAareaGroupCell),1);
for ii=1:numel(FAareaGroupCell)
    FAareaGroupCellSmall{ii,1} = ...
    quantile(FAareaGroupCell{ii},(1-perc)+(perc-0.01)*rand(1,round((perc-0.01)*sum(~isnan(FAareaGroupCell{ii})))));
end

h1=figure;
plotSuccess=boxPlotCellArray(FAareaGroupCellSmall,nameList,convertArea,false,true);
if plotSuccess
    ylabel(['Focal adhesion area (um^2)'])
    title(['Focal adhesion area (top ' num2str(percLT) ' percentile)'])
    ylim auto
    hgexport(h1,[figPath filesep 'faAreaTop' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faAreaTop' num2str(percLT)],'-v7.3')
    print(h1,[figPath filesep 'faAreaTop' num2str(percLT)],'-dtiff')
    faAreaTop20=table(FAareaGroupCellSmall,'RowNames',nameList);
    writetable(faAreaTop20,[dataPath filesep 'faAreaTop' num2str(percLT) '.csv'],'WriteRowNames',true)    
end
%% FA length - 
FAlengthGroupCell = cellfun(@(x) cell2mat(x),FAlengthGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(FAlengthGroupCell,nameList,convertL,false,true);
if plotSuccess
    ylabel(['Focal adhesion length (\mum)'])
    title(['Focal adhesion length'])
    hgexport(h1,[figPath filesep 'faLength'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faLength'],'-v7.3')
    print(h1,[figPath filesep 'faLength'],'-dtiff')

    tableFALength=table(FAlengthGroupCell,'RowNames',nameList);
    writetable(tableFALength,[dataPath filesep 'faLength.csv'],'WriteRowNames',true)
end
%% FA length - bar plot
h1=figure;
barPlotCellArray(FAlengthGroupCell,nameList,convertL);
ylabel(['Focal adhesion length (\mum)'])
title(['Focal adhesion length - bar'])
hgexport(h1,[figPath filesep 'faLengthBar'],hgexport('factorystyle'),'Format','eps')
hgsave(h1,[figPath filesep 'faLengthBar'],'-v7.3')
print(h1,[figPath filesep 'faLengthBar'],'-dtiff')

%% NA density
NADensityGroupCell = cellfun(@(x) cell2mat(x),NADensityGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(NADensityGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel(['NA density (#/\mum^2)'])
    title(['NA density'])
    hgexport(h1,[figPath filesep 'naDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'naDensity'],'-v7.3')
    print(h1,[figPath filesep 'naDensity'],'-dtiff')

    tableNADensity=table(NADensityGroup,'RowNames',nameList);
    writetable(tableNADensity,[dataPath filesep 'naDensity.csv'],'WriteRowNames',true)
end
%% FA density
FADensityGroupCell = cellfun(@(x) cell2mat(x),FADensityGroup,'unif',false);
h1=figure; 
%     FADensityGroup=cellfun(@(x,y) x./y,numPureFAsGroup,cellAreaGroup,'unif',false);
plotSuccess=boxPlotCellArray(FADensityGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel(['FA density (#/um^2)'])
    title(['FA density'])
    hgexport(h1,[figPath filesep 'faDensity'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'faDensity'],'-v7.3')
    print(h1,[figPath filesep 'faDensity'],'-dtiff')

    tableFADensity=table(FADensityGroupCell,'RowNames',nameList);
    writetable(tableFADensity,[dataPath filesep 'faDensity.csv'],'WriteRowNames',true)
end
%% Life times
lifeTimeAllGroupCell = cellfun(@(x) cell2mat(x'),lifeTimeAllGroup,'unif',false);
h1=figure; 
%     FADensityGroup=cellfun(@(x,y) x./y,numPureFAsGroup,cellAreaGroup,'unif',false);
plotSuccess=boxPlotCellArray(lifeTimeAllGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel(['Life time (sec)'])
    title(['Life time of all NAs'])
    hgexport(h1,[figPath filesep 'lifeTimeAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeAll'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeAll'],'-dtiff')

    lifeTimeAll=table(lifeTimeAllGroupCell,'RowNames',nameList);
    writetable(lifeTimeAll,[dataPath filesep 'lifeTimeAll.csv'],'WriteRowNames',true)    
end
%% lifeTimeFailingNAsGroup
lifeTimeFailingNAsGroupCell = cellfun(@(x) cell2mat(x'),lifeTimeFailingNAsGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(lifeTimeFailingNAsGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel(['Life time (sec)'])
    title(['Life time of non-maturing NAs'])
    hgexport(h1,[figPath filesep 'lifeTimeFailing'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeFailing'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeFailing'],'-dtiff')

    lifeTimeFailing=table(lifeTimeFailingNAsGroupCell,'RowNames',nameList);
    writetable(lifeTimeFailing,[dataPath filesep 'lifeTimeFailing.csv'],'WriteRowNames',true)    
end
%% lifeTimeMaturingNAsGroup
lifeTimeMaturingNAsGrouppCell = cellfun(@(x) cell2mat(x'),lifeTimeMaturingNAsGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(lifeTimeMaturingNAsGrouppCell,nameList,1,false,true);
if plotSuccess
    ylabel(['Life time (sec)'])
    title(['Life time of maturing NAs'])
    hgexport(h1,[figPath filesep 'lifeTimeMaturing'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeMaturing'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeMaturing'],'-dtiff')

    lifeTimeMaturing=table(lifeTimeMaturingNAsGrouppCell,'RowNames',nameList);
    writetable(lifeTimeMaturing,[dataPath filesep 'lifeTimeMaturing.csv'],'WriteRowNames',true)    
end

% lifeTimeMaturingNAsGrouppCellSmall{1} = ...
% quantile(lifeTimeMaturingNAsGrouppCell{1},0.8+0.19*rand(1,round(0.1*sum(~isnan(lifeTimeMaturingNAsGrouppCell{1})))));
% lifeTimeMaturingNAsGrouppCellSmall{2,1} = ...
% quantile(lifeTimeMaturingNAsGrouppCell{2},0.8+0.19*rand(1,round(0.1*sum(~isnan(lifeTimeMaturingNAsGrouppCell{2})))));
lifeTimeMaturingNAsGrouppCellSmall=cell(numConditions,1)
for ii=1:numConditions
    lifeTimeMaturingNAsGrouppCellSmall{ii,1} = ...
    quantile(lifeTimeMaturingNAsGrouppCell{2},0.8+0.19*rand(1,round(0.1*sum(~isnan(lifeTimeMaturingNAsGrouppCell{ii})))));
end
h1=figure; 
plotSuccess=boxPlotCellArray(lifeTimeMaturingNAsGrouppCellSmall,nameList,1,false,true);
if plotSuccess
    ylabel(['Life time (sec)'])
    title(['Life time of maturing NAs (top 20 percentile)'])
    ylim auto
    hgexport(h1,[figPath filesep 'lifeTimeMaturingTop20'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeMaturingTop20'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeMaturingTop20'],'-dtiff')
    lifeTimeMaturingT20=table(lifeTimeMaturingNAsGrouppCellSmall,'RowNames',nameList);
    writetable(lifeTimeMaturingT20,[dataPath filesep 'lifeTimeMaturingT20.csv'],'WriteRowNames',true)    
end
    %% numPureFAsGroup
h1=figure; 
plotSuccess=boxPlotCellArray(numPureFAsGroup,nameList,1,false,true);
if plotSuccess
    ylabel(['The number of pure FAs in a cell (#/cell)'])
    title(['The number of pure FAs in a cell'])
    hgexport(h1,[figPath filesep 'numFAs'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'numFAs'],'-v7.3')
    print(h1,[figPath filesep 'numFAs'],'-dtiff')

    tableNumPureFAs=table(numPureFAsGroup,'RowNames',nameList);
    writetable(tableNumPureFAs,[dataPath filesep 'numFAs.csv'],'WriteRowNames',true)
end
%% maturingRatioGroup
h1=figure; 
plotSuccess=boxPlotCellArray(maturingRatioGroup,nameList,1,false,true);
if plotSuccess
    ylabel(['Maturing ratio (1)'])
    title(['Maturing ratio of NAs to FCs'])
    hgexport(h1,[figPath filesep 'maturingRatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatio'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatio'],'-dtiff')

    tableMaturingRatio=table(maturingRatioGroup,'RowNames',nameList);
    writetable(tableMaturingRatio,[dataPath filesep 'maturingRatio.csv'],'WriteRowNames',true)
end
%% maturingRatioNAtoFAGroup
h1=figure; 
plotSuccess=boxPlotCellArray(maturingRatioNAtoFAGroup,nameList,1,false,true);
if plotSuccess
    ylabel(['Maturing ratio from NA to FA (1)'])
    title(['Maturing ratio of NAs to FAs'])
    hgexport(h1,[figPath filesep 'maturingRatioNAtoFA'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatioNAtoFA'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatioNAtoFA'],'-dtiff')

    tableMaturingRatioNAtoFAGroup=table(maturingRatioNAtoFAGroup,'RowNames',nameList);
    writetable(tableMaturingRatioNAtoFAGroup,[dataPath filesep 'maturingRatioNAtoFA.csv'],'WriteRowNames',true)
end
%% maturingRatioFCtoFAGroup
h1=figure; 
plotSuccess=boxPlotCellArray(maturingRatioFCtoFAGroup,nameList,1,false,true);
if plotSuccess
    ylabel('Maturing ratio from FCs to FAs (1)')
    title(['Maturing ratio of FCs to FAs'])
    hgexport(h1,[figPath filesep 'maturingRatioFCtoFA'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'maturingRatioFCtoFA'],'-v7.3')
    print(h1,[figPath filesep 'maturingRatioFCtoFA'],'-dtiff')

    tableMaturingRatioFCtoFA=table(maturingRatioFCtoFAGroup,'RowNames',nameList);
    writetable(tableMaturingRatioFCtoFA,[dataPath filesep 'maturingRatioFCtoFA.csv'],'WriteRowNames',true)
end
%% stableNAFCratioGroup
h1=figure; 
plotSuccess=boxPlotCellArray(stableNAFCratioGroup,nameList,1,false,true);
if plotSuccess
    ylabel('The ratio of stable NAs and FCs (1)')
    title(['The ratio of stable NAs and FCs'])
    hgexport(h1,[figPath filesep 'stableNAFCratio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'stableNAFCratio'],'-v7.3')
    print(h1,[figPath filesep 'stableNAFCratio'],'-dtiff')

    tableStableNAFCratio=table(stableNAFCratioGroup,'RowNames',nameList);
    writetable(tableStableNAFCratio,[dataPath filesep 'stableNAFCratio.csv'],'WriteRowNames',true)
end
%% lifeTimeFAsGroup
lifeTimeFAsGroupCell = cellfun(@(x) cell2mat(x),lifeTimeFAsGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(lifeTimeFAsGroupCell,nameList,curMovie.timeInterval_/60,false,true);
%     barPlotCellArray(lifeTimeFAsGroupCell,nameList,curMovie.timeInterval_/60);
if plotSuccess
    ylabel('Life time of FAs (min)')
    title(['Life time of FAs'])
    hgexport(h1,[figPath filesep 'lifeTimeFAs'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeFAs'],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeFAs'],'-dtiff')

    tableStableNAFCratio=table(stableNAFCratioGroup,'RowNames',nameList);
    writetable(tableStableNAFCratio,[dataPath filesep 'lifeTimeFAs.csv'],'WriteRowNames',true)
end
%% Percentage of life time
percLT=60;
perc=percLT/100;
lifeTimeFAsGroupCellSmall=cell(numConditions,1);
for ii=1:numConditions
    lifeTimeFAsGroupCellSmall{ii,1} = ...
    quantile(lifeTimeMaturingNAsGrouppCell{ii},(1-perc)+(perc-0.01)*rand(1,round((perc-0.01)*sum(~isnan(lifeTimeMaturingNAsGrouppCell{ii})))));
end
% lifeTimeFAsGroupCellSmall{1} = ...
% quantile(lifeTimeFAsGroupCell{1},(1-perc)+(perc-0.01)*rand(1,round((perc-0.01)*sum(~isnan(lifeTimeFAsGroupCell{1})))));
% lifeTimeFAsGroupCellSmall{2,1} = ...
% quantile(lifeTimeMaturingNAsGrouppCell{2},(1-perc)+(perc-0.01)*rand(1,round((perc-0.01)*sum(~isnan(lifeTimeMaturingNAsGrouppCell{2})))));
h1=figure; 
plotSuccess=boxPlotCellArray(lifeTimeFAsGroupCellSmall,nameList,curMovie.timeInterval_/60,false,true);
if plotSuccess
    ylabel(['Life time (min)'])
    title(['Life time of FAs (top ' num2str(percLT) ' percentile)'])
    ylim auto
    hgexport(h1,[figPath filesep 'lifeTimeFATop' num2str(percLT)],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'lifeTimeFATop' num2str(percLT)],'-v7.3')
    print(h1,[figPath filesep 'lifeTimeFATop' num2str(percLT)],'-dtiff')
    lifeTimeFATop20=table(lifeTimeFAsGroupCellSmall,'RowNames',nameList);
    writetable(lifeTimeFATop20,[dataPath filesep 'lifeTimeFATop' num2str(percLT) '.csv'],'WriteRowNames',true)    
end
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
plotSuccess=boxPlotCellArray(intensityGroupAll,nameListAdhComb,1,false,true);
if plotSuccess
    ylabel(['Fluorescence Intensity (a.u.)'])
    title(['F.I. of the other channel'])
    hgexport(h1,[figPath filesep 'intenTheOtherChannel'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'intenTheOtherChannel'],'-v7.3')
    print(h1,[figPath filesep 'intenTheOtherChannel'],'-dtiff')

    tableIntensityTheOther=table(intensityGroupAll,'RowNames',nameListAdhComb);
    writetable(tableIntensityTheOther,[dataPath filesep 'intenTheOtherChannel.csv'],'WriteRowNames',true)    
end
%% amplitude of the other channel at each adhesion type
h1=figure; 
amplitudeInNAsGroupCell = cellfun(@(x) cell2mat(x),amplitudeInNAsGroup,'unif',false);
amplitudeInFCsGroupCell = cellfun(@(x) cell2mat(x),amplitudeInFCsGroup,'unif',false);
amplitudeInFAsGroupCell = cellfun(@(x) cell2mat(x),amplitudeInFAsGroup,'unif',false);
amplitudeGroupAll=cell(numel(amplitudeInNAsGroup)*3,1);
amplitudeGroup={amplitudeInNAsGroupCell, amplitudeInFCsGroupCell, amplitudeInFAsGroupCell};

for ii=1:numel(intensitiesInNAsGroup)
    p=ii-1;
    for jj=1:3
        nameListAdhComb{3*p+jj,1} = [nameList{ii} '-' nameListAdh{jj}];
        amplitudeGroupAll{3*p+jj,1} = amplitudeGroup{jj}{ii};
    end

end
plotSuccess=boxPlotCellArray(amplitudeGroupAll,nameListAdhComb,1,false,true);
if plotSuccess
    ylabel('Fluorescence amplitude (a.u.)')
    title('Background-subtracted amplitude of the other channel')
    hgexport(h1,[figPath filesep 'amplitudeTheOther'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'amplitudeTheOther'],'-v7.3')
    print(h1,[figPath filesep 'amplitudeTheOther'],'-dtiff')

    tableAmplitudeTheOther=table(amplitudeGroupAll,'RowNames',nameListAdhComb);
    writetable(tableAmplitudeTheOther,[dataPath filesep 'amplitudeTheOther.csv'],'WriteRowNames',true)    
end
%% nucleatingNARatioGroup
nucleatingNARatioGroupCell = cellfun(@(x) cell2mat(x),nucleatingNARatioGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(nucleatingNARatioGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel(['nucleating NA Ratio (1)'])
    title(['newly nucleating NA Ratio'])
    hgexport(h1,[figPath filesep 'nucleatingNARatio'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'nucleatingNARatio'],'-v7.3')
    print(h1,[figPath filesep 'nucleatingNARatio'],'-dtiff')

    tableNucleatingNARatio=table(nucleatingNARatioGroup,'RowNames',nameList);
    writetable(tableNucleatingNARatio,[dataPath filesep 'nucleatingNARatio.csv'],'WriteRowNames',true)
end
%% Profile plot  - ampMature
ampArrayAllMature = cellfun(@(x) cell2mat(x),ampMaturingArrayGroup,'unif',false);
h1=figure; hold on
colorList{1} = [255/255, 204/255, 103/255];
colorList{2} = [181/255, 217/255, 148/255];
colorList{3} = [20/255, 110/255, 255/255];
colorList{4} = [255/255, 128/255, 255/255];

colorMeanList{1} = [248/255, 152/255, 56/255];
colorMeanList{2} = [13/255, 159/255, 73/255];
colorMeanList{3} = [10/255, 60/255, 200/255];
colorMeanList{4} = [180/255, 72/255, 180/255];

nSampleFrames = size(ampArrayAllMature{1},2);
t = (0:nSampleFrames-1)*curMovie.timeInterval_;


for ii=1:numel(ampArrayAllMature)
    plot(t,ampArrayAllMature{ii}','Color',colorList{ii})
end
% For mean
for ii=1:numel(ampArrayAllMature)
    curAmpAvg = nanmean(ampArrayAllMature{ii},1)';
    pM(ii)=plot(t,curAmpAvg','Color',colorMeanList{ii},'LineWidth',2);
end
legend(pM,nameList')
ylabel('Amplitude (a.u.)')
xlabel('Time (sec)')
title({'Amplitude timeseries for'; ...
    ['first ' num2str(round(nSampleFrames*curMovie.timeInterval_/60)) ' min at maturing NAs']})
hgsave(h1,[figPath filesep 'ampArrayMature'],'-v7.3')
print(h1,[figPath filesep 'ampArrayMature'],'-dtiff')
close(h1)
%% Profile plot  - ampFail
ampArrayAllFail = cellfun(@(x) cell2mat(x),ampFailingArrayGroup,'unif',false);
h1=figure; hold on

for ii=1:numel(ampArrayAllFail)
    plot(t,ampArrayAllFail{ii}','Color',colorList{ii})
end
% For mean
for ii=1:numel(ampArrayAllFail)
    curAmpAvg = nanmean(ampArrayAllFail{ii},1)';
    pM(ii)=plot(t,curAmpAvg','Color',colorMeanList{ii},'LineWidth',2);
end
legend(pM,nameList')
ylabel('Amplitude (a.u.)')
xlabel('Time (sec)')
title({'Amplitude timeseries for'; ...
    ['first ' num2str(round(nSampleFrames*curMovie.timeInterval_/60)) ' min at failing NAs']})
hgsave(h1,[figPath filesep 'ampArrayFail'],'-v7.3')
print(h1,[figPath filesep 'ampArrayFail'],'-dtiff')
close(h1)
%% Profile plot  - ampEmerg
ampArrayAllEmerge = cellfun(@(x) cell2mat(x),ampEmergingArrayGroup,'unif',false);
h1=figure; hold on

for ii=1:numel(ampArrayAllEmerge)
    plot(t,ampArrayAllEmerge{ii}','Color',colorList{ii})
end
% For mean
for ii=1:numel(ampArrayAllEmerge)
    curAmpAvg = nanmean(ampArrayAllEmerge{ii},1)';
    pM(ii)=plot(t,curAmpAvg','Color',colorMeanList{ii},'LineWidth',2);
end
legend(pM,nameList')
ylabel('Amplitude (a.u.)')
xlabel('Time (sec)')
title({'Amplitude timeseries for'; ...
    ['first ' num2str(round(nSampleFrames*curMovie.timeInterval_/60)) ' min at all emerging NAs']})
hgsave(h1,[figPath filesep 'ampArrayEmerge'],'-v7.3')
print(h1,[figPath filesep 'ampArrayEmerge'],'-dtiff')
close(h1)
%% fraction of each class
allNumGroup = cellfun(@(a,b,c,d,e,f,g,h,i) arrayfun(@(a2,b2,c2,d2,e2,f2,g2,h2,i2) ...
    a2+b2+c2+d2+e2+f2+g2+h2+i2,a,b,c,d,e,f,g,h,i),numG1Group,numG2Group,numG3Group,...
    numG4Group,numG5Group,numG6Group,numG7Group,numG8Group,numG9Group,'unif',false);
fractionGroupG1 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG1Group, allNumGroup,'unif',false);
fractionGroupG2 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG2Group, allNumGroup,'unif',false);
fractionGroupG3 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG3Group, allNumGroup,'unif',false);
fractionGroupG4 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG4Group, allNumGroup,'unif',false);
fractionGroupG5 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG5Group, allNumGroup,'unif',false);
fractionGroupG6 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG6Group, allNumGroup,'unif',false);
fractionGroupG7 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG7Group, allNumGroup,'unif',false);
fractionGroupG8 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG8Group, allNumGroup,'unif',false);
fractionGroupG9 = cellfun(@(a,b) arrayfun(@(a2,b2) ...
    a2/b2,a,b),numG9Group, allNumGroup,'unif',false);
avgFractionGroupG1=cellfun(@mean,fractionGroupG1);
avgFractionGroupG2=cellfun(@mean,fractionGroupG2);
avgFractionGroupG3=cellfun(@mean,fractionGroupG3);
avgFractionGroupG4=cellfun(@mean,fractionGroupG4);
avgFractionGroupG5=cellfun(@mean,fractionGroupG5);
avgFractionGroupG6=cellfun(@mean,fractionGroupG6);
avgFractionGroupG7=cellfun(@mean,fractionGroupG7);
avgFractionGroupG8=cellfun(@mean,fractionGroupG8);
avgFractionGroupG9=cellfun(@mean,fractionGroupG9);

stdFractionGroupG1=cellfun(@std,fractionGroupG1);
stdFractionGroupG2=cellfun(@std,fractionGroupG2);
stdFractionGroupG3=cellfun(@std,fractionGroupG3);
stdFractionGroupG4=cellfun(@std,fractionGroupG4);
stdFractionGroupG5=cellfun(@std,fractionGroupG5);
stdFractionGroupG6=cellfun(@std,fractionGroupG6);
stdFractionGroupG7=cellfun(@std,fractionGroupG7);
stdFractionGroupG8=cellfun(@std,fractionGroupG8);
stdFractionGroupG9=cellfun(@std,fractionGroupG9);

avgFractionGroupAll=[avgFractionGroupG1 avgFractionGroupG2 avgFractionGroupG3 avgFractionGroupG4 ...
    avgFractionGroupG5 avgFractionGroupG6 avgFractionGroupG7 avgFractionGroupG8 avgFractionGroupG9];
h1=figure; 
bar(avgFractionGroupAll,'stacked');legend('g1','g2','g3','g4','g5','g6','g7','g8','g9')
%% assemRateGroup for all
assemRateGroupCell = cellfun(@(x) cell2mat(x),assemRateGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(assemRateGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel('assembly rate (1/min)')
    title(['assembly rate all adhesions'])
    hgexport(h1,[figPath filesep 'assemRateAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'assemRateAll'],'-v7.3')
    print(h1,[figPath filesep 'assemRateAll'],'-dtiff')

    tableAssemRateGroup=table(assemRateGroup,'RowNames',nameList);
    writetable(tableAssemRateGroup,[dataPath filesep 'assemRateAll.csv'],'WriteRowNames',true)
end
%% disassemRateGroup for all
disassemRateGroupCell = cellfun(@(x) cell2mat(x),disassemRateGroup,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(disassemRateGroupCell,nameList,1,false,true);
if plotSuccess
    ylabel('disassembly rate (1/min)')
    title(['disassembly rate all adhesions'])
    hgexport(h1,[figPath filesep 'disassemRateAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'disassemRateAll'],'-v7.3')
    print(h1,[figPath filesep 'disassemRateAll'],'-dtiff')

    tableDisassemRateGroup=table(disassemRateGroupCell,'RowNames',nameList);
    writetable(tableDisassemRateGroup,[dataPath filesep 'disassemRateAll.csv'],'WriteRowNames',true)
end
%% assemRate2Group for all
assemRateGroup2Cell = cellfun(@(x) cell2mat(x'),assemRate2Group,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(assemRateGroup2Cell,nameList,1,false,true);
if plotSuccess
    ylabel('assembly rate (1/min)')
    title(['assembly rate all adhesions of the other channel'])
    hgexport(h1,[figPath filesep 'assemRate2All'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'assemRate2All'],'-v7.3')
    print(h1,[figPath filesep 'assemRate2All'],'-dtiff')

    tableAssemRate2Group=table(assemRate2Group,'RowNames',nameList);
    writetable(tableAssemRate2Group,[dataPath filesep 'assemRate2All.csv'],'WriteRowNames',true)
end
%% disassemRateGroup2 for all
disassemRate2GroupCell = cellfun(@(x) cell2mat(x'),disassemRate2Group,'unif',false);
h1=figure; 
plotSuccess=boxPlotCellArray(disassemRate2GroupCell,nameList,1,false,true);
if plotSuccess
    ylabel('disassembly rate (1/min)')
    title(['disassembly rate all adhesions of the other channel'])
    hgexport(h1,[figPath filesep 'disassemRate2All'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'disassemRate2All'],'-v7.3')
    print(h1,[figPath filesep 'disassemRate2All'],'-dtiff')

    tableDisassemRate2Group=table(disassemRate2GroupCell,'RowNames',nameList);
    writetable(tableDisassemRate2Group,[dataPath filesep 'disassemRate2All.csv'],'WriteRowNames',true)
end
%% Plotting each - peakGroup - all classes - usually not interesting: there is no lag.
if ~isempty(initRiseProc) && ~isempty(peakGroup{1}{1})
    nonZeroGroups=~cellfun(@isempty,assemRateEachGroup{1});
    for curGroup=find(nonZeroGroups)'
        try
          peakLagGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),peakGroup,'unif',false);
          h1=figure; 
          plotSuccess=boxPlotCellArray(peakLagGroupEach,nameList,1,false,true);
          if plotSuccess
              ylabel(['Peak lag in group ' num2str(curGroup) ' (sec)'])
              title(['Peak lag in group ' num2str(curGroup)])
              hgexport(h1,[figPath filesep 'peakLagG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
              hgsave(h1,[figPath filesep 'peakLagG' num2str(curGroup)],'-v7.3')
              print(h1,[figPath filesep 'peakLagG' num2str(curGroup)],'-dtiff')

              tablePeakLagGroupEach=table(peakLagGroupEach,'RowNames',nameList);
              writetable(tablePeakLagGroupEach,[dataPath filesep 'peakLagG' num2str(curGroup) '.csv'],'WriteRowNames',true)
          end
        catch
          disp(['No data in ' num2str(curGroup) 'th group in peakGroup'])
        end
    end
    %% Plotting each - endTimeGroup - all classes - usually not interesting: there is no lag.
    for curGroup=find(nonZeroGroups)'
        try
            endLagGroupEach = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x{curGroup}','unif',false)),endTimeGroup,'unif',false);
            h1=figure; 
            plotSuccess=boxPlotCellArray(endLagGroupEach,nameList,1,false,true);
            if plotSuccess
                ylabel(['Ending Time lag in group ' num2str(curGroup) ' (sec)'])
                title(['Ending Time lag in group ' num2str(curGroup)])
                hgexport(h1,[figPath filesep 'endLagG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
                hgsave(h1,[figPath filesep 'endLagG' num2str(curGroup)],'-v7.3')
                print(h1,[figPath filesep 'endLagG' num2str(curGroup)],'-dtiff')

                tableEndLagGroupEach=table(endLagGroupEach,'RowNames',nameList);
                writetable(tableEndLagGroupEach,[dataPath filesep 'endLagG' num2str(curGroup) '.csv'],'WriteRowNames',true)
            end
        catch
            disp(['No data in ' num2str(curGroup) 'th group in endTimeGroup'])
        end
    end
end
if ~isempty(initRiseProc) && ~isempty(assemRateEachGroup{1}{1})
    nonZeroGroups = ~cellfun(@isempty, assemRateEachGroup{1});
    %% assemRateGroup
    for curGroup=find(nonZeroGroups)'
        try
            assemRateGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),assemRateEachGroup,'unif',false);
            % Discarding not assemblying adhesions
            noAssembling = cell(numel(assemRateGroupEach),1);
            for ii=1:numel(assemRateGroupEach)
                noAssembling{ii} = assemRateGroupEach{ii}<=0;
                assemRateGroupEach{ii}(noAssembling{ii})=[];
            end
            h1=figure; 
            boxPlotCellArray(assemRateGroupEach,nameList,1,false,true);
            ylabel('assembly rate (1/min)')
            title(['assembly rate in group ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'assemRateG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'assemRateG' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'assemRateG' num2str(curGroup)],'-dtiff')

            tableAssemRate=table(assemRateGroupEach,'RowNames',nameList);
            writetable(tableAssemRate,[dataPath filesep 'assemRateG' num2str(curGroup) '.csv'],'WriteRowNames',true)
        catch
            disp(['Nothing in ' num2str(curGroup) 'th group'])
        end
    end
    %% assemRate for all - again from step 11 quantification
    assemRateAll = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x,'unif',false)),assemRateEachGroup,'unif',false); %need to work on this
    % Discarding not assemblying adhesions
    noAssembling = cell(numel(assemRateAll),1);
    h1=figure; 
    boxPlotCellArray(assemRateAll,nameList,1,false,true);
    ylabel('Assembly rate (1/min)')
    title(['Assembly rate in all groups'])
    hgexport(h1,[figPath filesep 'assemRateAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'assemRateAll'],'-v7.3')
    print(h1,[figPath filesep 'assemRateAll'],'-dtiff')

    tableAssemRate=table(assemRateAll,'RowNames',nameList);
    writetable(tableAssemRate,[dataPath filesep 'assemRateAll.csv'],'WriteRowNames',true)
    
    %% disassemRateGroup
    for curGroup=find(nonZeroGroups)'
        try
            disassemRateGroupEach = cellfun(@(x) cell2mat(x{curGroup}'),disassemRateEachGroup,'unif',false);
            % Discarding not-disassemblying adhesions
            noDisassembling = cell(numel(disassemRateGroupEach),1);
            for ii=1:numel(disassemRateGroupEach)
                noDisassembling{ii} = disassemRateGroupEach{ii}<=0;
                disassemRateGroupEach{ii}(noDisassembling{ii})=[];
            end
            h1=figure; 
            boxPlotCellArray(disassemRateGroupEach,nameList,1,false,true);
            ylabel('disassembly rate (1/min)')
            title(['disassembly rate in group ' num2str(curGroup)])
            hgexport(h1,[figPath filesep 'disassemRateG' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'disassemRateG' num2str(curGroup)],'-v7.3')
            print(h1,[figPath filesep 'disassemRateG' num2str(curGroup)],'-dtiff')

            tableDisassemRate=table(disassemRateGroupEach,'RowNames',nameList);
            writetable(tableDisassemRate,[dataPath filesep 'disassemRateG' num2str(curGroup) '.csv'],'WriteRowNames',true)
        catch
            disp(['Nothing in ' num2str(curGroup) 'th group'])
        end
    end
    %% disassemRate for all - again from step 11 quantification
    disassemRateAll = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x,'unif',false)),disassemRateEachGroup,'unif',false); %need to work on this
    h1=figure; 
    boxPlotCellArray(disassemRateAll,nameList,1,false,true);
    ylabel('disassembly rate (1/min)')
    title('disassembly rate in all groups ')
    hgexport(h1,[figPath filesep 'disassemRateAll'],hgexport('factorystyle'),'Format','eps')
    hgsave(h1,[figPath filesep 'disassemRateAll'],'-v7.3')
    print(h1,[figPath filesep 'disassemRateAll'],'-dtiff')

    tableDisassemRate=table(disassemRateAll,'RowNames',nameList);
    writetable(tableDisassemRate,[dataPath filesep 'disassemRateAll.csv'],'WriteRowNames',true)
    %% assemRate2Group
    if exist('assemRate2EachGroup','var') && any(~cellfun(@isempty,assemRate2EachGroup))
        try
            nonZeroGroups2 = ~cellfun(@isempty, assemRate2EachGroup{1});
        catch
            nonZeroGroups2 = ~cellfun(@isempty, assemRate2EachGroup); %This will anyway generates 0
        end

        for curGroup=find(nonZeroGroups2)
            %try
                assemRate2GroupEach = cellfun(@(x) cell2mat(x{curGroup}'),assemRate2EachGroup,'unif',false);
    %             % Discarding not assemblying adhesions - this doesn't make
    %             sense because the other channel might not necessarily
    %             correlate with the main adhesion channel.
                for ii=1:numel(assemRate2GroupEach)
                    assemRate2GroupEach{ii}(noAssembling{ii})=[];
                end
                h1=figure; 
                boxPlotCellArray(assemRate2GroupEach,nameList,1,false,true);
                ylabel('assembly rate (1/min)')
                title(['assembly rate in group ' num2str(curGroup) ' of the other channel'])
                hgexport(h1,[figPath filesep 'assemRate2_G' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
                hgsave(h1,[figPath filesep 'assemRate2_G' num2str(curGroup)],'-v7.3')
                print(h1,[figPath filesep 'assemRate2_G' num2str(curGroup)],'-dtiff')

                tableAssemRate2=table(assemRate2GroupEach,'RowNames',nameList);
                writetable(tableAssemRate2,[dataPath filesep 'assemRate2_G' num2str(curGroup) '.csv'],'WriteRowNames',true)
            %catch
                disp(['Nothing in ' num2str(curGroup) 'th group'])
            %end
        end
    
    %% assemRate for all - again from step 11 quantification
        assemRate2All = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x','unif',false)),assemRate2EachGroup,'unif',false); %need to work on this
        % Discarding not assemblying adhesions
        noAssembling = cell(numel(assemRate2All),1);
        h1=figure; 
        boxPlotCellArray(assemRate2All,nameList,1,false,true);
        ylabel('Assembly rate (1/min)')
        title(['Assembly rate in all groups of the other channel'])
        hgexport(h1,[figPath filesep 'assemRate2All'],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'assemRate2All'],'-v7.3')
        print(h1,[figPath filesep 'assemRate2All'],'-dtiff')

        tableAssemRate2=table(assemRate2All,'RowNames',nameList);
        writetable(tableAssemRate2,[dataPath filesep 'assemRate2All.csv'],'WriteRowNames',true)

        %% disassemRate2Group
        for curGroup=find(nonZeroGroups2)
            %try
                disassemRate2GroupEach = cellfun(@(x) cell2mat(x{curGroup}'),disassemRate2EachGroup,'unif',false);
    %             % Discarding not-disassemblying adhesions- this doesn't make
    %             sense because the other channel might not necessarily
    %             correlate with the main adhesion channel.
                for ii=1:numel(disassemRate2GroupEach)
                    disassemRate2GroupEach{ii}(noDisassembling{ii})=[];
                end
                h1=figure; 
                boxPlotCellArray(disassemRate2GroupEach,nameList,1,false,true);
                ylabel('disassembly rate (1/min)')
                title(['disassembly rate in group ' num2str(curGroup) ' of the other channel'])
                hgexport(h1,[figPath filesep 'disassemRate2_G' num2str(curGroup)],hgexport('factorystyle'),'Format','eps')
                hgsave(h1,[figPath filesep 'disassemRate2_G' num2str(curGroup)],'-v7.3')
                print(h1,[figPath filesep 'disassemRate2_G' num2str(curGroup)],'-dtiff')

                tableDisassem2Rate=table(disassemRateGroupEach,'RowNames',nameList);
                writetable(tableDisassem2Rate,[dataPath filesep 'disassemRate2_G' num2str(curGroup) '.csv'],'WriteRowNames',true)
            %catch
                disp(['Nothing in ' num2str(curGroup) 'th group'])
            %end
        end
        %% disassemRate for all - again from step 11 quantification
        disassemRate2All = cellfun(@(x) cell2mat(cellfun(@(y) cell2mat(y'),x','unif',false)),disassemRate2EachGroup,'unif',false); %need to work on this
        % Discarding not disassemblying adhesions
        nodisassembling = cell(numel(disassemRate2All),1);
        h1=figure; 
        boxPlotCellArray(disassemRate2All,nameList,1,false,true);
        ylabel('Disassembly rate (1/min)')
        title(['Disassembly rate in all groups of the other channel'])
        hgexport(h1,[figPath filesep 'disassemRate2All'],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'disassemRate2All'],'-v7.3')
        print(h1,[figPath filesep 'disassemRate2All'],'-dtiff')

        tableDisassemRate2=table(disassemRate2All,'RowNames',nameList);
        writetable(tableDisassemRate2,[dataPath filesep 'disassemRate2All.csv'],'WriteRowNames',true)
    end
end
%% save entire workspace for later
close all
save([dataPath filesep 'allData.mat'])
disp('Done plotting data!')