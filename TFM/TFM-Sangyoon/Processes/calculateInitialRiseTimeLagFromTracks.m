function [] = calculateInitialRiseTimeLagFromTracks(MD,varargin)
% calculateInitialRiseTimeLagFromTracks calculates time lag between the
% main indexed channel vs. channels read in TheOtherChannelReadingProcess
% or TractionForceReadingProcess.
% Sangyoon Han April 2013

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('MD', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(MD,varargin{:});
paramsIn=ip.Results.paramsIn;
%Parse input, store in parameter structure
%Get the indices of any previous threshold processes from this function   
try
    iProc = MD.getProcessIndex('InitialRiseTimeLagCalculationProcess',1,0);
    %If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(MD.processes_)+1;
        MD.addProcess(InitialRiseTimeLagCalculationProcess(MD));                                                                                                 
    end
    timeLagProc = MD.processes_{iProc};
catch
    iFAPack =  MD.getPackageIndex('FocalAdhesionPackage');
    FAPack = MD.getPackage(iFAPack);
    timeLagProc = FAPack.processes_{11};
end
p = parseProcessParams(timeLagProc,paramsIn);
iBeadChan = 1; % might need to be updated based on asking TFMPackage..
p.measureAmp2rates = true; % Should be replaced with a user input.

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(MD,showAllTracks,plotEachTrack,varargin{:});
% iChanMaster=p.ChannelIndex;

%% Reading tracks from master channel
% Use the previous analysis folder structure
% It might be good to save which process the tracksNA was obtained.
disp('Reading tracksNA ...')
iFAPack =  MD.getPackageIndex('FocalAdhesionPackage');
FAPack = MD.getPackage(iFAPack);

tic
try
    iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
    adhAnalProc = MD.getProcess(iAdhProc);
catch
    adhAnalProc = FAPack.processes_{7};
end
% Doublecheck the p.ChannelIndex
p.ChannelIndex = adhAnalProc.funParams_.ChannelIndex;
tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA');
% numChans = numel(p.ChannelIndex);
%% Doublecheck if amp has similar range as amp2 and read it again at least for the extra part
% find one examplary track that has intermediate starting point
startingFrameExtraAll = arrayfun(@(x) x.startingFrameExtra,tracksNA);
startingFrameExtraExtraAll = arrayfun(@(x) x.startingFrameExtraExtra,tracksNA);
intermedTrackIDs = find(startingFrameExtraExtraAll<(startingFrameExtraAll-10));
if ~isempty(intermedTrackIDs)
    trackInspected = tracksNA(intermedTrackIDs(1));
    % See if amp is NaN before startingFrameExtra (after
    % startingFrameExtraExtra)
    if isnan(trackInspected.amp(trackInspected.startingFrameExtraExtra+1))
        % then we need to read this amp again
        disp('Reading image stack again to read extra time regime of amp...')
        tic
        imgStack = getAnyStacks(MD);
        toc
        disp('Reading from tracks'); tic
        tracksNA = readIntensityFromTracks(tracksNA,imgStack,1,'extraReadingOnly',true); toc;
        disp('Saving the tracks...'); tic
        numTracks=numel(tracksNA);
        fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
        numStr = @(trackNum) num2str(trackNum,fString);
        trackFolderPath = [adhAnalProc.funParams_.OutputDirectory filesep 'trackIndividual'];
        trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
        
        parfor k=1:numTracks
            curTrack=tracksNA(k);
            parsave(feval(trackIndPath,k),curTrack)
        end
        toc
    end
else
    disp('All tracks have the same or similar starting points')
end

%% Now we have to combine this with readings from step 9 and 10
iTheOtherProc=9; 
theOtherReadProc=FAPack.processes_{iTheOtherProc};

if ~isempty(theOtherReadProc)
    ampObj = load(theOtherReadProc.outFilePaths_{1,p.ChannelIndex},'tracksAmpTotal'); % the later channel has the most information.
    tracksAmpTotal = ampObj.tracksAmpTotal;
%     tracksAmpTotal = theOtherReadProc.loadChannelOutput(p.ChannelIndex);
    
    if isfield(tracksAmpTotal,'ampTotal2')
        [tracksNA(:).ampTotal2] = tracksAmpTotal.ampTotal2;
    end
    if isfield(tracksAmpTotal,'amp2')
        [tracksNA(:).amp2] = tracksAmpTotal.amp2;
    end
    if isfield(tracksAmpTotal,'ampTotal3')
        [tracksNA(:).ampTotal3] = tracksAmpTotal.ampTotal3;
    end
end

%% Combinging from Step 10
iForceRead=10;
forceReadProc=FAPack.processes_{iForceRead};
if ~isempty(forceReadProc)
    forceReadObj = load(forceReadProc.outFilePaths_{1,p.ChannelIndex},'tracksForceMag'); % the later channel has the most information.
    tracksForceMag = forceReadObj.tracksForceMag;
    idxTracksObj = load(forceReadProc.outFilePaths_{2,p.ChannelIndex},'idxTracks');
    if ~isfield(idxTracksObj,'idxTracks')
        idxTracksObj = load(forceReadProc.outFilePaths_{6,p.ChannelIndex},'idxTracks');
    end
    idxTracks = idxTracksObj.idxTracks;
    tracksNA = tracksNA(idxTracks);
    if isfield(tracksForceMag,'forceMag')
        [tracksNA(:).forceMag] = tracksForceMag.forceMag;
    end
end
toc
%% Reading classes
disp('Reading idsClassified ...')
try
    try
        iClaProc = MD.getProcessIndex('AdhesionClassificationProcess');
        classProc = MD.getProcess(iClaProc);
    catch
        classProc = FAPack.processes_{8};
    end
    % numChans = numel(p.ChannelIndex);
    idsClassified = load(classProc.outFilePaths_{4,p.ChannelIndex});
    idGroup1 = idsClassified.idGroup1;
    idGroup2 = idsClassified.idGroup2;
    idGroup3 = idsClassified.idGroup3;
    idGroup4 = idsClassified.idGroup4;
    idGroup5 = idsClassified.idGroup5;
    idGroup6 = idsClassified.idGroup6;
    idGroup7 = idsClassified.idGroup7;
    idGroup8 = idsClassified.idGroup8;
    idGroup9 = idsClassified.idGroup9;
catch
    disp('No Classified groups. Using no classification...')
    % Potentially we can use the dedactically chosen classes here (e.g.
    % shown in the code used for Tristan's movie analysis.
    idGroup1 = []; idGroup2 = []; idGroup3 = []; idGroup4 = [];
    idGroup5 = []; idGroup6 = []; idGroup7 = []; idGroup8 = []; idGroup9 = [];
end    
if ~isempty(forceReadProc)
    idGroup1 = idGroup1(idxTracks);
    idGroup2 = idGroup2(idxTracks);
    idGroup3 = idGroup3(idxTracks);
    idGroup4 = idGroup4(idxTracks);
    idGroup5 = idGroup5(idxTracks);
    idGroup6 = idGroup6(idxTracks);
    idGroup7 = idGroup7(idxTracks);
    idGroup8 = idGroup8(idxTracks);
    idGroup9 = idGroup9(idxTracks);
else
    disp('Traction reading was not done. No further filtering...')
end
% idGroups = {idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
%% Data Set up
% Set up the output file path for master channel
outputFile = cell(6, numel(MD.channels_));
for i = p.ChannelIndex
    [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_Chan' num2str(i) '_tracksNA'];
    outputFile{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
end
mkClrDir(p.OutputDirectory);
dataPath=[p.OutputDirectory filesep 'data'];mkClrDir(dataPath);
epsPath=[p.OutputDirectory filesep 'eps'];mkClrDir(epsPath);
figPath=[p.OutputDirectory filesep 'figs'];mkClrDir(figPath);
tifPath=[p.OutputDirectory filesep 'tif'];mkClrDir(tifPath);

%% Initial force-increase time quantification! - time
tInterval = MD.timeInterval_;
potentialSlaves = {'forceMag','ampTotal2','ampTotal3'};
existingSlaveIDs = isfield(tracksNA,potentialSlaves);

for ii=1:numel(tracksNA)
    tracksNA(ii).ampTotal(tracksNA(ii).ampTotal==0)=NaN;
end
%% Another output setup
existingSlaveIDs=find(existingSlaveIDs);
for i = existingSlaveIDs
    [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_initialTimeDelayIntAgainst_' potentialSlaves{i}];
    outputFile{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
    outFilename = [chanDirName 'peakTimeIntAgainst_' potentialSlaves{i}];
    outputFile{3,i} = [p.OutputDirectory filesep outFilename '.mat'];
    outFilename = [chanDirName '_endingIntAgainst_' potentialSlaves{i}];
    outputFile{4,i} = [p.OutputDirectory filesep outFilename '.mat'];
    outFilename = [chanDirName '_crossCorrelationScore_' potentialSlaves{i}];
    outputFile{5,i} = [p.OutputDirectory filesep outFilename '.mat'];
    outFilename = [chanDirName '_crossCorrelationTimeLag_' potentialSlaves{i}];
    outputFile{6,i} = [p.OutputDirectory filesep outFilename '.mat'];
end
timeLagProc.setOutFilePaths(outputFile);

%% I. Time series analyses
p.numWinSize=0; %frames
p.preDetecPeriod=60; %sec
for jj=existingSlaveIDs
    curSlaveCell = potentialSlaves(jj);
    curSlave=curSlaveCell{1};
    % initial time estimation
    %     preDetecFactor=3/5; %for paxillin because diffuse signal before point-source-like signal
    %     preDetecFactor=1/5; %for vinculin 
    preDetecFactor=1/10; %for talin
    % 1. Initial force-increase time quantification! - time
    [curFirstIncreseTimeIntAgainstSlave,SlaveTransmittingAll{jj}...
    ,firstIncreseTimeIntAll{jj},firstIncreseTimeSlaveAll{jj},bkgMaxIntAll{jj},bkgMaxSlaveAll{jj},~,ampIncreasing] ...
        = calculateFirstIncreaseTimeTracks(tracksNA,p.numWinSize,p.preDetecPeriod,tInterval,'slaveSource',curSlave);
    firstIncreseTimeIntAgainstSlaveAll{jj}=curFirstIncreseTimeIntAgainstSlave;
    disp(['Median of firstIncreseTimeIntAgainst' curSlave 'All = ' num2str(nanmedian(firstIncreseTimeIntAgainstSlaveAll{jj}))])
    % 2. Peak intensity lag against force
    %   perform spline filter for ampTotal
    splineParam=0.1;
    [curPeakTimeIntAgainstSlave,peakTimeIntAll{jj},peakTimeSlaveAll{jj}]...
        = calculatePeakTimeLagFromTracks(tracksNA,splineParam,tInterval,'slaveSource',curSlave);
    peakTimeIntAgainstSlaveAll{jj}=curPeakTimeIntAgainstSlave;
    % statistics about peakTimeIntAgainstSlaveAll
%     figure, histogram(peakTimeIntAgainstSlaveAll{jj})
    disp(['Median of peakTimeIntAgainst' curSlave 'All = ' num2str(nanmedian(peakTimeIntAgainstSlaveAll{jj}))])
    % 3. Ending time lag: how much early the main signal ends against when the slave signal ends
    [curEndingIntAgainstSlave,endingTimeIntAll{jj},endingTimeSlaveAll{jj}]...
        = calculateEndingTimeLagFromTracks(tracksNA,splineParam,preDetecFactor,tInterval,'slaveSource',curSlave);
    endingIntAgainstSlaveAll{jj}=curEndingIntAgainstSlave;
    disp(['Median of endingIntAgainstSlaveAll' curSlave 'All = ' num2str(nanmedian(endingIntAgainstSlaveAll{jj}))])
    % 4. Cross-correlation: score and 5. lag
    for k=1:numel(tracksNA)
        presIdx = tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra;%logical(tracksNA(k).presence);
        maxLag = ceil(tracksNA(k).lifeTime/2);
        [curCC,~,curLag]  = nanCrossCorrelation(tracksNA(k).ampTotal(presIdx),...
            getfield(tracksNA(k),{1},curSlave,{presIdx}),'corrType','Pearson','maxLag',maxLag);
        [CCscoreMax(k), curMaxInd] = max(curCC);
        CCmaxLag(k)=curLag(curMaxInd)*tInterval; % in sec
    end
    CCscoreMaxAll{jj} = CCscoreMax;
    CClagAll{jj} = CCmaxLag;
    disp(['Median of CCmaxLag' curSlave 'All = ' num2str(median(CCmaxLag))])
    % Save these shortly
    save(outputFile{2,jj},'curFirstIncreseTimeIntAgainstSlave','-v7.3'); 
    save(outputFile{3,jj},'curPeakTimeIntAgainstSlave','-v7.3'); 
    save(outputFile{4,jj},'curEndingIntAgainstSlave','-v7.3'); 
    save(outputFile{5,jj},'CCscoreMax','-v7.3'); 
    save(outputFile{6,jj},'CCmaxLag','-v7.3'); 
end
%% II. Separate the processed data to all different classes
% Filtering for group1
% lateAmpSlopeG1 = arrayfun(@(x) x.lateAmpSlope,tracksNA(idGroup1));
% idxLateAmpSlopeG1 = lateAmpSlopeG1<0;
% We filter out tracks whose ampTotal is too high, bigger than mean value
% of ampTotal maxima
% meanAmpMaximum = mean(arrayfun(@(x) nanmax(x.ampTotal),tracksNA(idGroup1)));
% ampSlopeG1 = arrayfun(@(x) x.earlyAmpSlope,tracksNA(idGroup1));
% lateAmpTotalG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNA(idGroup1));
% initForceG1 = arrayfun(@(x) x.forceMag(x.startingFrame),tracksNA(idGroup1));
% idxIncreasingAmpG1 = ampSlopeG1>0;
% idxLowInitForceG1= initForceG1<500;
% idxLateAmpLow = lateAmpTotalG1<meanAmpMaximum;
% idGroup1f = idxLateAmpLow & idxIncreasingAmpG1 & idxLowInitForceG1;

iForceReadProc = 10;
forceReadProc = FAPack.processes_{iForceReadProc};
if ~isempty(forceReadProc)
    % Filtering for group1
    idGroup1FT = getForceTransmittingG1(idGroup1,tracksNA(idGroup1));
    % Filtering for group2
    idGroup2long = getForceTransmittingG2(idGroup2,tracksNA(idGroup2),MD.timeInterval_);

    % Potential G2 in G1
    idGroup2fromG1 = getForceTransmittingG2(idGroup1,tracksNA(idGroup1),MD.timeInterval_);
    % Potential G1 in G2
    idGroup1fromG2 = idGroup2 & ~idGroup2long;

    idGroup2f = idGroup2long | idGroup2fromG1;
    idGroup1f = (idGroup1FT | idGroup1fromG2) & ~idGroup2f;


    idGroups = {idGroup1f,idGroup2f,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
else
    idGroups = {idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
end
save([dataPath filesep 'idGroups.mat'],'idGroups');

%% Get the fraction of how many adhesions are
% increasing the force during their lifetime.
if ismember(1,existingSlaveIDs)
    isForceTransmitting = SlaveTransmittingAll{1};
    isAmpIncreasing = ampIncreasing;
%     idGroupsBeforeFiltering = {idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
    numEachGroup = cellfun(@(x) sum(x & isAmpIncreasing),idGroups); %idGroupsBeforeFiltering);
    numEachGroup(3:9) = cellfun(@(x) sum(x),idGroups(3:9)); %G3-9 do not necessarily have to increase amplitudes during their lifetime
    numEachGroupForceTransmitting = cellfun(@(x) sum(x & isForceTransmitting),idGroups);
    fractionForceTransmitting = numEachGroupForceTransmitting./numEachGroup;
    groupLabel = categorical(arrayfun(@(x,y) ['G' num2str(x) '(N=' num2str(y) ')'],1:9,numEachGroup,'unif',false));
    
    h3=figure; 
    bar(groupLabel,fractionForceTransmitting); title('Fraction of force-transmitting NAs per group')
    ylabel('Fraction (1)'); nameTitle = 'fractionForceTransmitting';
    hgexport(h3,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h3,strcat(figPath,filesep,nameTitle))

    save([dataPath filesep 'forceTransmitting.mat'],'fractionForceTransmitting','numEachGroup','numEachGroupForceTransmitting','groupLabel');
else
    disp('No force channel, skipping quantification of force-transmitting NA fraction...')
end
% idGroups = {idGroup1f,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
%% For loop
numClasses = numel(idGroups);
numSlaves = numel(existingSlaveIDs);
initialLagGroups=cell(numClasses,numSlaves);
peakLagGroups=cell(numClasses,numSlaves);
endingLagGroups=cell(numClasses,numSlaves);
ccScoreGroups=cell(numClasses,numSlaves);
ccLagGroups=cell(numClasses,numSlaves);
firstTimeAboveZeroBccAll=cell(2,numSlaves);
firstTimeAboveHalfBccAll=cell(2,numSlaves);
firstTimeAboveOneBccAll=cell(2,numSlaves);
BccGroups=cell(2,numSlaves);
avgBccGroups=cell(2,numSlaves);
indexValidBccInTracks=cell(2,numSlaves);
tStartingFrame=zeros(2,numSlaves);
mainTimeToPeakGroup=cell(2,numSlaves);
mainBccPeakValuesGroup=cell(2,numSlaves);
sideTimeToPeakGroup=cell(2,numSlaves);
sideBccPeakValuesGroup=cell(2,numSlaves);
for jj=existingSlaveIDs
    for k=1:numClasses
        initialLagGroups{k,jj}=firstIncreseTimeIntAgainstSlaveAll{jj}(idGroups{k});
        % For inspection: at first, I need to look at each time series (at
        % least for G1 and G2) 
        if ismember(k,[1 2])
            disp(['Current Class: ' num2str(k) ', Current Slave: ' potentialSlaves{jj}])
            [BccGroups{k,jj},avgBccGroups{k,jj},tStartingFrame(k,jj),indexValidBccInTracks{k,jj},firstTimeAboveZeroBccAll{k,jj},...
                firstTimeAboveHalfBccAll{k,jj},firstTimeAboveOneBccAll{k,jj}] = ...
                inspectNATimeSeries(tracksNA(idGroups{k}),tInterval,'slaveSource',potentialSlaves{jj}); %,p.OutputDirectory);
        end
        peakLagGroups{k,jj}=peakTimeIntAgainstSlaveAll{jj}(idGroups{k});
        endingLagGroups{k,jj}=endingIntAgainstSlaveAll{jj}(idGroups{k});
        ccScoreGroups{k,jj}=CCscoreMaxAll{jj}(idGroups{k});
        ccLagGroups{k,jj}=CClagAll{jj}(idGroups{k});
    end
end
%% III. Save time-related information per each class
save([dataPath filesep 'initialLagGroups.mat'],'initialLagGroups');
save([dataPath filesep 'peakLagGroups.mat'],'peakLagGroups');
save([dataPath filesep 'endingLagGroups.mat'],'endingLagGroups');
save([dataPath filesep 'ccScoreGroups.mat'],'ccScoreGroups');
save([dataPath filesep 'ccLagGroups.mat'],'ccLagGroups');
save([dataPath filesep 'BccGroups.mat'],'BccGroups','indexValidBccInTracks');
%% IV. Plot these wisely: if there are more than one slaves, realign events based on p.mainSlave
% If there are only ampTotal and forceMag, then it would be one box
% plotted. If there are also ampTotal2, then everything will be aligned to
% the p.mainSlave
% In case of only forceMag
for k=1:numClasses
    initialLagTogether = initialLagGroups(k,existingSlaveIDs);
    peakLagTogether = peakLagGroups(k,existingSlaveIDs);
    endingLagTogether = endingLagGroups(k,existingSlaveIDs);
    ccLagTogether = ccLagGroups(k,existingSlaveIDs);
    if ismember(k,[1 2])
        zeroBccTogether = firstTimeAboveZeroBccAll(k,existingSlaveIDs);
        halfBccTogether = firstTimeAboveHalfBccAll(k,existingSlaveIDs);
        oneBccTogether = firstTimeAboveOneBccAll(k,existingSlaveIDs);
    end
    
    nameList = potentialSlaves(existingSlaveIDs); nameList2=nameList;
    % Readjustment
    kk=0;
    initialLagTogetherAdjusted=cell(1,numel(initialLagTogether));
    peakLagTogetherAdjusted=cell(1,numel(peakLagTogether));
    endingLagTogetherAdjusted=cell(1,numel(endingLagTogether));
    ccLagTogetherAdjusted=cell(1,numel(ccLagTogether));
    for jj=existingSlaveIDs
        kk=kk+1;
        if jj==p.mainSlave
            initialLagTogetherAdjusted{kk} = initialLagTogether{existingSlaveIDs==p.mainSlave};
            peakLagTogetherAdjusted{kk} = peakLagTogether{existingSlaveIDs==p.mainSlave};
            endingLagTogetherAdjusted{kk} = endingLagTogether{existingSlaveIDs==p.mainSlave};
            ccLagTogetherAdjusted{kk} = ccLagTogether{existingSlaveIDs==p.mainSlave};
            nameList2{kk} = ['ampTotal - ' nameList{existingSlaveIDs==p.mainSlave}];
            if ismember(k,[1 2])
                zeroBccTogetherAdjusted{kk} = zeroBccTogether{existingSlaveIDs==p.mainSlave};
                halfBccTogetherAdjusted{kk} = halfBccTogether{existingSlaveIDs==p.mainSlave};
                oneBccTogetherAdjusted{kk} = oneBccTogether{existingSlaveIDs==p.mainSlave};
            end
        else
            initialLagTogetherAdjusted{kk} = initialLagTogether{existingSlaveIDs==p.mainSlave}-initialLagTogether{existingSlaveIDs==jj};
            peakLagTogetherAdjusted{kk} = peakLagTogether{existingSlaveIDs==p.mainSlave}-peakLagTogether{existingSlaveIDs==jj};
            endingLagTogetherAdjusted{kk} = endingLagTogether{existingSlaveIDs==p.mainSlave}-endingLagTogether{existingSlaveIDs==jj};
            ccLagTogetherAdjusted{kk} = ccLagTogether{existingSlaveIDs==p.mainSlave}-ccLagTogether{existingSlaveIDs==jj};
            nameList2{kk} = [potentialSlaves{existingSlaveIDs==jj} ' - ' nameList{existingSlaveIDs==p.mainSlave}];
            if ismember(k,[1 2])
                zeroBccTogetherAdjusted{kk} = zeroBccTogether{existingSlaveIDs==p.mainSlave}-zeroBccTogether{existingSlaveIDs==jj};
                halfBccTogetherAdjusted{kk} = halfBccTogether{existingSlaveIDs==p.mainSlave}-halfBccTogether{existingSlaveIDs==jj};
                oneBccTogetherAdjusted{kk} = oneBccTogether{existingSlaveIDs==p.mainSlave}-oneBccTogether{existingSlaveIDs==jj};
            end
        end
    end
    
    if all(cellfun(@isempty,initialLagTogetherAdjusted)) || all(cellfun(@(x) all(isnan(x)),initialLagTogetherAdjusted))
        disp(['The class ' num2str(k) ' has no data to analyze. Skipping...'])
        continue
    end
    
    h2=figure; ax=axes(h2);
    boxPlotCellArray(initialLagTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['initialLag Class' num2str(k)];
    title(ax,nameTitle); ylabel(ax,'Time lag (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'initialLagTogetherAdjusted','nameList2');    
    close(h2)

    h2=figure; ax=axes(h2);
    boxPlotCellArray(zeroBccTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['zeroBccTogetherAdjusted' num2str(k)];
    title(ax,nameTitle); ylabel(ax,'Time (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'zeroBccTogetherAdjusted','nameList2');    
    close(h2)

    h2=figure; ax=axes(h2);
    boxPlotCellArray(halfBccTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['halfBccTogetherAdjusted' num2str(k)];
    title(ax,nameTitle); ylabel(ax,'Time (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'halfBccTogetherAdjusted','nameList2');    
    close(h2)

    h2=figure; ax=axes(h2);
    boxPlotCellArray(oneBccTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['oneBccTogetherAdjusted' num2str(k)];
    title(ax,nameTitle); ylabel(ax,'Time (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'oneBccTogetherAdjusted','nameList2');    
    close(h2)

%     tShift = 36; % This is same as tStartingFrame = 30 + tFluc/2 + 1;
%     if ismember(k,[1 2])
%         if length(existingSlaveIDs)==1
%         	jj=existingSlaveIDs;
%             slaveSource = potentialSlaves{jj};
%             curBcc = BccGroups{k,jj};
%             if ~isempty(curBcc)
%                 idGroupIndex=find(idGroups{k});
%                 curTrackIndex = idGroupIndex(indexValidBccInTracks{k,jj});
%     %             curAvgBcc = avgBccGroups{k,jj};
%                 tRange = (1:size(curBcc,2))-tShift;
%                 h2=figure;
%                 plot(tRange,curBcc); hold on
%                 plot(tRange,nanmean(curBcc,1),'k','Linewidth',3)
%                 line([36 36],[min(min(curBcc(:,:))),max(max(curBcc(:,:)))],'LineStyle','--','Color',[0.5 .5 .5],'linewidth',4)
%                 nameTitle=['Bcc ' num2str(k) 'th class'];
%                 title(nameTitle); ylabel('Time lag (s)')
%                 hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%                 hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%                 close(h2)
%                 % Filtering done previously
%                 % Try clustering
%                 % #1 clustering with the features from my feature collection
%                 meas=extractGeneralFeatureTimeSeries(curBcc); %time of peak is the main feature
%                 numClusters = min(10,max(1,round(length(meas)/3))); 
%                 clusterArray=1:numClusters;
%                 T=kmeans(meas,numClusters);
% 
%                 % Sort T based on timeToMax
%                 meanTimeToMaxPerClus = arrayfun(@(x) mean(meas(T==x)),clusterArray);
%                 numPerClusters = arrayfun(@(x) sum(T==x),clusterArray);
%                 % Discard the cluster that has too early time to max
%                 [~,indClusterEarlyTimetoMax] = min(meanTimeToMaxPerClus);
%                 % Discard clusters whose sample sizes are small
%                 thresClusterSize = mean(numPerClusters)-0.5*std(numPerClusters);
%                 indSmallClusters = find(numPerClusters < thresClusterSize);
%                 % Choose the main cluster
%                 mainClusters = find(numPerClusters > max(numPerClusters)*0.95);
%                 sideClusters = setdiff(clusterArray,[mainClusters indSmallClusters indClusterEarlyTimetoMax]);
% 
%                 % Found that They have very similar life time per cluster
%                 % Checking lifetime... They are not necessarily the same!
%                 % Inspect
%                 mainBccPeakValues=[];
%                 mainTimeToPeak=[];
%                 sideBccPeakValues=[];
%                 sideTimeToPeak=[];
%                 for curT=[mainClusters sideClusters]
%                     h2=figure; 
%                     plot(1:size(curBcc,2),curBcc(T==curT,:),'o-')
%                     hold on
%                     plot(1:size(curBcc,2),nanmean(curBcc(T==curT,:),1),'Linewidth',3,'color','k')
%                     line([36 36],[min(min(curBcc(T==curT,:))),max(max(curBcc(T ==curT,:)))],'LineStyle','--','Color',[0.5 .5 .5],'linewidth',4)
%                     hold off
%                     if ismember(curT,mainClusters)
%                         nameTitle=['Amp-' slaveSource ' Bcc ' num2str(k) 'th class, ' num2str(curT) 'th cluster, main cluster'];
%                         mainBccPeakValues=[mainBccPeakValues; nanmax(curBcc(T==curT,:),[],2)];
%                         mainTimeToPeak = [mainTimeToPeak; (meas(T==curT)-tStartingFrame(k,jj))*tInterval];
%                     else
%                         nameTitle=['Amp-' slaveSource ' Bcc ' num2str(k) 'th class, ' num2str(curT) 'th cluster, side cluster'];
%                         sideBccPeakValues=[sideBccPeakValues; nanmax(curBcc(T==curT,:),[],2)];
%                         sideTimeToPeak = [sideTimeToPeak; (meas(T==curT)-tStartingFrame(k,jj))*tInterval];
%                     end
%                     title(nameTitle); 
%                     xlabel('Time (frame, 36th frame = time zero)')
%                     hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%                     hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%                     print(h2,'-dtiff',strcat(tifPath,filesep,nameTitle,'.tif'))
% 
%                     curTind = find(T==curT);
% 
%                     iii=0;
%                     for ii=1:sum(T==curT)
%                         iii=iii+1;
%                         if iii==1
%                             h=figure;  
%                         end
%                         curCurBcc = curBcc(curTind(ii),:);
%                         % corresponding index in tracksNA
%                         curTrackInd = curTrackIndex((curTind(ii)));
%                         curTrack = tracksNA(curTrackInd);
%                         tShift = tStartingFrame(k,jj) + curTrack.startingFrameExtraExtra-curTrack.startingFrameExtra;
%                         tFluc = 10; % this is actually in frame
%                         lastFrameCC = curTrack.endingFrameExtraExtra-tFluc; 
%                         sd=curTrack.ampTotal;
%                         subplot(2,2,1), hold off, plot(tShift-curTrack.startingFrameExtraExtra+(curTrack.startingFrameExtraExtra:lastFrameCC),sd(curTrack.startingFrameExtraExtra:lastFrameCC),'o-'), %hold on, plot(tRange(firstIncreaseTimeInt)-curTrack.startingFrameExtraExtra+1,sd(firstIncreaseTimeInt),'ro'); 
%                         title(['Fluorescent intensity (ii: ' num2str(ii) ')'])
%                         sCurForce_sd=curTrack.(slaveSource);
%                         subplot(2,2,3), hold off, plot(tShift-curTrack.startingFrameExtraExtra+(curTrack.startingFrameExtraExtra:lastFrameCC),sCurForce_sd(curTrack.startingFrameExtraExtra:lastFrameCC),'o-'), %hold on, plot(tRange(firstIncreaseTimeForce)-curTrack.startingFrameExtraExtra+1,sCurForce_sd(firstIncreaseTimeForce),'ro')    
%                         title([slaveSource ' (ii: ' num2str(ii) ')'])
%                         subplot(2,2,2), hold off, plot(tShift-curTrack.startingFrameExtraExtra+(curTrack.startingFrameExtraExtra:lastFrameCC),curTrack.distToEdge(curTrack.startingFrameExtraExtra:lastFrameCC),'o-'); 
%                         title({['Distance to edge' ' (ii: ' num2str(ii) ')']; ['Class' num2str(k)]})
%                         subplot(2,2,4), plot(1:length(curCurBcc),curCurBcc,'o-'), hold on, 
%                         curName = [slaveSource '-G-' num2str(k) '-Clst-' num2str(curT) '-numClus-' num2str(ii) '-tNum',num2str(curTrackInd)];
%                         title({curName; ['co-variance' ' (ii: ' num2str(ii) ')']})
%                         hold off
%                         hgsave(h,strcat(figPath,filesep,curName),'-v7.3')
%                         print(h,'-dtiff',strcat(tifPath,filesep,curName,'.tif'))
% %                         print(h,'-depsc','-loose',[epsPath filesep curName '.eps']);
%                         save(strcat(dataPath,filesep,curName,'.mat'),'sd','sCurForce_sd','curCurBcc')
%     %                     uiwait(); 
%                     end
%                     close(h)
%                 end
%                 save([dataPath filesep 'mainBccPeakValues-G' num2str(k) curSlave '.mat'],'mainBccPeakValues')
%                 save([dataPath filesep 'mainTimeToPeak-G' num2str(k) curSlave '.mat'],'mainTimeToPeak')
%                 save([dataPath filesep 'sideBccPeakValues-G' num2str(k) curSlave '.mat'],'sideBccPeakValues')
%                 save([dataPath filesep 'sideTimeToPeak-G' num2str(k) curSlave '.mat'],'sideTimeToPeak')
% 
%                 close(h2)
% 
%                 h2=figure;
%                 plot(mainTimeToPeak,mainBccPeakValues,'k.','MarkerSize',14)
%                 nameTitle=['mainTimePeak-Bcc-Class-' num2str(k)];
%                 title(nameTitle); xlabel('Time (s)'); ylabel('Bcc'); xlim('auto')
%                 hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%                 hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%                 close(h2)
% 
%                 mainTimeToPeakGroup{k,jj}=mainTimeToPeak;
%                 mainBccPeakValuesGroup{k,jj}=mainBccPeakValues;
%                 sideTimeToPeakGroup{k,jj}=sideTimeToPeak;
%                 sideBccPeakValuesGroup{k,jj}=sideBccPeakValues;
%             end
%             h2=figure;
%             plot(mainTimeToPeak,mainBccPeakValues,'k.','MarkerSize',14)
%             nameTitle=['mainTimePeak-Bcc-Class-' num2str(k)];
%             title(nameTitle); xlabel('Time (s)'); ylabel('Bcc'); xlim('auto')
%             hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%             hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%             close(h2)
%         elseif length(existingSlaveIDs)>1
            % Here I'll draw Bcc of each slave pair separately
%             curBccGroup = BccGroups(k,existingSlaveIDs);
%             if ~isempty(curBccGroup)
% %                 idGroupIndex=find(idGroups{k});
% %                 curTrackIndex = idGroupIndex(indexValidBccInTracks{k,jj});
%     %             curAvgBcc = avgBccGroups{k,jj};
%                 h2=figure;
%                 for jj=existingSlaveIDs
%                     curBcc = curBccGroup{jj};
%                     curSlave = potentialSlaves{jj};
% 
%                     tRange = (1:size(curBcc,2))-tShift;
%                     subplot(2,1,jj)
%                     plot(tRange,curBcc); hold on
%                     plot(tRange,nanmedian(curBcc,1),'k','Linewidth',3)
%                     line([36 36],[min(min(curBcc(:,:))),max(max(curBcc(:,:)))],'LineStyle','--','Color',[0.5 .5 .5],'linewidth',4)
%                     nameTitle=['Bcc ' num2str(k) 'th class, amp-' curSlave];
%                     title(nameTitle);                     
%                     xlabel('Time lag (s)')
%                     ylabel('Bcc (a.u.)')
%                 end
%                 hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%                 hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%                 close(h2)
%                 % Filtering done previously
%                 % Try clustering
%                 % #1 clustering with the features from my feature collection
%                 meas=extractGeneralFeatureTimeSeries(curBcc); %time of peak is the main feature
%                 numClusters = min(10,max(1,round(length(meas)/5))); 
%                 clusterArray=1:numClusters;
%                 T=kmeans(meas,numClusters);
% 
%                 % Sort T based on timeToMax
%                 meanTimeToMaxPerClus = arrayfun(@(x) mean(meas(T==x)),clusterArray);
%                 numPerClusters = arrayfun(@(x) sum(T==x),clusterArray);
%                 % Discard the cluster that has too early time to max
%                 [~,indClusterEarlyTimetoMax] = min(meanTimeToMaxPerClus);
%                 % Discard clusters whose sample sizes are small
%                 thresClusterSize = mean(numPerClusters)-0.5*std(numPerClusters);
%                 indSmallClusters = find(numPerClusters < thresClusterSize);
%                 % Choose the main cluster
%                 mainClusters = find(numPerClusters > max(numPerClusters)*0.95);
%                 sideClusters = setdiff(clusterArray,[mainClusters indSmallClusters indClusterEarlyTimetoMax]);
%                 % Found that They have very similar life time per cluster
%                 % Checking lifetime... They are not necessarily the same!
%                 % Inspect
%                 mainBccPeakValues=[];
%                 mainTimeToPeak=[];
%                 sideBccPeakValues=[];
%                 sideTimeToPeak=[];
%                 for curT=[mainClusters sideClusters]
%                     h2=figure; 
%                     for jj=existingSlaveIDs
%                         curBcc = curBccGroup{jj};
%                         slaveSource = potentialSlaves{jj};
%                         subplot(2,1,jj)
%                         plot((1:size(curBcc,2))-tShift,abs(curBcc(T==curT,:)),'o-')
%                         hold on
%                         plot((1:size(curBcc,2))-tShift,nanmedian(abs(curBcc(T==curT,:)),1),'Linewidth',3,'color','k')
%                         line([0 0],[min(min(curBcc(T==curT,:))),max(max(curBcc(T ==curT,:)))],'LineStyle','--','Color',[0.5 .5 .5],'linewidth',4)
%                         hold off
%                         if ismember(curT,mainClusters)
%                             nameTitle=['Amp-' slaveSource ' Bcc ' num2str(k) 'th class, ' num2str(curT) 'th cluster, main cluster'];
%                             mainBccPeakValues=[mainBccPeakValues; nanmax(curBcc(T==curT,:),[],2)];
%                             mainTimeToPeak = [mainTimeToPeak; (meas(T==curT)-tStartingFrame(k,jj))*tInterval];
%                         else
%                             nameTitle=['Amp-' slaveSource ' Bcc ' num2str(k) 'th class, ' num2str(curT) 'th cluster, side cluster'];
%                             sideBccPeakValues=[sideBccPeakValues; nanmax(curBcc(T==curT,:),[],2)];
%                             sideTimeToPeak = [sideTimeToPeak; (meas(T==curT)-tStartingFrame(k,jj))*tInterval];
%                         end
%                         title(nameTitle); 
%                         xlabel('Time (frame, 0th frame = time zero)')
%                     end
%                     hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
%                     hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
%                     print(h2,'-dtiff',strcat(tifPath,filesep,nameTitle,'.tif'))
% 
%                     curTind = find(T==curT);
%                     idGroupIndex=find(idGroups{k});
%                     curTrackIndex = idGroupIndex(indexValidBccInTracks{k,jj});
%                     iii=0;
%                     numTotalSubplots = 3+length(existingSlaveIDs);
%                     for ii=1:sum(T==curT)
%                         iii=iii+1;
%                         if iii==1
%                             h=figure('Position',[100,50,300,1200]);  
%                         end
%                         % corresponding index in tracksNA
% 
%                         curTrackInd = curTrackIndex((curTind(ii)));
%                         curTrack = tracksNA(curTrackInd);
%                         tShift = tStartingFrame(k,jj) + curTrack.startingFrameExtraExtra-curTrack.startingFrameExtra;
%                         tFluc = 10; % this is actually in frame
%                         lastFrameCC = curTrack.endingFrameExtraExtra-tFluc; 
%                         sd=curTrack.ampTotal;
%                         subplot(numTotalSubplots,1,1), hold off, plot(tShift-curTrack.startingFrameExtraExtra+(curTrack.startingFrameExtraExtra:lastFrameCC),sd(curTrack.startingFrameExtraExtra:lastFrameCC),'o-'), %hold on, plot(tRange(firstIncreaseTimeInt)-curTrack.startingFrameExtraExtra+1,sd(firstIncreaseTimeInt),'ro'); 
%                         title(['Fluorescent intensity (ii: ' num2str(ii) ')']); ylabel('F.I. (a.u.)')
%                         xlabel('Time (frame)')
%                         for jj=existingSlaveIDs
%                             curBcc = curBccGroup{jj};
%                             curSlave = potentialSlaves{jj};
%                         
%                             sCurForce_sd=curTrack.(curSlave);
%                             subplot(numTotalSubplots,1,1+(jj-1)*2+1), hold off, 
%                             plot(tShift-curTrack.startingFrameExtraExtra+(curTrack.startingFrameExtraExtra:lastFrameCC),sCurForce_sd(curTrack.startingFrameExtraExtra:lastFrameCC),'o-'), %hold on, plot(tRange(firstIncreaseTimeForce)-curTrack.startingFrameExtraExtra+1,sCurForce_sd(firstIncreaseTimeForce),'ro')    
%                             title([curSlave ' (ii: ' num2str(ii) ')']); ylabel(curSlave)
%                             xlabel('Time (frame)')
% 
%                             curCurBcc = curBcc(curTind(ii),:);
%                             subplot(numTotalSubplots,1,2+(jj-1)*2+1), plot(1:length(curCurBcc),curCurBcc,'o-'), hold on, 
%                             curName = [curSlave '-G-' num2str(k) '-Clst-' num2str(curT) '-numClus-' num2str(ii) '-tNum',num2str(curTrackInd)];
%                             title({curName; 'co-variance'}); ylabel('co-variance (a.u.)')
%                             xlabel('Time (frame)')
%                             hold off
%                         end
%                         hgsave(h,strcat(figPath,filesep,curName),'-v7.3')
%                         print(h,'-dtiff',strcat(tifPath,filesep,curName,'.tif'))
% %                         print(h,'-depsc','-loose',[epsPath filesep curName '.eps']);
%                         save(strcat(dataPath,filesep,curName,'.mat'),'sd','sCurForce_sd','curCurBcc')
% %                         uiwait(); 
%                     end
%                 end
%             end
%         end
%     end
%     if iii>0 && ishandle(h)
%         close(h)
%     end
%     save([dataPath filesep 'mainBccPeakValuesGroup.mat'],'mainBccPeakValuesGroup')
%     save([dataPath filesep 'mainTimeToPeakGroup.mat'],'mainTimeToPeakGroup')
%     save([dataPath filesep 'sideBccPeakValuesGroup.mat'],'sideBccPeakValuesGroup')
%     save([dataPath filesep 'sideTimeToPeakGroup.mat'],'sideTimeToPeakGroup')
%  
    
    h2=figure; ax = axes(h2); 
    boxPlotCellArray(peakLagTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['peakLag Class' num2str(k)];
    title(nameTitle); ylabel('Time lag (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'peakLagTogetherAdjusted','nameList2');    
    close(h2)

    h2=figure; ax = axes(h2);
    boxPlotCellArray(endingLagTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['endingLag Class' num2str(k)];
    title(nameTitle); ylabel('Time lag (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'endingLagTogetherAdjusted','nameList2');    
    close(h2)
    
    h2=figure; ax = axes(h2);
    boxPlotCellArray(ccLagTogetherAdjusted,nameList2,1,false,true,false,5,'ax',ax);
    nameTitle=['ccLag Class' num2str(k)];
    title(nameTitle); ylabel('Time lag (s)')
    hgexport(h2,strcat(figPath,filesep,nameTitle),hgexport('factorystyle'),'Format','eps')
    hgsave(h2,strcat(figPath,filesep,nameTitle),'-v7.3')
    save([dataPath filesep nameTitle],'ccLagTogetherAdjusted','nameList2');    
    close(h2)
end
%% V. All the other features
%% Festure statistics
    %% Look at feature difference per each group
    pixSize=MD.pixelSize_/1000; % in um
    distToEdge{1} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroups{1}));
    distToEdge{2} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroups{2}));
    distToEdge{3} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup3));
    distToEdge{4} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup4));
    distToEdge{5} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup5));
    distToEdge{6} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup6));
    distToEdge{7} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup7));
    distToEdge{8} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup8));
    distToEdge{9} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup9));
    groupNameList={'g1','g2','g3','g4','g5','g6','g7','g8','g9'};
    h2=figure; ax = axes(h2);
    boxPlotCellArray(distToEdge,groupNameList,pixSize,false,true,false,5,'ax',ax);
    title(ax,'Distance to edge')
    ylabel(ax,'Distance to edge (um)')
    save([p.OutputDirectory filesep 'data' filesep 'distToEdge.mat'],'distToEdge','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'distToEdgeForAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/distToEdgeForAllGroups'),'-v7.3'); close(h2)
    %% Look at feature difference per each group - advanceDist
    advanceDist{1} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroups{1}));
    advanceDist{2} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroups{2}));
    advanceDist{3} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup3));
    advanceDist{4} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup4));
    advanceDist{5} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup5));
    advanceDist{6} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup6));
    advanceDist{7} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup7));
    advanceDist{8} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup8));
    advanceDist{9} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup9));

    h2=figure; ax = axes(h2);
    boxPlotCellArray(advanceDist,groupNameList,1,false,true,false,5,'ax',ax);
    title(ax,'Adhesion advancement forward')
    ylabel(ax,'Adhesion advancement (um)')
    save([p.OutputDirectory filesep 'data' filesep 'advanceDist.mat'],'advanceDist','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'advanceDistAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/advanceDistAllGroups'),'-v7.3'); close(h2)
    %% Look at feature difference per each group - ampTotal
    ampTotal{1} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroups{1}));
    ampTotal{2} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroups{2}));
    ampTotal{3} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup3));
    ampTotal{4} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup4));
    ampTotal{5} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup5));
    ampTotal{6} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup6));
    ampTotal{7} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup7));
    ampTotal{8} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup8));
    ampTotal{9} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup9));
    h2=figure; ax = axes(h2);
    boxPlotCellArray(ampTotal,groupNameList,1,false,true,false,5,'ax',ax);

    title(ax,'ampTotal')
    ylabel(ax,'Fluorescence intensity (A.U.)')
    save([p.OutputDirectory filesep 'data' filesep 'ampTotal.mat'],'ampTotal','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'ampTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/ampTotalAllGroups'),'-v7.3'); close(h2)
    %% Look at feature difference per each group - assemRate
    assemRate{1} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroups{1}));
    assemRate{2} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroups{2}));
    assemRate{3} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup3));
    assemRate{4} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup4));
    assemRate{5} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup5));
    assemRate{6} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup6));
    assemRate{7} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup7));
    assemRate{8} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup8));
    assemRate{9} =arrayfun(@(x) nanmean(x.assemRate),tracksNA(idGroup9));
    h2=figure; ax = axes(h2);
    boxPlotCellArray(assemRate,groupNameList,1,false,true,false,5,'ax',ax);

    title(ax,'assemRate')
    ylabel(ax,'Assembly Rate(1/min)')
    save([p.OutputDirectory filesep 'data' filesep 'assemRate.mat'],'assemRate','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'assemRateAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/assemRateAllGroups'),'-v7.3'); close(h2)
    %% Look at feature difference per each group - assemRate
%     if ~isfield(tracksNA,'disassemRate')
%         tracksNA(end).disassemRate=[];
%     end
%     disp('Calculating disassembly rates'); tic
%     parfor kk=1:numTracks
%         curTrack=tracksNA(kk);
%         curDisassemRate = getDisassemRateFromTracks(curTrack,tIntMin);
%         curTrack.disassemRate = curDisassemRate;
%         tracksNA(kk)=curTrack;
% %         progressText(kk/numTracks,'Calculating disassembly rates') % Update text
%     end
%     toc
    disassemRate{1} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroups{1}));
    disassemRate{2} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroups{2}));
    disassemRate{3} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup3));
    disassemRate{4} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup4));
    disassemRate{5} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup5));
    disassemRate{6} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup6));
    disassemRate{7} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup7));
    disassemRate{8} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup8));
    disassemRate{9} =arrayfun(@(x) nanmean(x.disassemRate),tracksNA(idGroup9));
    h2=figure; ax = axes(h2);
    boxPlotCellArray(disassemRate,groupNameList,1,false,true,false,5,'ax',ax);

    title(ax,'disassemRate')
    ylabel(ax,'Disssembly Rate(1/min)')
    save([p.OutputDirectory filesep 'data' filesep 'disassemRate.mat'],'disassemRate','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'disassemRateAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/disassemRateAllGroups'),'-v7.3'); close(h2)
    %% Look at feature difference per each group - ampTotal2
    if isfield(tracksNA,'ampTotal2')
        ampTotal2{1} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroups{1}));
        ampTotal2{2} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroups{2}));
        ampTotal2{3} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup3));
        ampTotal2{4} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup4));
        ampTotal2{5} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup5));
        ampTotal2{6} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup6));
        ampTotal2{7} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup7));
        ampTotal2{8} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup8));
        ampTotal2{9} =arrayfun(@(x) nanmean(x.ampTotal2),tracksNA(idGroup9));
        h2=figure; ax = axes(h2);
        boxPlotCellArray(ampTotal2,groupNameList,1,false,true,false,5,'ax',ax);

        title(ax,'ampTotal2')
        ylabel(ax,'Fluorescence intensity (A.U.)')
        save([p.OutputDirectory filesep 'data' filesep 'ampTotal2.mat'],'ampTotal2','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'ampTotal2AllGroups.eps']);% histogramPeakLagVinVsTal -transparent
        hgsave(strcat(figPath,'/ampTotal2AllGroups'),'-v7.3'); 
        
        %% Assembly rate for ampTotal2
        if p.measureAmp2rates
            numTracks=numel(tracksNA);
            tIntMin=MD.timeInterval_/60;
            if ~isfield(tracksNA,'assemRate2')
                tracksNA(end).assemRate2=[];
            end
            if ~isfield(tracksNA,'disassemRate2')
                tracksNA(end).disassemRate2=[];
            end
            disp('Calculating assembly/disassembly rates'); tic
            parfor kk=1:numTracks
                curTrack=tracksNA(kk);
                curAssemRate = getAssemRateFromTracks(curTrack,tIntMin,'amp2');
                curTrack.assemRate2 = curAssemRate;
                curDisassemRate = getDisassemRateFromTracks(curTrack,tIntMin,'amp2');
                curTrack.disassemRate2 = curDisassemRate;
                tracksNA(kk)=curTrack;
            end
            toc
            
            assemRate2{1} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroups{1}));
            assemRate2{2} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroups{2}));
            assemRate2{3} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup3));
            assemRate2{4} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup4));
            assemRate2{5} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup5));
            assemRate2{6} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup6));
            assemRate2{7} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup7));
            assemRate2{8} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup8));
            assemRate2{9} =arrayfun(@(x) nanmean(x.assemRate2),tracksNA(idGroup9));
            h2=figure; ax = axes(h2);
            boxPlotCellArray(assemRate2,groupNameList,1,false,true,false,5,'ax',ax);

            title(ax,'assemRate2')
            ylabel(ax,'Assembly Rate 2 (1/min)')
            save([p.OutputDirectory filesep 'data' filesep 'assemRate2.mat'],'assemRate2','-v7.3')
            
            print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'assemRate2AllGroups.eps']);% histogramPeakLagVinVsTal -transparent
            hgsave(strcat(figPath,'/assemRate2AllGroups'),'-v7.3');
        %% disassembly rate for ampTotal2
            disassemRate2{1} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroups{1}));
            disassemRate2{2} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroups{2}));
            disassemRate2{3} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup3));
            disassemRate2{4} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup4));
            disassemRate2{5} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup5));
            disassemRate2{6} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup6));
            disassemRate2{7} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup7));
            disassemRate2{8} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup8));
            disassemRate2{9} =arrayfun(@(x) nanmean(x.disassemRate2),tracksNA(idGroup9));
            h2=figure; ax = axes(h2);
            boxPlotCellArray(disassemRate2,groupNameList,1,false,true,false,5,'ax',ax);

            title(ax,'disassemRate2')
            ylabel(ax,'Disssembly Rate 2 (1/min)')
            save([p.OutputDirectory filesep 'data' filesep 'disassemRate2.mat'],'disassemRate2','-v7.3')
            print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'disassemRate2AllGroups.eps']);% histogramPeakLagVinVsTal -transparent
            hgsave(strcat(figPath,'/disassemRate2AllGroups'),'-v7.3'); 
        end
    end
    %%
    if isfield(tracksNA,'ampTotal2')
        close(h2)
    end
    %% Look at feature difference per each group - starting ampTotal
    startingAmpTotal{1} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroups{1}));
    startingAmpTotal{2} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroups{2}));
    startingAmpTotal{3} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup3));
    startingAmpTotal{4} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup4));
    startingAmpTotal{5} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup5));
    startingAmpTotal{6} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup6));
    startingAmpTotal{7} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup7));
    startingAmpTotal{8} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup8));
    startingAmpTotal{9} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup9));
    hFig=figure; ax=axes(hFig);
    boxPlotCellArray(startingAmpTotal,groupNameList,1,false,true,false,5,'ax',ax);
    title(ax,'Starting Amplitude')
    ylabel(ax,'Fluorescence intensity (A.U.)')
    save([p.OutputDirectory filesep 'data' filesep 'startingAmpTotal.mat'],'startingAmpTotal','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'startingAmpTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/startingAmpTotalAllGroups'),'-v7.3'); close(hFig)
    %% Look at feature difference per each group - starting ampTotal2
    if isfield(tracksNA,'ampTotal2')
        startingAmpTotal2{1} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroups{1}));
        startingAmpTotal2{2} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroups{2}));
        startingAmpTotal2{3} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup3));
        startingAmpTotal2{4} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup4));
        startingAmpTotal2{5} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup5));
        startingAmpTotal2{6} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup6));
        startingAmpTotal2{7} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup7));
        startingAmpTotal2{8} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup8));
        startingAmpTotal2{9} =arrayfun(@(x) x.ampTotal2(x.startingFrameExtra),tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(startingAmpTotal2,groupNameList,1,false,true,false,5,'ax',ax);

        title(ax,'startingAmpTotal2')
        ylabel(ax,'Fluorescence intensity (A.U.)')
        save([p.OutputDirectory filesep 'data' filesep 'startingAmpTotal2.mat'],'startingAmpTotal2','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'startingAmpTotal2.eps']);% histogramPeakLagVinVsTal -transparent
        hgsave(strcat(figPath,'/startingAmpTotal2'),'-v7.3'); 
    end
    %%
    if isfield(tracksNA,'ampTotal2')
        close(hFig)
    end
   
    %% Look at feature difference per each group - starting edgeAdvanceDistChange
    edgeAdvanceDistChange{1} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroups{1}));
    edgeAdvanceDistChange{2} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroups{2}));
    edgeAdvanceDistChange{3} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup3));
    edgeAdvanceDistChange{4} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup4));
    edgeAdvanceDistChange{5} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup5));
    edgeAdvanceDistChange{6} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup6));
    edgeAdvanceDistChange{7} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup7));
    edgeAdvanceDistChange{8} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup8));
    edgeAdvanceDistChange{9} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup9));
    hFig=figure; ax=axes(hFig);
    boxPlotCellArray(edgeAdvanceDistChange,groupNameList,1,false,true,false,5,'ax',ax);
    title(ax,'edgeAdvanceDistChange at the end of tracks (to see g7 has nearly zero edge advance)')
    ylabel(ax,'edgeAdvanceDistChange (um)')
    save([p.OutputDirectory filesep 'data' filesep 'edgeAdvanceDistChange.mat'],'edgeAdvanceDistChange','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'edgeAdvanceDistChangeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/edgeAdvanceDistChangeAllGroups'),'-v7.3'); close(hFig)
    %% Look at feature difference per each group - starting forceMag
    if isfield(tracksNA,'forceMag')  
        startingForceMag{1} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroups{1}));
        startingForceMag{2} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroups{2}));
        startingForceMag{3} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup3));
        startingForceMag{4} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup4));
        startingForceMag{5} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup5));
        startingForceMag{6} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup6));
        startingForceMag{7} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup7));
        startingForceMag{8} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup8));
        startingForceMag{9} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(startingForceMag,groupNameList,1,false,true,false,5,'ax',ax);
        title(ax,'startingForceMag (to see g1,2,3,7 start nearly at similar force)')
        ylabel(ax,'startingForceMag (Pa)')
        save([p.OutputDirectory filesep 'data' filesep 'startingForceMag.mat'],'startingForceMag','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'startingForceMagAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
        hgsave(strcat(figPath,'/startingForceMagAllGroups'),'-v7.3'); 
        close(hFig)
    end
    %%
    %% Look at feature difference per each group - mainTimeToPeakGroup
%     nameTwoGroups={'G1','G2'};
%     numIdGroups=cellfun(@sum,idGroups);
%     nameTwoGroups(numIdGroups(1:2)<1)=[];
%     for jj=existingSlaveIDs    
%         curSlaveCell = potentialSlaves(jj);
%         curSlave=curSlaveCell{1};
%         
%         figure;
%         boxPlotCellArray(mainTimeToPeakGroup(:,jj)',nameTwoGroups);
%         title(['mainTimeToPeakBccGroup, amp-' curSlave])
%         ylabel('mainTimeToPeakBccGroup (sec)')
%         print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'mainTimeToPeakBccGroup' curSlave '.eps']);% histogramPeakLagVinVsTal -transparent
%         hgsave(strcat(figPath,'/mainTimeToPeakBccGroup', curSlave),'-v7.3'); 
% %         %% Look at feature difference per each group - mainBccPeakValuesGroup
% %         figure;
% %         boxPlotCellArray(mainBccPeakValuesGroup(:,jj)',nameTwoGroups);
% %         title(['mainBccPeakValuesGroup, amp-' curSlave])
% %         ylabel('mainBccPeakValuesGroup (AU)')
% %         print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'mainBccPeakValuesGroup' curSlave '.eps']);% histogramPeakLagVinVsTal -transparent
% %         hgsave(strcat(figPath,'/mainBccPeakValuesGroup', curSlave),'-v7.3'); 
% %         %% Look at feature difference per each group - mainBccPeakValuesGroup
% %         figure;
% %         boxPlotCellArray(sideBccPeakValuesGroup(:,jj)',nameTwoGroups);
% %         title(['sideBccPeakValuesGroup, amp-' curSlave])
% %         ylabel('sideBccPeakValuesGroup (AU)')
% %         print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'sideBccPeakValuesGroup' curSlave '.eps']);% histogramPeakLagVinVsTal -transparent
% %         hgsave(strcat(figPath,'/sideBccPeakValuesGroup', curSlave),'-v7.3'); 
% %         %% Look at feature difference per each group - mainBccPeakValuesGroup
% %         figure;
% %         boxPlotCellArray(sideTimeToPeakGroup(:,jj)',nameTwoGroups);
% %         title(['sideTimeToPeakGroup, amp-' curSlave])
% %         ylabel('sideTimeToPeakGroup (AU)')
% %         print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'sideTimeToPeakGroup' curSlave '.eps']);% histogramPeakLagVinVsTal -transparent
% %         hgsave(strcat(figPath,'/sideTimeToPeakGroup', curSlave),'-v7.3'); 
%     end
    %% recalculate force slope 
    disp('Calculating slopes of amplitude and force again...')
    tic
    tracksNA = calculateTrackSlopes(tracksNA,tInterval);
    toc
    %% assembly rate and force growth rate in curIndices
    if isfield(tracksNA,'forceMag')  
        forceSlopeG1 =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup1f));
        earlyAmpSlopeG1 =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup1f));

        save([p.OutputDirectory filesep 'data' filesep 'assemblyRateForceSlopes.mat'],'forceSlopeG1','earlyAmpSlopeG1')
    end
    %% Look at feature difference per each group - earlyAmpSlope
    earlyAmpSlope{1} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroups{1}));
    earlyAmpSlope{2} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroups{2}));
    earlyAmpSlope{3} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup3));
    earlyAmpSlope{4} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup4));
    earlyAmpSlope{5} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup5));
    earlyAmpSlope{6} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup6));
    earlyAmpSlope{7} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup7));
    earlyAmpSlope{8} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup8));
    earlyAmpSlope{9} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup9));
    %     [lengthLongestSlope]=max(cellfun(@(x) length(x),earlyAmpSlope));
    % 
    %     matrixEarlyAmpSlope = NaN(lengthLongestSlope,9);
    %     for ii=1:9
    %         matrixEarlyAmpSlope(1:length(earlyAmpSlope{ii}),ii) = earlyAmpSlope{ii};
    %     end
    %     boxWidth=0.5;
    %     figure
    %     boxplot(matrixEarlyAmpSlope,'orientation','vertical','whisker',0.5,'notch','on',...
    %         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    %     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    %     set(findobj(gca,'tag','Median'),'LineWidth',2)
    % %     ylim([-2 50])
    hFig=figure; ax=axes(hFig);
    boxPlotCellArray(earlyAmpSlope,groupNameList,1,false,true,false,5,'ax',ax);
    title(ax,'earlyAmpSlope (this shows that g7 has the same early slope as g3).')
    ylabel(ax,'earlyAmpSlope (A.U./min)')
    hgsave(strcat(figPath,'/earlyAmpSlopeAllGroups'),'-v7.3')
    save([p.OutputDirectory filesep 'data' filesep 'earlyAmpSlopeAllGroups.mat'],'earlyAmpSlope','-v7.3')
    print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'earlyAmpSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %%
    close(hFig)
    %% Look at feature difference per each group - earlyAmpSlope2
    if isfield(tracksNA,'ampTotal2')
        earlyAmpSlope2{1} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroups{1}));
        earlyAmpSlope2{2} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroups{2}));
        earlyAmpSlope2{3} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup3));
        earlyAmpSlope2{4} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup4));
        earlyAmpSlope2{5} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup5));
        earlyAmpSlope2{6} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup6));
        earlyAmpSlope2{7} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup7));
        earlyAmpSlope2{8} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup8));
        earlyAmpSlope2{9} =arrayfun(@(x) x.earlyAmpSlope2,tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(earlyAmpSlope2,groupNameList,1,false,true,false,5,'ax',ax);

        title(ax,'earlyAmpSlope2')
        ylabel(ax,'AmpSlope (A.U./min)')
        save([p.OutputDirectory filesep 'data' filesep 'earlyAmpSlope2.mat'],'earlyAmpSlope2','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'earlyAmpSlope2.eps']);% histogramPeakLagVinVsTal -transparent
        hgsave(strcat(figPath,'/earlyAmpSlope2'),'-v7.3'); 
    end
    %% Look at feature difference per each group - ampSlope2
    if isfield(tracksNA,'ampTotal2')
        ampSlope2{1} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroups{1}));
        ampSlope2{2} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroups{2}));
        ampSlope2{3} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup3));
        ampSlope2{4} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup4));
        ampSlope2{5} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup5));
        ampSlope2{6} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup6));
        ampSlope2{7} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup7));
        ampSlope2{8} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup8));
        ampSlope2{9} =arrayfun(@(x) x.ampSlope2,tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(ampSlope2,groupNameList,1,false,true,false,5,'ax',ax);

        title(ax,'ampSlope2')
        ylabel(ax,'AmpSlope (A.U./min)')
        save([p.OutputDirectory filesep 'data' filesep 'ampSlope2.mat'],'ampSlope2','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'ampSlope2.eps']);% histogramPeakLagVinVsTal -transparent
        hgsave(strcat(figPath,'/ampSlope2'),'-v7.3'); 
    end
    %%
    if isfield(tracksNA,'ampTotal2')
        close(hFig)
    end
    %% Look at feature difference per each group - force slope
    if isfield(tracksNA,'forceMag')  
        forceSlope{1} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroups{1}));
        forceSlope{2} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroups{2}));
        forceSlope{3} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup3));
        forceSlope{4} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup4));
        forceSlope{5} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup5));
        forceSlope{6} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup6));
        forceSlope{7} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup7));
        forceSlope{8} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup8));
        forceSlope{9} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(forceSlope,groupNameList,1,false,true,false,5,'ax',ax);
        title(ax,'forceSlope (all nascent adhesions (g1,2,3,7) show the same force slopes).')
        ylabel(ax,'forceSlope (Pa/min)')
        hgsave(strcat(figPath,'/forceSlopeAllGroups'),'-v7.3')
        save([p.OutputDirectory filesep 'data' filesep 'forceSlopeAllGroups.mat'],'forceSlope','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'forceSlopeAllGroups.eps']); 
        %%
        close(hFig)
        %% Look at feature difference per each group - force slope
        earlyForceSlope{1} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroups{1}));
        earlyForceSlope{2} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroups{2}));
        earlyForceSlope{3} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup3));
        earlyForceSlope{4} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup4));
        earlyForceSlope{5} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup5));
        earlyForceSlope{6} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup6));
        earlyForceSlope{7} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup7));
        earlyForceSlope{8} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup8));
        earlyForceSlope{9} =arrayfun(@(x) (x.earlyForceSlope),tracksNA(idGroup9));
        hFig=figure; ax=axes(hFig);
        boxPlotCellArray(earlyForceSlope,groupNameList,1,false,true,false,5,'ax',ax);
        title(ax,'earlyForceSlope')
        ylabel(ax,'earlyForceSlope (Pa/min)')
        hgsave(strcat(figPath,'/earlyForceSlopeAllGroups'),'-v7.3')
        save([p.OutputDirectory filesep 'data' filesep 'earlyForceSlopeAllGroups.mat'],'earlyForceSlope','-v7.3')
        print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'earlyForceSlopeAllGroups.eps']); 
        %%
        close(hFig)
    end
%% Distributing to each group (after filtering)
    if isfield(tracksNA,'forceMag')  
        %% drawing group1
        fileStore = [epsPath filesep 'ampForcePlotG1.eps'];
        %     plotIntensityForce(tracksNA(idGroup1f),fileStore,false,false)
        [~,h]=plotIntensityForce(tracksNA(idGroup1f),fileStore,false,false); if ~isempty(h); close(h); end
        fileStore = [epsPath filesep 'ampForcePlotG1cohorts.eps'];
        [~,h]=plotIntensityForce(tracksNA(idGroup1f),fileStore,false,false,'plotCohorts',true); if ~isempty(h); close(h); end
        fileStore = [epsPath filesep 'ampForcePlotG1color.eps'];
        [~,h]=plotIntensityForce(tracksNA(idGroup1f),fileStore,false,true); if ~isempty(h); close(h); end

        %% group 2
        fileStoreG2 = [epsPath filesep 'ampForcePlotG2.eps'];
        [~,h]= plotIntensityForce(tracksNA(idGroup2f),fileStoreG2,false,true); if ~isempty(h); close(h); end
        fileStoreG2 = [epsPath filesep 'ampForcePlotG2cohorts.eps'];
        [~,h]= plotIntensityForce(tracksNA(idGroup2f),fileStoreG2,false,true,'plotCohorts',true); if ~isempty(h); close(h); end
        %% group 3 plotting
        fileStoreG3 = [epsPath filesep 'ampForcePlotG3.eps'];
        [~,h]=plotIntensityForce(tracksNA(idGroup3),fileStoreG3,false,false); if ~isempty(h); close(h); end
        %% group4 plotting
        fileStoreG4 = [epsPath filesep 'ampForcePlotG4.eps'];
        [~,h]=plotIntensityForce(tracksNA(idGroup4),fileStoreG4,false,false); if ~isempty(h); close(h); end
    %     %% group5 plotting
    %     fileStoreG5 = [epsPath filesep 'ampForcePlotG5.eps'];
    %     [~,h]=plotIntensityForce(tracksNA(idGroup5),fileStoreG5,false,false); if ~isempty(h); close(h); end
    %     %% group6 plotting
    %     fileStoreG6 = [epsPath filesep 'ampForcePlotG6.eps'];
    %     [~,h]=plotIntensityForce(tracksNA(idGroup6),fileStoreG6,false,false); if ~isempty(h); close(h); end
    %     %% group7 plotting
    %     fileStoreG7 = [epsPath filesep 'ampForcePlotG7.eps'];
    %     [~,h]=plotIntensityForce(tracksNA(idGroup7),fileStoreG7,false,false); if ~isempty(h); close(h); end
    %     %% group8 plotting
    %     fileStoreG8 = [epsPath filesep 'ampForcePlotG8.eps'];
    %     [~,h]=plotIntensityForce(tracksNA(idGroup8),fileStoreG8,false,false); if ~isempty(h); close(h); end
    %     %% group9 plotting
    %     fileStoreG9 = [epsPath filesep 'ampForcePlotG9.eps'];
    %     [~,h]=plotIntensityForce(tracksNA(idGroup9),fileStoreG9,false,false); if ~isempty(h); close(h); end
    end
%% G3 vs. G7 comparison
    %% export tracksG1, G2, G3 and G7 separately
%     tracksG1 = tracksNA(idGroup1);
%     tracksG2 = tracksNA(idGroup2f);
%     tracksG3 = tracksNA(idGroup3);
%     tracksG7 = tracksNA(idGroup7);
%     save([p.OutputDirectory filesep 'data' filesep 'tracksG1.mat'],'tracksG1','-v7.3')
%     save([p.OutputDirectory filesep 'data' filesep 'tracksG2.mat'],'tracksG2','-v7.3')
%     save([p.OutputDirectory filesep 'data' filesep 'tracksG3.mat'],'tracksG3','-v7.3')
%     save([p.OutputDirectory filesep 'data' filesep 'tracksG7.mat'],'tracksG7','-v7.3')
%     %% Filtering G7 for those that increase edge advance 
%     timeInterval = MD.timeInterval_/60; % in min
%     numG7=numel(tracksG7);
%     risingAndStoppingG7=[];
%     edgeRisingStoppingG7{1}=[];
%     adhAdvanceRisingStoppingG7{1}=[];
%     divPointRS=[];
%     tractionRisingStoppingG7{1}=[];
%     risingSlopeG7=[];
%     forceIncrease1Half=[];
%     forceIncrease2Half=[];
%     edgeIncrease1Half=[];
%     edgeIncrease2Half=[];
%     adhIncrease1Half=[];
%     adhIncrease2Half=[];
%     p=0;
%     q=0;
%     for kk=1:numG7
%         curG7=tracksG7(kk);
%         if curG7.lifeTime<5
%             continue
%         end
%         tRange = 1:length(curG7.edgeAdvanceDist);
%         if curG7.startingFrameExtra>1
%             curG7.edgeAdvanceDist(1:curG7.startingFrameExtra-1)=NaN;
%             curG7.forceMag(1:curG7.startingFrameExtra-1)=NaN;
%             curG7.advanceDist(1:curG7.startingFrameExtra-1)=NaN;
%         end
%         % fit with partial sloped line + flat line in the later
%         x=(0:(sum(~isnan(curG7.edgeAdvanceDist))-1))*timeInterval;
%         y=curG7.edgeAdvanceDist(~isnan(curG7.edgeAdvanceDist));
%         fun=@(z) [(z(1)*x(x<=z(3))+z(2)) (z(1)*z(3)+z(2))*ones(1,sum(x>z(3)))]-y;
%         z0=[0.1,0,20*timeInterval];
%     %     [x0out,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0i);
%         z=lsqnonlin(fun,z0);
%         % Since the purpose is to separate out those that increase about a
%         % half and stay the rest of half, I'll filter out a<=0,
%         % xt/lifetime>0.6.
%         if z(1)>0.1 && z(3)/x(end)<0.8 && z(3)/x(end)>0.2
%             p=p+1;
%             q=q+1;
%             risingAndStoppingG7=[risingAndStoppingG7 kk];
%             curG7Advance=curG7.advanceDist(~isnan(curG7.advanceDist));
%             figure, subplot(2,1,1),plot(x,curG7Advance,'o--','Color',[1,140/255,0]), hold on,
%             plot(x,y,'ko--'), 
%             plot(x,[(z(1)*x(x<=z(3))+z(2)) (z(1)*z(3)+z(2))*ones(1,sum(x>z(3)))],'k-','Linewidth',2)
%             title(['Edge Protrusion, ID:' num2str(kk)]); ylabel('Protrusion distance (um)')
%             xt = round(z(3)/timeInterval);
%             edgeRisingStoppingG7{p}=y;
%             adhAdvanceRisingStoppingG7{p}=curG7Advance;
%             divPointRS=[divPointRS xt];
%             risingSlopeG7=[risingSlopeG7 z(1)];
% 
%             try
%                 sForce= csaps(tRange,curG7.forceMag,splineParam);
%             catch
%                 curG7.forceMag(curG7.endingFrameExtra+1:end)=[];
%                 sForce= csaps(tRange,curG7.forceMag,splineParam);
%             end
%             sdForce=ppval(sForce,tRange);
%             sdForce(isnan(curG7.forceMag))=NaN;
%             curG7Force = sdForce(~isnan(curG7.edgeAdvanceDist));
%             curG7ForceRaw = curG7.forceMag(~isnan(curG7.edgeAdvanceDist));
%             tractionRisingStoppingG7{p}=curG7ForceRaw;
%         %     [~,m_first]=regression(x(1:xt),y(1:xt))
%         %     [~,m_second]=regression(x(xt:end),y(xt:end))
%         %     [~,m_first]=regression(x(1:xt),curG7Force(1:xt))
%         %     [~,m_second]=regression(x(xt:end),curG7Force(xt:end))
%             forceIncrease1Half(p)=(max(curG7ForceRaw(xt-3:xt))-min(curG7ForceRaw(1:3)))/z(3);
%             forceIncrease2Half(q)=(max(curG7ForceRaw(end-5:end))-min(curG7ForceRaw(xt-1:xt+1)))/(x(end)-z(3));
%             edgeIncrease1Half(p)=(max(y(xt-3:xt))-min(y(1:3)))/z(3);
%             edgeIncrease2Half(q)=(max(y(end-5:end))-min(y(xt-1:xt+1)))/(x(end)-z(3));
%             adhIncrease1Half(p)=(max(curG7Advance(xt-3:xt))-min(curG7Advance(1:3)))/z(3);
%             adhIncrease2Half(q)=(max(curG7Advance(end-5:end))-min(curG7Advance(xt-1:xt+1)))/(x(end)-z(3));
%             subplot(2,1,2),plot(x,curG7ForceRaw,'ko-'),hold on,plot(x,curG7Force,'k-','LineWidth',2)
%             title('Traction'); ylabel('Traction Magnitude (Pa)'); xlabel('Time (min)');
%             hgsave(strcat(figPath,'/EdgeAndForceInG7',num2str(kk)),'-v7.3')
%             print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'EdgeAndForceInG7' num2str(kk) '.eps']);% histogramPeakLagVinVsTal -transparent
%         else
%             q=q+1;
%             curG7ForceRaw = curG7.forceMag(~isnan(curG7.edgeAdvanceDist));
%             forceIncrease2Half(q)=(max(curG7ForceRaw(end-5:end))-min(curG7ForceRaw(1:3)))/(x(end)-x(1));
%         end
%     end
%     save(strcat(dataPath,'/EdgeAndForceInG7selected.mat'),'risingAndStoppingG7','risingSlopeG7',...
%         'forceIncrease1Half','forceIncrease2Half','edgeIncrease1Half','edgeIncrease2Half','adhIncrease1Half','adhIncrease2Half',...
%         'edgeRisingStoppingG7','adhAdvanceRisingStoppingG7','divPointRS','tractionRisingStoppingG7','-v7.3')
%     %% Filtering G3 for those that increase edge advance 
%     numG3=numel(tracksG3);
%     risingAndStoppingG3=[];
%     edgeRisingStoppingG3{1}=[];
%     adhAdvanceRisingStoppingG3{1}=[];
%     forceIncreaseG3=[];
%     edgeIncreaseG3=[];
%     adhIncreaseG3=[];
%     tractionRisingStoppingG3{1}=[];
%     risingSlopeG3=[];
%     p=0;
%     for kk=1:numG3
%         curG3=tracksG3(kk);
%         if curG3.lifeTime<5
%             continue
%         end
%         tRange = 1:length(curG3.edgeAdvanceDist);
%         if curG3.startingFrameExtra>1
%             curG3.edgeAdvanceDist(1:curG3.startingFrameExtra-1)=NaN;
%             curG3.forceMag(1:curG3.startingFrameExtra-1)=NaN;
%             curG3.advanceDist(1:curG3.startingFrameExtra-1)=NaN;
%         end
%         % fit with partial sloped line + flat line in the later
%         x=(0:(sum(~isnan(curG3.edgeAdvanceDist))-1))*timeInterval;
%         y=curG3.edgeAdvanceDist(~isnan(curG3.edgeAdvanceDist));
%         fun=@(z) (z(1)*x+z(2))-y;
%         z0=[0.1,0];
%     %     [x0out,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0i);
%         z=lsqnonlin(fun,z0);
%         % Since the purpose is to separate out those that increase about a
%         % half and stay the rest of half, I'll filter out a<=0,
%         % xt/lifetime>0.6.
%         if z(1)>0.1
%             p=p+1;
%             risingAndStoppingG3=[risingAndStoppingG3 kk];
%             curG3Advance=curG3.advanceDist(~isnan(curG3.advanceDist));
%             figure, subplot(2,1,1),plot(x,curG3Advance,'o--','Color',[1,140/255,0]), hold on,
%             plot(x,y,'ko--'), 
%             plot(x,(z(1)*x+z(2)),'k-','Linewidth',2)
%             title(['Edge Protrusion, ID:' num2str(kk)]); ylabel('Protrusion distance (um)')
%             edgeRisingStoppingG3{p}=y;
%             adhAdvanceRisingStoppingG3{p}=curG3Advance;
%             risingSlopeG3=[risingSlopeG3 z(1)];
% 
%             try
%                 sForce= csaps(tRange,curG3.forceMag,splineParam);
%             catch
%                 curG3.forceMag(curG3.endingFrameExtra+1:end)=[];
%                 sForce= csaps(tRange,curG3.forceMag,splineParam);
%             end
%             sdForce=ppval(sForce,tRange);
%             sdForce(isnan(curG3.forceMag))=NaN;
%             curG3Force = sdForce(~isnan(curG3.edgeAdvanceDist));
%             curG3ForceRaw = curG3.forceMag(~isnan(curG3.edgeAdvanceDist));
%             tractionRisingStoppingG3{p}=curG3ForceRaw;
% 
%             forceIncreaseG3(p)=(max(curG3Force(end-2:end))-min(curG3Force(1:3)))/x(end);
%             edgeIncreaseG3(p)=(max(y(end-2:end))-min(y(1:3)))/x(end);
%             adhIncreaseG3(p)=(max(curG3Advance(end-2:end))-min(curG3Advance(1:3)))/x(end);
%             subplot(2,1,2),plot(x,curG3ForceRaw,'ko-'),hold on,plot(x,curG3Force,'k-','LineWidth',2)
%             title('Traction'); ylabel('Traction Magnitude (Pa)'); xlabel('Time (min)');
%             hgsave(strcat(figPath,'/EdgeAndForceInG3',num2str(kk)),'-v7.3')
%             print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'EdgeAndForceInG3' num2str(kk) '.eps']);% histogramPeakLagVinVsTal -transparent
%         end
%     end
%     save(strcat(dataPath,'/EdgeAndForceInG3selected.mat'),'risingAndStoppingG3','edgeRisingStoppingG3',...
%         'adhAdvanceRisingStoppingG3','forceIncreaseG3','edgeIncreaseG3','adhIncreaseG3','tractionRisingStoppingG3','-v7.3')
%     %% Last, compare force increase among forceIncreaseG3,forceIncrease1Half and forceIncrease2Half
%     forceIncreaseG3G7Cell={ forceIncreaseG3,forceIncrease1Half, forceIncrease2Half};
%     nameG3G7={'G3', 'G7-protruding phase', 'G7-stalling phase'};
%     figure;
%     barPlotCellArray(forceIncreaseG3G7Cell,nameG3G7)
%     title('Increase in force in G3 and two separate phases in G7')
%     ylabel('Change in force (Pa/min)')
%     hgsave(strcat(figPath,'/forceIncreaseG3G7'),'-v7.3')
%     save([p.OutputDirectory filesep 'data' filesep 'forceIncreaseG3G7.mat'],'forceIncreaseG3G7Cell','nameG3G7','-v7.3')
%     print('-depsc','-loose',[p.OutputDirectory filesep 'eps' filesep 'forceIncreaseG3G7.eps']);% histogramPeakLagVinVsTal -transparent
% 
%     edgeAdhIncreaseG3G7Cell={ edgeIncreaseG3,adhIncreaseG3,edgeIncrease1Half,adhIncrease1Half,...
%         edgeIncrease2Half,adhIncrease2Half};
%     nameG3G7edgeAdh={'G3-edge','Ge-adh', 'G7-edge-protruding phase','G7-adh-protruding phase',...
%         'G7-edge-stalling phase','G7-adh-stalling phase'};
%     figure;
%     barPlotCellArray(edgeAdhIncreaseG3G7Cell,nameG3G7edgeAdh,pixSize)
%     title('Advance in edge and adhesion in G3 and two separate phases in G7')
%     ylabel('Change in edge (um/min)')
%     hgsave(strcat(figPath,'/edgeAdhIncreaseG3G7'),'-v7.3')
%     save([p.OutputDirectory filesep 'data' filesep 'edgeAdhIncreaseG3G7.mat'],'edgeAdhIncreaseG3G7Cell','nameG3G7edgeAdh','-v7.3')


%% Saving
% disp('Saving...')
% save(outputFile{1,p.ChannelIndex},'tracksNA','-v7.3'); % the later channel has the most information.
disp('Initial Rise Time Lab Process Done!')
end
%
