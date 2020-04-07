function []=classifyMovieNascentAdhesions(MD)
% 
% this function classifyMoiveNascentAdhesions(MD) reads from analyzeAdhesionMaturation 
% and classify NA tracks based on fluorescence signal amplitude evolution
% and distance to the edge. 
% group 1: NAs that form at the edge but say there as the edge protrude and
% turn over:
% group 2: NAs that form and mature as the edge protrude: 
% group 3: NAs that move along the edge 

% This function does not rely on colocalizationAdhesionWithTFM because it
% should basically not require force value for classification.

%% input reading
ip = inputParser;
ip.addRequired('MD', @(x)(isa(x,'MovieData')));
ip.parse(MD);
%Get the indices of any previous processes
iProc = MD.getProcessIndex('AdhesionClassificationProcess', 'nDesired', 1, 'askUser', false);
%Check if process exists
if isempty(iProc)
    iProc = numel(MD.processes_)+1;
    MD.addProcess(AdhesionClassificationProcess(MD,MD.outputDirectory_));                                                                                                 
end
%Parse input, store in parameter structure
adhClassProc = MD.processes_{iProc};
p = parseProcessParams(adhClassProc);
p.useHomogeneity = false;
% p.startingDist = 3; % in micron

MD=ip.Results.MD;
iChan = p.ChannelIndex;
sampleFolders=p.labelData;
%% Load processed data
nFrames=MD.nFrames_;

iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
adhAnalProc = MD.getProcess(iAdhProc);
% numChans = numel(p.ChannelIndex);
% tracksNA = load(adhAnalProc.outFilePaths_{1,iChan},'tracksNA');
% tracksNA = tracksNA.tracksNA;

% try
    tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA');
% catch
%     % Check if the outFilePath has tableTracksNA
%     disp('Checking if the outFilePath has tracksNA...')
%     s = load(adhAnalProc.outFilePaths_{1,iChan},'tracksNA');
%     if isfield(s,'tracksNA')
%         disp('Found the old format. Resaving this with the new format...')
%         % Saving with each track
%         tracksNA = s.tracksNA;
%         trackFolderPath = [adhAnalProc.funParams_.OutputDirectory filesep 'trackIndividual'];
%         mkdir(trackFolderPath)
%         numTracks = numel(tracksNA);
%         fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
%         numStr = @(trackNum) num2str(trackNum,fString);
%         trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
% 
%         for ii=1:numTracks
%             curTrack = tracksNA(ii);
%             if iscell(curTrack.state)
%                 curTrack.state = strcmp(curTrack.state,'BA')+2*strcmp(curTrack.state,'NA')+...
%                     3*strcmp(curTrack.state,'FC')+4*strcmp(curTrack.state,'FA')+...
%                     5*strcmp(curTrack.state,'ANA')+6*strcmp(curTrack.state,'Out_of_Band');
%             end
%             save(trackIndPath(ii),'curTrack')
%             progressText((ii)/numTracks,'Saving individual tracksNA') % Update text
%         end
%         % Saving the metaTrackData
%         metaTrackData.numTracks = numTracks;
%         metaTrackData.trackFolderPath = trackFolderPath;
%         metaTrackData.eachTrackName = 'curTrack';
%         metaTrackData.fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
%         metaTrackData.numStr = @(trackNum) num2str(trackNum,fString);
%         metaTrackData.trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
%         save(adhAnalProc.outFilePaths_{1,iChan},'metaTrackData')
%     end
% end

numTracks = numel(tracksNA);
if iscell(tracksNA(1).state) %in case the tracksNA.state is in cell format
    for ii=1:numTracks
        curState = tracksNA(ii).state;
        curState = strcmp(curState,'BA')+2*strcmp(curState,'NA')+...
                3*strcmp(curState,'FC')+4*strcmp(curState,'FA')+...
                5*strcmp(curState,'ANA')+6*strcmp(curState,'Out_of_Band');
        tracksNA(ii).state = curState;
    end
end

% This case SDC was not used and first img frame was used.
paxImage=MD.getChannel(iChan).loadImage(1); 
[h,w] = size(paxImage);

imgMap = zeros(h,w,nFrames);
for ii=1:nFrames
    paxImage=MD.getChannel(iChan).loadImage(ii); 
    imgMap(:,:,ii) = paxImage;
end
%% Backup previous Analysis output - This created too much memory. Deleting...
% if exist(p.OutputDirectory,'dir')
%     if p.backupOldResults
%         disp('Backing up the original data')
%         ii = 1;
%         backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
%         while exist(backupFolder,'dir')
%             backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
%             ii=ii+1;
%         end
%         copyfile(p.OutputDirectory, backupFolder,'f')
%     end
% end
mkdir(p.OutputDirectory); %mkClrDir(p.OutputDirectory);

%% Output setup
outputPath = [p.OutputDirectory filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
dataPath = [p.OutputDirectory filesep 'data'];
if ~exist(dataPath,'dir')
    mkdir(dataPath);
end
epsPath=[p.OutputDirectory filesep 'eps'];
tifPath=[p.OutputDirectory filesep 'tif'];
figPath=[p.OutputDirectory filesep 'figs'];
if ~exist(epsPath,'dir')
    mkdir(epsPath);
end
if ~exist(tifPath,'dir')
    mkdir(tifPath);
end
if ~exist(figPath,'dir')
    mkdir(figPath);
end

% outputFile = cell(5, 1);
% outputFile{1,1} = [p.OutputDirectory filesep 'data' filesep 'selectedGroups.mat'];
% outputFile{2,1} = [p.OutputDirectory filesep 'figs' filesep 'featureSpace.fig'];
% outputFile{3,1} = [p.OutputDirectory filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified.fig'];
% outputFile{4,1} = [p.OutputDirectory filesep 'data' filesep 'idsClassified.mat'];
% outputFile{5,1} = [p.OutputDirectory filesep 'data' filesep 'tracksNA.mat'];

% Set up the output files
outFilePaths = cell(5, numel(MD.channels_));
for i = p.ChannelIndex
    [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_Chan' num2str(i) '_selectedGroups'];
    outFilePaths{1,i} = [p.OutputDirectory filesep outFilename '.mat'];

    outFilename = [chanDirName '_Chan' num2str(i) '_featureSpace'];
    outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '.fig'];

    outFilename = [chanDirName '_Chan' num2str(i) '_FluorescenceChannelWithIdsClassified'];
    outFilePaths{3,i} = [p.OutputDirectory filesep outFilename '.fig'];

    outFilename = [chanDirName '_Chan' num2str(i) '_idsClassified'];
    outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '.mat'];

    outFilename = [chanDirName '_Chan' num2str(i) '_presence'];
    outFilePaths{5,i} = [p.OutputDirectory filesep outFilename '.mat'];
end

adhClassProc.setOutFilePaths(outFilePaths);

% numFrames = size(imgMap,3);
% startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
% endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
% movieData to find out pixel size
% pixSize = MD.pixelSize_; % nm/pixel
% tInterval = MD.timeInterval_; % time interval in sec
% scaleBar = 1; %micron
%% Look at distToEdge and ampTotal
% % How to find a threshold for distToEdge?
% % % Plot histogram
% tIntervalMin = tInterval/60; % in min
% periodMin = 1;
% periodFrames = floor(periodMin/tIntervalMin); % early period in frames
% % disp(['Extrapolate by ' num2str(periodMin) ' min...'])
% % tic
% % for k=1:numel(tracksNA)
% %     tracksNA(k) = readIntensityFromTracks(tracksNA(k),imgMap,1,'extraLength',10);
% %     tracksNA(k) = readIntensityFromTracks(tracksNA(k),tMap,2,'extraLength',10);
% % end
% % toc
if p.useSimpleClassification
    disp({'Using the simple classification... Classes will be ';
        '1) non-maturing NAs at the edge; 2) maturing NAs to FAs';
        '5) further growing FAs from existing FA state; and';
        '4) decaying FAs in terms of area, and 6) all the others';
        'group 3 is empty in order to be consistent with 9-class classification'});
    [~,indAll]=distinguishFocalAdhesions(tracksNA,MD,[]);
    idGroup1=indAll{1};
    idGroup2=indAll{2};
    idGroup3=indAll{3};
    idGroup4=indAll{4};
    idGroup5=indAll{5};
    idGroup6=indAll{6};
    idGroup7=indAll{7};
    idGroup8=indAll{8};
    idGroup9=indAll{9};
    save(outFilePaths{1,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3') %This is temporary remedy
    save(outFilePaths{2,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3') %This is temporary remedy
    save(outFilePaths{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3')
%     tableTracksNA = struct2table(tracksNA);
%     save(outFilePaths{5,iChan},'tracksNA','tableTracksNA','-v7.3')
else
    %% Integration of existing classifier(s)
    nTrainingSets = numel(sampleFolders);
    if ~isempty(sampleFolders)
        importSelectedGroups=true;
        for jj=1:nTrainingSets
    %         curImportFilePath = fullfile(PathName,FileName);
            disp(['Loading ' sampleFolders{jj} '...'])
            [idGroupSelectedAll{jj},MDAll{jj},curTracksNA{jj}]=loadSampleFolder(sampleFolders{jj});
        end
    end

    if p.useAutomaticallySelectedData
        % Here is automatically pre-selected samples for G1 (5/18/17)
        % 1. Edge criterion
        % 1-1. Adhesions with zero edge movement should be G6 (not
        % interested) - by looking at edge velocity and edge MSD
        disp('Automatic labeling for G1 ...'); tic;
        edgeVelAll=arrayfun(@(x) regress(x.edgeAdvanceDist(x.startingFrameExtra:x.endingFrameExtra)',(x.startingFrameExtra:x.endingFrameExtra)'), tracksNA);
        edgeStdAll=arrayfun(@(x) nanstd(x.edgeAdvanceDist(x.startingFrameExtra:x.endingFrameExtra)),tracksNA);
        edgeDistAll=arrayfun(@(x) (x.edgeAdvanceDist(x.endingFrameExtra)-x.edgeAdvanceDist(x.startingFrameExtra))', tracksNA);
        advDistAll=arrayfun(@(x) (x.advanceDist(x.endingFrameExtra)-x.advanceDist(x.startingFrameExtra))', tracksNA);

        % 
        % 1. edge should protrude - should see overall positive edge
        % G1 should have decently high edge vel.
        thresEdgeVelG1 = mean(edgeVelAll); %+std(edgeVelAll);
        % In case the mean value is negative, we make it zero as an
        % insurance
        if thresEdgeVelG1<0
            thresEdgeVelG1=0;
        end
        indEdgeVelG1 = edgeVelAll>thresEdgeVelG1;

        % 2. Relative distance from edge should increase significantly.
        distToEdgeVelAll=arrayfun(@(x) regress(x.distToEdge(x.startingFrameExtra:x.endingFrameExtra)',(x.startingFrameExtra:x.endingFrameExtra)'), tracksNA);
        thresRelEdgeVel = max((50/MD.pixelSize_)/(60/MD.timeInterval_), quantile(distToEdgeVelAll,0.05)); %10 nm/min is the absolute vel.
        indRelEdgeVelG1 = distToEdgeVelAll>thresRelEdgeVel;

        % 3. Should start from relatively close to an edge
        try
            distToEdgeFirstAll = arrayfun(@(x) x.distToEdgeNaive(x.startingFrameExtra),tracksNA);
        catch
            distToEdgeFirstAll = arrayfun(@(x) x.distToEdge(x.startingFrameExtra),tracksNA);
        end
        thresStartingDistG1 = p.startingDist*1000/MD.pixelSize_;%quantile(distToEdgeFirstAll,0.25);
        indCloseStartingEdgeG1 = distToEdgeFirstAll < thresStartingDistG1;

        % 5. Clean rising phase
        assemRateAll = arrayfun(@(y) y.assemRate, tracksNA);
    %         thresAssemRateG1G2 = quantile(assemRateAll,0.25);
        indCleanRisingG1G2 = ~isnan(assemRateAll); %assemRateAll>thresAssemRateG1G2;        

        % 6. Clean decaying phase
        disassemRateAll = arrayfun(@(y) y.lateAmpSlope, tracksNA);
        thresDisassemRateG1 = nanmean(disassemRateAll);
        if thresDisassemRateG1<-0.5 || thresDisassemRateG1>0
            thresDisassemRateG1=-0.5;
        end
        indCleanDecayingG1 = disassemRateAll<thresDisassemRateG1;  
    %         disassemRateAll = arrayfun(@(y) y.disassemRate, tracksNA);
    %         thresDisassemRateG1 = quantile(disassemRateAll,0.01);
    %         indCleanDecayingG1 = disassemRateAll>thresDisassemRateG1;  

        % 7. maximum point location compared to the life time (G2 has
        % maximum point at the later phase).
        timeToMaxInten=zeros(numel(tracksNA),1);
        maxIntenAll=zeros(numel(tracksNA),1);
        splineParam=0.01;
        for ii=1:numel(tracksNA)
            curFrameRange = tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra;
            d = tracksNA(ii).ampTotal(curFrameRange);
            tRange = tracksNA(ii).iFrame(curFrameRange);
            warning('off','SPLINES:CHCKXYWP:NaNs')
            d(d==0)=NaN;
            try
                sd_spline= csaps(tRange,d,splineParam);
            catch
                d = tracksNA(ii).amp;
                d(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra) = ...
                    tracksNA(ii).ampTotal(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra);
                sd_spline= csaps(tRange,d,splineParam);
            end
            sd=ppval(sd_spline,tRange);
            %         tRange = [NaN(1,numNan) tRange];
        %     sd = [NaN(1,numNan) sd];
            sd(isnan(d))=NaN;
            %         sd(isnan(d)) = NaN;
            % Find the maximum
            [maxIntenAll(ii),curFrameMaxAmp]=nanmax(sd);
            curFrameMaxAmp = curFrameRange(curFrameMaxAmp);
            timeToMaxInten(ii) = curFrameMaxAmp-tracksNA(ii).startingFrameExtraExtra;
        end
        lifeTimesAll = arrayfun(@(y) y.endingFrameExtra-y.startingFrameExtra, tracksNA);
        relMaxPoints = timeToMaxInten./lifeTimesAll;
        thresRelMax = max(0.8,mean(relMaxPoints));
        indEarlyMaxPointG1 = relMaxPoints<thresRelMax;
        % 8. life time

        % Inspect each (temporary)
    %         figure; hold on
    %         for mm=indexG1'
    %             plot(tracksNA(mm).ampTotal)
    % %             plot(tracksNA(mm).forceMag)
    %         end
    % %         nn=nn+1;
    %         showSingleAdhesionTrackSummary(MD,tracksNA(indexG1(end-3)),imgMap,tMap,indexG1(end-3));

        % G2

        % G2-1. FA segmentation overlapping (for G2)
        faAssocAll = arrayfun(@(x) any(x.state==4 | x.state==3),tracksNA); %any(strcmp(x.state,'FA') | strcmp(x.state,'FC')),tracksNA);
        % Have to think about having to start with NA state, and FC vs. FA

        % G2-2. Area should be increasing overall
        slopeArea=-100*ones(size(faAssocAll));
        meanAreaAll= zeros(size(faAssocAll));
        maxAreaAll= zeros(size(faAssocAll));
        minAreaAll= zeros(size(faAssocAll));
        for k=find(faAssocAll')
            curArea = tracksNA(k).area;
            meanAreaAll(k)=nanmean(curArea);
            maxAreaAll(k)=nanmax(curArea);
            minAreaAll(k)=nanmin(curArea);
            curArea = curArea(~isnan(curArea));
            tRange = tracksNA(k).iFrame(~isnan(curArea));
            slopeArea(k)=regress(curArea',tRange');
        end
        indIncreasingArea = slopeArea>0;

        % G2-3. Edge should not retract
        edgeDistAll=edgeVelAll.*lifeTimesAll;
        indEdgeNotRetracting = edgeDistAll>(10/MD.pixelSize_); %-0.1*std(edgeDistAll);

        indNonDecayingG2 = disassemRateAll>thresDisassemRateG1 | isnan(disassemRateAll);
    %         indNonDecayingG2 = disassemRateAll<thresDisassemRateG1 | isnan(disassemRateAll);
        indLateMaxPointG2 = relMaxPoints>thresRelMax;
    % 
    % %     % Also, G2's max point should be higher than that of G1
    % %     meanIntenG1 = mean(maxIntenAll(indAbsoluteG1));
    % %     stdIntenG1 = std(maxIntenAll(indAbsoluteG1));
    % %     indHighEnoughMaxAmp = maxIntenAll>(meanIntenG1+0.1*stdIntenG1);
    % 
        % G1 should not have too big area
        thresMatureArea=500/MD.pixelSize_*200/MD.pixelSize_;
        smallEnoughArea = meanAreaAll<thresMatureArea; %mean(meanAreaAll)+std(meanAreaAll);
        
%         % Should exclude the ones near the interior boundary. Checking for
%         % edgeVariationNaive
%         if ~isfield(tracksNA,'edgeAdvanceDistNaive') && isfield(tracksNA,'distToEdgeNaive')
%             tracksNANew=getFeaturesFromTracksNA(tracksNA,MD.timeInterval_,1,true);
%             tracksNA = tracksNANew;
%         end
%         if isfield(tracksNA,'edgeAdvanceDistNaive')
%             % edgeVel!=0 will be discarded because it means it is not moving
%             edgeVariation = arrayfun(@(x) min(nanstd(x.closestBdPointNaive(:,1)),nanstd(x.closestBdPoint(:,2))),tracksNA);
%             indEdgeVary = edgeVariation~=0;
%         else
%             indEdgeVary=true(numel(tracksNA),1);
%         end
        
%         % And this G1 should be sliding backward or forward while it is
%         % close to the edge. Thus it's closestBd (in the direction of
%         % sliding) should be the same as closestBDNaive.
%         progressText(0,'Adhesion''s main movement direction', 'Adhesion Classification');
%         directionalityAll=zeros(numel(tracksNA),1);
%         for k=1:numTracks
%             curTrack = tracksNA(k);
%             sF=curTrack.startingFrameExtra; eF=curTrack.endingFrameExtra;
%             try
%                 [~, gof] = fit(curTrack.xCoord(sF:eF)',curTrack.yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
%             catch
%                 % This means fitting a line with the track has failed, highly
%                 % likely due to too short track lifetime or too variable
%                 % lacations or having NaNs
%                 curX = curTrack.xCoord(sF:eF); curY=curTrack.yCoord(sF:eF);
%                 indNaNX = isnan(curX);
%                 t = 1:length(curX);
%                 t_nn = t(~indNaNX);
%                 curX2 = interp1(t_nn,curX(~indNaNX),t,'linear','extrap');
%                 curY2 = interp1(t_nn,curY(~indNaNX),t,'linear','extrap');
% 
%                 [~, gof] = fit(curX2(~isnan(curX2))',curY2(~isnan(curY2))','poly1'); % this is an average linear line fit of the adhesion track
%             end
%             directionalityAll(k)=gof.adjrsquare;
%             tracksNA(k).directionality=gof.adjrsquare;
%             progressText(k/numTracks);
%         end   
%         directionalityAll=arrayfun(@(x) x.directionality,tracksNA);
% %         distBetweenNaiveAndProj = arrayfun(@(x) nanmean(x.distToEdge-x.distToEdgeNaive), tracksNA);
%         indDirectional = directionalityAll>.6;
        
        % We don't want the noise which should belong to G6. 
        meanAmpAll = arrayfun(@(y) nanmean(y.amp), tracksNA);
        maxMeanAmp = max(meanAmpAll);
        meanAmpNorm = meanAmpAll/maxMeanAmp;
        try
            thresMeanAmpNorm = graythresh(meanAmpNorm);
            thresMeanAmp = thresMeanAmpNorm*maxMeanAmp;
        catch
            thresMeanAmp = mean(meanAmpAll)+0.5*std(meanAmpAll);
        end
        lowAmpPopul = meanAmpAll(meanAmpAll<thresMeanAmp);
        thresLowAmpG1 = quantile(lowAmpPopul,.5);
        maxAmpAll = arrayfun(@(y) nanmax(y.ampTotal), tracksNA);
        indHighEnoughMaxInten = maxAmpAll>thresLowAmpG1;
        
        % Summing all those for G1
        indAbsoluteG1 = indEdgeVelG1 & indRelEdgeVelG1 & indCloseStartingEdgeG1 & ...
            indCleanRisingG1G2 & indEarlyMaxPointG1 & smallEnoughArea ...
            & indHighEnoughMaxInten & indCleanDecayingG1; % & indDirectional;
        disp([num2str(sum(indAbsoluteG1)) ' tracks'])
        toc
        
        disp('Automatic labeling for G2...'); tic;
        % At the same time, G2 should start from a decently low amplitude
        % as in G1
        initIntenAll = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),tracksNA);
        if sum(indAbsoluteG1)>1
            initIntenG1 = initIntenAll(indAbsoluteG1);
            indInitIntenG2 = initIntenAll<(mean(initIntenG1)+0.5*std(initIntenG1));
        else
            indInitIntenG2 = initIntenAll<(mean(initIntenAll)-0.5*std(initIntenAll));
        end
        
        % G1 summarizing -> below at Line 505

        % G2 adhesions should have sufficient period where they are in
        % their NA state
        lengthNA = floor(10/MD.timeInterval_); %at least 10 sec
        sufficientInitialNAState = false(numel(tracksNA),1);
        for ii=1:numel(tracksNA)
            curTrack = tracksNA(ii);
            indFirstFA = find(curTrack.state==3 | curTrack.state==4, 1);
            if ~isempty(indFirstFA)
                sufficientInitialNAState(ii) = sum(curTrack.state==2 & ...
                        curTrack.iFrame<indFirstFA)>=lengthNA;
            end
        end
        
        % And G2 should end with FA or FC state in the end
        lastFAFrame = arrayfun(@(x) find(x.state==3 | x.state==4, 1, 'last'),tracksNA,'unif',false);
        lastFAFrame(cellfun(@isempty,lastFAFrame))={0};
        lastFAFrame = cell2mat(lastFAFrame);
        FAfinishing = arrayfun(@(x) x.endingFrameExtra,tracksNA)-lastFAFrame;
        nFrames = nanmax(arrayfun(@(x) x.endingFrameExtraExtra,tracksNA));
        FAfinishing = (nFrames - FAfinishing)/nFrames; %This means that the FA whose state is FA at the end of endingFrameExtra has 1.
        endingWithFAState = FAfinishing>0.95;
        
    %     indAbsoluteG2G1 = indEdgeVelG1 & indRelEdgeVelG1 & indCloseStartingEdgeG1G2 & ...
    %         indCleanRisingG1G2 & indNonDecayingG2 & faAssocAll & indLateMaxPointG2 & indInitIntenG2;
    %     additionalG1 = indAbsoluteG2G1 & ~indHighEnoughMaxAmp;
    %     indAbsoluteG2 = indAbsoluteG2G1 & indHighEnoughMaxAmp;
        thresStartingDistG2 = 4000/MD.pixelSize_;%quantile(distToEdgeFirstAll,0.25);
        indCloseStartingEdgeG2 = distToEdgeFirstAll < thresStartingDistG2;
        
        % Re-considering FA area criteria: I quite don't believe this
        % criterion used in step 7. Basically the G2 adhesions should  ...
        % evolve to a decently large area.  - SH 180422 
        bigEnoughArea = maxAreaAll>thresMatureArea; %mean(meanAreaAll)+std(meanAreaAll);
        
        % One more feature to add to distinguish the maturing adhesion from
        % clumps of nascent adhesion at the cell edge. We should expect
        % that the G2 should end up being a way behind the cell edge at the
        % end of the day, e.g. at least 3 um away
        distToEdgeLastAll = arrayfun(@(x) x.distToEdge(x.endingFrameExtraExtra),tracksNA);
        thresDistAwayFromEdge = 2000/MD.pixelSize_;
        awayfromEdgeLater = distToEdgeLastAll > thresDistAwayFromEdge;
        
        % Yet another feature to distinguish the maturing adhesion from
        % clumps of nascent adhesion is texture, specifically the spatial
        % homogeneity. This needs GLCM matrix calculation of the FA
        % segmentation. We'll need to do this for only when the segmented
        % area is largest to reduce computational burdon. First, let's find
        % the FA segmentation -Sangyoon Nov 19 2018
        % Finding area
        areaAll = arrayfun(@(x) x.area,tracksNA,'unif',false);
        iEmptyArea = cellfun(@isempty,areaAll);
        iFrameMaxArea = zeros(numel(iEmptyArea),1);
        actualMaxArea = zeros(numel(iEmptyArea),1);
        [actualMaxArea(~iEmptyArea),iFrameMaxArea(~iEmptyArea)]=cellfun(@(x) nanmax(x),areaAll(~iEmptyArea));
        framesToInspect = unique(iFrameMaxArea);
        framesToInspect(framesToInspect==0)=[]; %delete frame 0
        % Go over each frame and find the tracks that has the 
        iiformat = ['%.' '3' 'd'];
        iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
        FAPackage=MD.packages_{iFAPack}; iSDCProc=1;
        SDCProc=FAPackage.processes_{iSDCProc};      
        if p.useHomogeneity
            homogeneityAll=zeros(numel(tracksNA),1);
            if ~isfield(tracksNA,'FAtextureHomogeneity')
                tracksNA(end).FAtextureHomogeneity=0;
                [tracksNA(:).FAtextureHomogeneity] = deal(0);
            end
            progressText(0,'Calculating FA texture homogeneity', 'Adhesion Classification');
            for jj=framesToInspect'
                %Find the tracks
                relevantTrackIDs = iFrameMaxArea==jj & actualMaxArea>0;
                if isempty(relevantTrackIDs)
                    continue
                else
                    %Load the image
                    if ~isempty(SDCProc)
                        curImg = SDCProc.loadOutImage(iChan,jj);
                    else
                        curImg = MD.channels_(iChan).loadImage(jj);
                    end
                    % Load the label
                    pAnal=adhAnalProc.funParams_;
                    labelTifPath = [pAnal.OutputDirectory filesep 'labelTifs'];
                    maskAdhesion = imread(strcat(labelTifPath,'/label',num2str(jj,iiformat),'.tif'));
                    labelAdhesion = bwlabel(maskAdhesion>0,4); % strongly assumes each has only one boundary
                    % Per track ID, get the FA image
                    for kk=find(relevantTrackIDs')
                        curFAID = tracksNA(kk).faID(jj);
                        if curFAID==0
                            continue
                        else
                            % Get the individual FA image
                            curFAstruct=regionprops(labelAdhesion==curFAID,'BoundingBox');
                            if ~isempty(curFAstruct)
                                curBoundingBox=curFAstruct.BoundingBox;
                                try
                                    curFAImg = imcrop(curImg,curBoundingBox);
                                catch
                                    curBoundingBox=[tracksNA(kk).xCoord(jj)-5,tracksNA(kk).yCoord(jj)-5,10,10];
                                    curFAImg = imcrop(curImg,curBoundingBox);
                                end                                
                                %Get the gray-level co-occurerence matrix
                                glcm = graycomatrix((curFAImg),'Offset',[1 1],'NumLevels',8,'GrayLimits',[]);
                                stats = graycoprops(glcm);
                                homogeneityAll(kk) = stats.Homogeneity;
                                tracksNA(kk).FAtextureHomogeneity= stats.Homogeneity;
                            end
                        end
                    end
                end
                progressText(jj/framesToInspect(end));
            end
            % See the distribution and determine the threshold figure,histogram(homogeneityAll) 
            thresHomogeneity = max(0.3,mean(homogeneityAll(homogeneityAll>0))+1*std(homogeneityAll(homogeneityAll>0)));
            homogeneousEnoughWhenMatured = homogeneityAll>thresHomogeneity;
        else
            homogeneousEnoughWhenMatured = true(size(indIncreasingArea));
        end
        indAbsoluteG2 = indIncreasingArea & indEdgeNotRetracting & indInitIntenG2 & ...
            sufficientInitialNAState & endingWithFAState & indCloseStartingEdgeG2 & ...
            bigEnoughArea & indCleanRisingG1G2 & awayfromEdgeLater & homogeneousEnoughWhenMatured; % & indNonDecayingG2 & indLateMaxPointG2;
        disp([num2str(sum(indAbsoluteG2)) ' tracks'])
        toc
        disp('Automatic labeling for G3...'); tic;        
        % G3
        thresEdgeVelG3 = max((50/MD.pixelSize_)/(60/MD.timeInterval_),quantile(edgeVelAll,0.75)); %Edge should protrude at least 0.05 um/min
        indEdgeVelG3 = edgeVelAll>thresEdgeVelG3;
        indRelEdgeVelG3 = distToEdgeVelAll<3*thresRelEdgeVel;
        earlyEdgeVelAll = arrayfun(@(x) regress(x.edgeAdvanceDist(x.startingFrameExtra:round((x.startingFrameExtra+x.endingFrameExtra)/2))',(x.startingFrameExtra:round((x.startingFrameExtra+x.endingFrameExtra)/2))'), tracksNA);
        indEarlyEdgeVelG3 = earlyEdgeVelAll>0.3*thresEdgeVelG3;

        adhVelAll=arrayfun(@(x) regress(x.advanceDist(x.startingFrameExtra:x.endingFrameExtra)',(x.startingFrameExtra:x.endingFrameExtra)'), tracksNA);
        thresAdhVelG3 = max(0.001,mean(adhVelAll)+0*std(adhVelAll));
        indForwardVelG3 = adhVelAll>thresAdhVelG3;

        indAbsoluteG3 = indEdgeVelG3 & indRelEdgeVelG3 & indCloseStartingEdgeG1 & indForwardVelG3 & indEarlyEdgeVelG3;
        disp([num2str(sum(indAbsoluteG3)) ' tracks'])
        toc
        % G4 - retracting, strong FAs
        disp('Automatic labeling for G4...'); tic;        
        thresEdgeVelG4 = min((-50/MD.pixelSize_)/(60/MD.timeInterval_), quantile(edgeVelAll(edgeVelAll<0),0.5));%Edge should retract at least 0.1 um/min
        thresEdgeDistG4 = min((-500/MD.pixelSize_), quantile(edgeDistAll(edgeDistAll<0),0.5));%Edge should retract at least 0.1 um/min
        thresAdvDistG4 = min((-500/MD.pixelSize_), quantile(advDistAll(advDistAll<0),0.5));%Edge should retract at least 0.1 um/min
        indEdgeVelG4 = edgeVelAll<thresEdgeVelG4;
        indEdgeDistG4 = edgeDistAll<thresEdgeDistG4;
        indAdvDistG4 = advDistAll<thresAdvDistG4;
        indAdvVelG4 = adhVelAll<thresEdgeVelG4;
        indAbsoluteG4 = indEdgeVelG4 & indAdvVelG4 & indCloseStartingEdgeG1 & faAssocAll & indEdgeDistG4 & indAdvDistG4;
        % Found that there are cases adhesions stay at the place while edge
        % is retracting. This case the adhesion is decaying its intensity.
        % Adding them... Nov 9 2018 Sangyoon
        ampSlopeAll = arrayfun(@(x) x.ampSlope, tracksNA);
        indDecayingAmp = ampSlopeAll<0;
        indAdditionalG4 = indEdgeVelG4 & indDecayingAmp & indCloseStartingEdgeG1 & faAssocAll & indEdgeDistG4 & indAdvDistG4;
        disp([num2str(sum(indAbsoluteG4)) ' tracks'])
        toc
        % G5 - stable at the edge
        disp('Automatic labeling for G5...'); tic;
        % edge doesn't move much
        if isfield(tracksNA,'edgeAdvanceDistNaive') %override if edgeAdvanceDistNaive exists.
            edgeVelAll=arrayfun(@(x) regress(x.edgeAdvanceDistNaive(x.startingFrameExtra:x.endingFrameExtra)',(x.startingFrameExtra:x.endingFrameExtra)'), tracksNA);
        end
        posEdgeVel = edgeVelAll(edgeVelAll>=0);
        negEdgeVel = edgeVelAll(edgeVelAll<=0);
        lowPosEVel = max((5/MD.pixelSize_)/(60/MD.timeInterval_),quantile(posEdgeVel,0.1)); %At most 5 nm/min
        lowNegEVel = max((5/MD.pixelSize_)/(60/MD.timeInterval_),-quantile(negEdgeVel,0.8)); %At most 5 nm/min
        thresEdgeVelG5 = (lowPosEVel+lowNegEVel)/2;
        indEdgeVelG5 = abs(edgeVelAll)<thresEdgeVelG5;

        thresEdgeStdG5 = nanmedian(edgeStdAll); %quantile(edgeStdAll,0.1);
        indEdgeStdG5 = edgeStdAll<thresEdgeStdG5;
        % long life time
        thresLifetimeG5 = quantile(lifeTimesAll,0.75);
        indLifetimeG5 = lifeTimesAll>=thresLifetimeG5;
        % high-enough intensity
        
        % G5 - 2 Should be close to an edge
        try
            distToEdgeMeanAll = arrayfun(@(x) nanmean(x.distToEdgeNaive),tracksNA);
            distToEdgeStdAll = arrayfun(@(x) nanstd(x.distToEdgeNaive),tracksNA);
        catch
            distToEdgeMeanAll = arrayfun(@(x) nanmean(x.distToEdge),tracksNA);
            distToEdgeStdAll = arrayfun(@(x) nanstd(x.distToEdge),tracksNA);
        end
        thresStartingDistG5 = 2000/MD.pixelSize_;%quantile(distToEdgeFirstAll,0.25);
        indCloseEdgeG5 = distToEdgeMeanAll < thresStartingDistG5;
        thresEdgeStd = median(distToEdgeStdAll); %G5 should stay at edge
        indStayAtEdge = distToEdgeStdAll<thresEdgeStd;
        indHighAmpG5 = meanAmpAll>thresMeanAmp;
        indAbsoluteG5 = indStayAtEdge & indLifetimeG5 & indHighAmpG5 & indCloseEdgeG5; %indEdgeVelG5 & indEdgeStdG5
        disp([num2str(sum(indAbsoluteG5)) ' tracks'])
        toc
        % G6 : noise. There are several types of noises, or uninterested
        % tracks
        disp('Automatic labeling for G6...'); tic;
        % G6-1. too short tracks or tracks that are unfinished
        thresShortLifeG6 = min(100/MD.timeInterval_, quantile(lifeTimesAll,0.05));
        indShortLifeG6 = lifeTimesAll<=thresShortLifeG6;
        % G6-2. amplitude too small
        thresLowAmpG6 = mean(lowAmpPopul)-1*std(lowAmpPopul);
        indLowAmpG6 = meanAmpAll<thresLowAmpG6;
        % G6-3. OR, tracks near the image borders (zero edge movement or std
        % I am not sure if assigning two differently positioned labels work
        % for classification - so I'll do this later after classification
        thresStartingDistG6 = 1000/MD.pixelSize_;%quantile(distToEdgeFirstAll,0.25);
        indInsideG6 = distToEdgeFirstAll > thresStartingDistG6; % decided to exclude ones at very edge
        indAbsoluteG6 = (indShortLifeG6 | indLowAmpG6) & indInsideG6;
        disp([num2str(sum(indAbsoluteG6)) ' tracks'])
        toc
        % G7 NAs at stalling edge: big difference from G5 is that it has
        % some early history of edge protrusion & relative weak signal
        % (ampTotal)
        disp('Automatic labeling for G7...'); tic;
        edgeAdvanceDistLastChangeNAs =  arrayfun(@(x) x.edgeAdvanceDistChange2min(x.endingFrameExtra),tracksNA); %this should be negative for group 5 and for group 7
        indLastAdvanceG7 = edgeAdvanceDistLastChangeNAs<max(0.0001, nanmean(edgeAdvanceDistLastChangeNAs));

        thresEdgeVelG7 = max(0.2*thresEdgeVelG3, quantile(edgeVelAll(edgeVelAll>0),0.1));
        indEdgeVelG7 = earlyEdgeVelAll > thresEdgeVelG7;
        indLowAmpG7 = meanAmpAll < thresMeanAmp;

        indAbsoluteG7 = indLastAdvanceG7 & indEdgeVelG7 & indRelEdgeVelG3 & indLowAmpG7;
        disp([num2str(sum(indAbsoluteG7)) ' tracks'])
        toc
        % G8 strong inside FAs
%         distToEdgeMean = arrayfun(@(x) mean(x.distToEdge(x.startingFrameExtra:x.endingFrameExtra)),tracksNA);
%         indInsideG8G9 = distToEdgeMean > thresStartingDistG1;        
        disp('Automatic labeling for G8...'); tic;
        indInsideG8G9 = distToEdgeFirstAll > thresStartingDistG1; % decided to 
        
        indAbsoluteG8 = indHighAmpG5 & indInsideG8G9 & indLifetimeG5;
        disp([num2str(sum(indAbsoluteG8)) ' tracks'])
        toc
        % G9 weak inside NAs
        disp('Automatic labeling for G9...'); tic;
        indMinLifeG9 = lifeTimesAll>2*thresShortLifeG6;
        indAbsoluteG9 = indLowAmpG7 & indInsideG8G9 & indMinLifeG9;
        disp([num2str(sum(indAbsoluteG9)) ' tracks'])
        toc
        %% I decided to get mutually exclusive indices
        indexG1 = find(indAbsoluteG1 & indInitIntenG2 & ~indAbsoluteG2); 
        % Inspect each (temporary)
        indexG2 = find(indAbsoluteG2 & ~indAbsoluteG8 & ~indAbsoluteG9); 
    %     indexG1 = [indexG1; find(additionalG1)];
    %         figure; hold on
    %         for mm=indexG2'
    % %             plot(tracksNA(mm).ampTotal)
    %             plot(tracksNA(mm).forceMag)
    %         end

    %         nn=nn+1; close; showSingleAdhesionTrackSummary(MD,tracksNA(indexG2(nn)),imgMap,tMap,indexG2(nn));
        indexG3 = find(indAbsoluteG3 & ~indAbsoluteG7 & ~indAbsoluteG2 & ~(indAbsoluteG1 & indInitIntenG2)); 
        indexG4 = find(indAbsoluteG4 | indAdditionalG4); 
        indexG5 = find(indAbsoluteG5 & ~(indAbsoluteG4 | indAdditionalG4) & ~indAbsoluteG2 & ~(indAbsoluteG1 & indInitIntenG2)); 
        indexG6 = find(indAbsoluteG6 & ~(indAbsoluteG1 & indInitIntenG2) & ~indAbsoluteG3 & ~indAbsoluteG7...
                        & ~indAbsoluteG9); 
        indexG7 = find(indAbsoluteG7 & ~indAbsoluteG3 & ~indAbsoluteG6 & ~indAbsoluteG9); 
        indexG8 = find(indAbsoluteG8 & ~indAbsoluteG2 & ~(indAbsoluteG1 & indInitIntenG2)...
                    & ~indAbsoluteG3 & ~(indAbsoluteG4 | indAdditionalG4) & ~indAbsoluteG6 ...
                    & ~indAbsoluteG7 & ~indAbsoluteG9 & ~indAbsoluteG5); 
        indexG9 = find(indAbsoluteG9 & ~indAbsoluteG6 & ~indAbsoluteG7 & ~indAbsoluteG8 & ~(indAbsoluteG1 & indInitIntenG2) ...
                        & ~indAbsoluteG2 & ~indAbsoluteG3); 
        

        %% Putting together integrated labels
        % I have to tone down the number of the large label group according to
        % groups with smaller number (especially indexG2)
        indexAll={indexG1, indexG2, indexG3, indexG4, indexG5, indexG6, indexG7, indexG8, indexG9};
        meanSampleNum = round(median(cellfun(@numel,indexAll))); %mean([numel(indexG1) numel(indexG2) numel(indexG3) numel(indexG5)]));
        for ii=1:9
            for jj=ii+1:9
                indCommon = intersect(indexAll{ii},indexAll{jj});
                if numel(indexAll{ii})<numel(indexAll{jj}) %priority to smaller group
                    indexAll{jj} = setdiff(indexAll{jj},indCommon);
                else
                    indexAll{ii} = setdiff(indexAll{ii},indCommon);
                end
            end
        end
        for ii=1:9
            numMax = numel(indexAll{ii});
            if numMax>meanSampleNum
                randNumInt = ceil(numMax*rand(round(meanSampleNum*(1+log(numMax/meanSampleNum))),1));
                indexAll{ii} = indexAll{ii}(randNumInt);
    %             indexAll{ii} = indexAll{ii}(1:meanSampleNum);
            elseif numMax>1
                %Need oversampling up to two thirds of the mean sample
                %number
                while numel(indexAll{ii})<meanSampleNum*(1-exp(-3*numMax/meanSampleNum))
                    indexAll{ii}=[indexAll{ii}; indexAll{ii}];
                end
            end
        end
        idTracksAdditionalAuto = [indexAll{1}' indexAll{2}' indexAll{3}' indexAll{4}' ...
            indexAll{5}' indexAll{6}' indexAll{7}' indexAll{8}' indexAll{9}'];
        iGroupAdditionalAuto = [1*ones(size(indexAll{1}')) 2*ones(size(indexAll{2}')) 3*ones(size(indexAll{3}')) ...
            4*ones(size(indexAll{4}')) 5*ones(size(indexAll{5}')) 6*ones(size(indexAll{6}')) ...
            7*ones(size(indexAll{7}')) 8*ones(size(indexAll{8}')) 9*ones(size(indexAll{9}'))]; % On top of this, we can definitely add other group samples

    %     idTracksAdditionalAuto = [indexG1' indexG2' indexG3' indexG4' indexG5' indexG6' indexG7' indexG8' indexG9'];
    %     iGroupAdditionalAuto = [1*ones(size(indexG1')) 2*ones(size(indexG2')) 3*ones(size(indexG3')) ...
    %         4*ones(size(indexG4')) 5*ones(size(indexG5')) 6*ones(size(indexG6')) ...
    %         7*ones(size(indexG7')) 8*ones(size(indexG8')) 9*ones(size(indexG9'))]; % On top of this, we can definitely add other group samples
    % 
        idTracksAdditionalManual=[]; iGroupAdditionalManual=[]; % I decided to remove this manual selection because it takes forever - 052517
    %         [idTracksAdditionalManual, iGroupAdditionalManual] = showAdhesionTracks(p.OutputDirectory,'all',...
    %             'tracksNA',tracksNA,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD, 'autoIndex', [idTracksAdditionalAuto; iGroupAdditionalAuto]);

        idTracksAdditional = [idTracksAdditionalAuto idTracksAdditionalManual];
        iGroupAdditional = [iGroupAdditionalAuto iGroupAdditionalManual];

        [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
            idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
            sortIDTracks(idTracksAdditional,iGroupAdditional);
        idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        nTrainingSets=nTrainingSets+1;
        curTracksNA{nTrainingSets}=tracksNA;
        idGroupSelectedAll{nTrainingSets}=idGroupSelected;
        MDAll{nTrainingSets}=MD;
    %         [curT,allData,meas] = extractFeatureNA(tracksNA,idGroupSelected,2,MD);
    %         T=[T; curT];
    elseif p.manualLabeling
        fh2=figure;
        fh2.Color=[0 0 0];
        fh2.Position(3)=280;
        fh2.Position(4)=200;
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

        disp('Click tracks that belong to each group ...')
    %     newTracksNA=tracksNA(~idxMatureNAs);
    %     idNAs = find(~idxMatureNAs);
        [idTracksAdditional, iGroupAdditional] = showAdhesionTracks(p.OutputDirectory,'all','tracksNA',tracksNA,'trainedData',T,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD);
    %     [idTracks, iGroup] = showAdhesionTracks(p.OutputDirectory,'all','tracksNA',tracksNA,'iChan',iChan,'iChanSlave',iChanSlave,'movieData',MD);
    %     [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
    %         idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
    %         sortIDTracks(idTracks,iGroup);
        [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
            idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
            sortIDTracks(idTracksAdditional,iGroupAdditional);
        idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
        bigEnoughGroups=cellfun(@length,idGroupSelected);
        bigEnoughGroups=find(bigEnoughGroups>=5);
        idGroupFiltered = idGroupSelected;
        idGroupFiltered(setdiff(1:9,bigEnoughGroups))={[]};
        nTrainingSets=1;
        curTracksNA{nTrainingSets}=tracksNA;
        idGroupSelectedAll{nTrainingSets}=idGroupFiltered;
        MDAll{nTrainingSets}=MD;
    %         [T,allData]=extractFeatureNA(tracksNA,idGroupFiltered,2,MD);
    end
        % figure, plot(asymTracks',MSDall','.')
        % hold on
        % figure,plot(advanceDistNAs(idGroup1)',distToEdgeLastNAs(idGroup1)','ro'), hold on
        % xlabel('asymmetry')
        % ylabel('MSD')
    %     display('Click tracks that belong to group 3 (moving along with protruding edge)')
    %     [idGroup3Selected] = showAdhesionTracks(p.OutputDirectory,'all','tracksNA',newTracksNA);
        % figure, plot(distToEdgeLastNAs(idGroup3)',distToEdgeChangeNAs(idGroup3)','go')
        % plot(advanceDistNAs(idGroup3)',distToEdgeLastNAs(idGroup3)','go')
        % figure, histogram(advanceDistNAs(idGroup1),1:30); hold on; histogram(advanceDistNAs(idGroup3),1:30)
        % figure, plot(asymTracks(idGroup3)',MSDall(idGroup3)','go')
        % create linear discriminator
    save(outFilePaths{1,i},'idGroup1Selected',...
        'idGroup2Selected','idGroup3Selected','idGroup4Selected','idGroup5Selected','idGroup6Selected',...
        'idGroup7Selected','idGroup8Selected','idGroup9Selected');

    %% feature extraction
    %% visualize feature space and p-dist (similarity) matrix
    % features = meas;

    %% classifier training or import
    useDefinedClassifier='n'; %input('Do you want to use already defined classifier? [y/(n)]: ','s');
    if isempty(useDefinedClassifier)
        useDefinedClassifier='n';
    end

    if strcmp(useDefinedClassifier,'n')
    %     if importSelectedGroups && strcmp(reuseSelectedGroups, 'u')
    %         bigEnoughGroups=cellfun(@length,idGroupSelected);
    %         bigEnoughGroups=find(bigEnoughGroups>=5);
    %         idGroupFiltered = idGroupSelected;
    %         idGroupFiltered(setdiff(1:9,bigEnoughGroups))={[]};
    %         [T,allData]=extractFeatureNA(tracksNA,idGroupFiltered,2,MD);
    %     end
        T = table();kk=2;
        for jj=1:nTrainingSets
            T=[T; extractFeatureNA(curTracksNA{jj},idGroupSelectedAll{jj},kk,MDAll{jj})];
        %         T=extractFeatureNA(curTracksNA{jj},idGroupSelected{jj},ii,curMD);
        end

        disp(['Training SVM classifier (' num2str(size(T,1)) ' labeled data)...'])
        tic
        [trainedClassifierSVM, validationAccuracySVM, CSVM, orderSVM] = trainClassifierNA(T);
        toc
    %     [trainedClassifierKNN, validationAccuracyKNN, CKNN, orderKNN] = trainClassifierKNN(T);
        % I will use SVM no matter what, because it will be compatible with
        % multiple different data
    %     if validationAccuracySVM>validationAccuracyKNN
        trainedClassifier=trainedClassifierSVM;
        validationAccuracy=validationAccuracySVM;
        C=CSVM;
        order=orderSVM;
        classifierInfo = fopen([p.OutputDirectory filesep 'data' filesep 'trainedClassifier is from SVM.txt'],'w');
        fprintf(classifierInfo, 'This is from quadratic SVM. \n');
        fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracy) '. \n']);
    %     fprintf(classifierInfo, ['The other (Weighted KNN) was ' num2str(validationAccuracyKNN) '. \n']);
    %     else
    %         trainedClassifier=trainedClassifierKNN;
    %         validationAccuracy=validationAccuracyKNN;
    %         C=CKNN;
    %         order=orderKNN;
    %         classifierInfo = fopen([p.OutputDirectory filesep 'data' filesep 'trainedClassifier is from KNN.txt'],'w');
    %         fprintf(classifierInfo, 'This is from Weighted KNN. \n');
    %         fprintf(classifierInfo, ['Validation accuracy is ' num2str(validationAccuracy) '. \n']);
    %         fprintf(classifierInfo, ['The other (SVM) was ' num2str(validationAccuracySVM) '. \n']);
    %     end
    %     [trainedClassifier, validationAccuracy, C, order,validationPredictions, validationScores] = trainClassifierNA(T);
    %     classificationLearner
        disp(['Validation accuracy is ' num2str(validationAccuracy) '.'])
        fclose(classifierInfo);
        save([p.OutputDirectory filesep 'data' filesep 'trainedClassifier.mat'],'trainedClassifier')

        % normalize confusion matrix
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
    %     while validationAccuracy<0.6 && strcmp(reuseSelectedGroups,'u')
    %         disp('Validation accuracy was low. Group reduction is needed.')
    %         interestedGroups = input('Which groups are in this specific movie? Use bracket form...');
    %         % Re-formatting T with interestedGroups...
    %         idGroupFiltered = idGroupSelected;
    %         idGroupFiltered(setdiff(1:9,interestedGroups))={[]};
    %         [T,allData]=extractFeatureNA(tracksNA,idGroupFiltered,2,MD);
    %         [trainedClassifier, validationAccuracy, C, order] = trainClassifierNA(T);
    %         disp(['New validation accuracy is ' num2str(validationAccuracy) '.'])
    %         % normalize confusion matrix
    %         for ii=1:size(C,1)
    %             C(ii,:) = C(ii,:)/sum(C(ii,:));
    %         end
    %         response = T.Group;
    %         % Get the unique resonses
    %         totalGroups = unique(response);
    %         
    %         figure; confAxis=axes; imagesc(C); title(['Confusion Matrix: Validation accuracy: ' num2str(validationAccuracy)])
    %         set(confAxis,'xticklabel',totalGroups')
    %         set(confAxis,'yticklabel',totalGroups')
    %         c = colorbar;
    %         c.Label.String = 'normalized prediction';
    %         print('-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'confusionMatrix.eps']);
    %         savefig([p.OutputDirectory filesep 'figs' filesep 'confusionMatrix.fig'])
    %         print('-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'eps' filesep 'confusionMatrix.tif'])
    %     end

        T = sortrows(T,size(T,2));
        features =table2array(T(:,1:end-1));
        species = table2array(T(:,end));
        nGroups = length(totalGroups);
        % normalize features
        for ii = 1 : size(features,2)
            features(:,ii) = (features(:,ii) - min(features(:,ii)))./(max(features(:,ii)) - min(features(:,ii)));
        end
        featFig=figure; imagesc(features');hold on
        c = colorbar;
        c.Label.String = 'feature value';
        print(featFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'featureSpace.eps']);
        savefig(featFig, outFilePaths{2,p.ChannelIndex})
        print(featFig,'-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'featureSpace.tif'])

        D = pdist(features);
        D1 =  squareform(D);
        similDataFig=figure; imagesc(D1);
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
        print(similDataFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'similarityAmongTrainedData.eps']);
    %     savefig([p.OutputDirectory filesep 'figs' filesep 'similarityAmongTrainedData.fig'])
        print(similDataFig,'-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'similarityAmongTrainedData.tif'])
%         close

        Dfeats = pdist(features');
        Dfeats1 =  squareform(Dfeats);
        dfeatFig=figure; imagesc(Dfeats1); title('similarityAmongFeatures')
        shading flat
        c = colorbar;
        c.Label.String = 'p-dist';

        print(dfeatFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'similarityAmongFeatures.eps']);
    %     savefig([p.OutputDirectory filesep 'figs' filesep 'similarityAmongFeatures.fig'])
        print(dfeatFig, '-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'similarityAmongFeatures.tif'])
%         close
        [~,allData] = extractFeatureNA(tracksNA,[],2,MD);

    %     disp('The order is :')
    %     disp(order)
    elseif strcmp(useDefinedClassifier,'y')
        [fileDefinedClassifier, pathDefinedClassifier]=uigetfile('*.mat', 'Select the mat file for trained classifier');
        trainedClassifier=load(fullfile(pathDefinedClassifier,fileDefinedClassifier));
        trainedClassifier=trainedClassifier.trainedClassifier;
        disp(['Classifier at ' fullfile(pathDefinedClassifier,fileDefinedClassifier) ' is used for cross-validation.'])
        % Check the Class names between classifier and T, and adjust T
        idxTinClassifier=ismember(T.Group,trainedClassifier.ClassNames);
        T= T(idxTinClassifier,:);
        [validationAccuracy, C, order] = validateClassifier(trainedClassifier,T);
        disp(['Cross-validation accuracy is ' num2str(validationAccuracy) '.'])

        % normalize confusion matrix
        for ii=1:size(C,1)
            C(ii,:) = C(ii,:)/sum(C(ii,:));
        end
        figure; imagesc(C); title('Confusion Matrix')
        c = colorbar;
        c.Label.String ='normalized prediction';
        print('-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'confusionMatrix_otherClassifier.eps']);
        savefig([p.OutputDirectory filesep 'figs' filesep 'confusionMatrix_otherClassifier.fig'])
        print('-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'eps' filesep 'confusionMatrix_otherClassifier.tif'])
        disp('The order is :')
    %     disp(order)
        [~,allData] = extractFeatureNA(tracksNA,[],2,MD);
    end

    allDataClass = predict(trainedClassifier,allData);

    % figure, imshow(imgMap(:,:,end),[])
    % hold on
    % colors = distinguishable_colors(9,'k');
    % htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1));
    % htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup2),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup2));
    % htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup3),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup3));
    % htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(4,:)),tracksNA(idGroup4),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(4,:)),tracksNA(idGroup4));
    % htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(5,:)),tracksNA(idGroup5),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(5,:)),tracksNA(idGroup5));
    % htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(9,:)),tracksNA(idGroup6),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(9,:)),tracksNA(idGroup6));
    % htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(7,:)),tracksNA(idGroup7),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(7,:)),tracksNA(idGroup7));
    % htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(8,:)),tracksNA(idGroup8),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(8,:)),tracksNA(idGroup8));
    % htrackG9=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(6,:)),tracksNA(idGroup9),'UniformOutput',false);
    % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(6,:)),tracksNA(idGroup9));
    % legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1} htrackG9{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
    %     'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
    % legend('boxoff')
    idGroup1 = strcmp(allDataClass,'Group1');
    idGroup2 = strcmp(allDataClass,'Group2');
    idGroup3 = strcmp(allDataClass,'Group3');
    idGroup4 = strcmp(allDataClass,'Group4');
    idGroup5 = strcmp(allDataClass,'Group5');
    idGroup6 = strcmp(allDataClass,'Group6');
    idGroup7 = strcmp(allDataClass,'Group7');
    idGroup8 = strcmp(allDataClass,'Group8');
    idGroup9 = strcmp(allDataClass,'Group9');
    
    %% minority group redemption
    % There is a case where idgroup with small number of labels is totally
    % neglected during classifier building. This case we'll put the labels
    % back to the idGroup
    if sum(idGroup2)<numel(indexG2)
        idGroup2(indexG2)=true;
    end
    if sum(idGroup1)<numel(indexG1)
        idGroup1(indexG1)=true;
    end
    if sum(idGroup3)<numel(indexG3)
        idGroup3(indexG3)=true;
    end
    
    %% Adhesion boundary homogenization
    % This is one last step where we make potentially different classes
    % associated with the same chunk of adhesion boundary the same.
    % We base on the labels obtained in Step 6, adhesion segmentation,
    % which is stored in tracksNA.faID
    faIDsCell = arrayfun(@(x) nanmax(x.faID),tracksNA,'unif',false);
    % Make the ones with empty faID to ID zero
    indEmpty = cellfun(@isempty,faIDsCell);
    faIDsCell(indEmpty)={0};
    indFaIDs = ~cellfun(@(x) isnan(x) || (x==0),faIDsCell);
    faIDsAll = cellfun(@(x) uint16(x), faIDsCell(indFaIDs));
    [faIDsUnique]=unique(faIDsAll); %,iAlltoUniq,iUniqToAll
    indexFaIDs=find(indFaIDs);
    idGroupLabel= 1*idGroup1 + 2*idGroup2 + 3*idGroup3 + ...
                4*idGroup4 + 5*idGroup5 + 6*idGroup6 + ...
                7*idGroup7 + 8*idGroup8 + 9*idGroup9;    
    edgeVariationAll = allData(:,end-2);
    edgeAdvanceSpeedAll = allData(:,2);
    for k=faIDsUnique'
        % Find the tracks associated with k, one of the faIDs
        curTrackIDs = indexFaIDs(faIDsAll==k); % This can be multiple
        % find all classes associated with these tracks
        curClasses = idGroupLabel(curTrackIDs);
        % Apply the rule for the representative Class
        curDistToEdges = arrayfun(@(x) x.distToEdge(x.startingFrameExtra), tracksNA(curTrackIDs));
        curEdgeVariation = edgeVariationAll(curTrackIDs);
        curEdgeAdvanceSpeed = edgeAdvanceSpeedAll(curTrackIDs);
        curSufficientInitialNAState = sufficientInitialNAState(curTrackIDs);
        curRepClass = mode(curClasses);
        if mean(curEdgeAdvanceSpeed)<=0 && min(curDistToEdges)>1.5*thresStartingDistG1 % If the edge is protruding, we'll leave this chuck to be regarded as NAs
            if (ismember(8,curClasses) && (mean(curDistToEdges)>=2*thresStartingDistG1) || sum(curEdgeVariation==0)/length(curEdgeVariation)>0.5)
                curRepClass = 8;
            elseif ismember(2,curClasses) && min(curDistToEdges)<1.5*thresStartingDistG1 && mode(curSufficientInitialNAState)
                curRepClass = 2;
%             elseif ismember(4,curClasses)
%                 curRepClass = 4; % G4 has the highest priority
            elseif ismember(5,curClasses) && min(curDistToEdges)<1.5*thresStartingDistG1 && sum(curEdgeVariation==0)<1
                curRepClass = 5;
            elseif curRepClass==9 && nanmean(arrayfun(@(x) nanmean(x.area),tracksNA(curTrackIDs)))>(mean(meanAreaAll(indexG1))+0.5*std(meanAreaAll(indexG1)))
                if min(curDistToEdges)<1*thresStartingDistG1 && mode(curSufficientInitialNAState)
                    curRepClass = 2;
                else
                    curRepClass = 8;
                end
            end
            % Putting this back to idGroupLabel
            idGroupLabel(curTrackIDs) = curRepClass;
        end
    end
    % Putting back to idGroups
    idGroup1=idGroupLabel==1; idGroup2=idGroupLabel==2; idGroup3=idGroupLabel==3; 
    idGroup4=idGroupLabel==4; idGroup5=idGroupLabel==5; idGroup6=idGroupLabel==6; 
    idGroup7=idGroupLabel==7; idGroup8=idGroupLabel==8; idGroup9=idGroupLabel==9; 
    %% Sample drawing
    allDataClass(idGroup1)={'Group1'};
    allDataClass(idGroup2)={'Group2'};
    allDataClass(idGroup3)={'Group3'};
    allDataClass(idGroup4)={'Group4'};
    allDataClass(idGroup5)={'Group5'};
    allDataClass(idGroup6)={'Group6'};
    allDataClass(idGroup7)={'Group7'};
    allDataClass(idGroup8)={'Group8'};
    allDataClass(idGroup9)={'Group9'};
    iFrameInterest=round(0.8*MD.nFrames_);
    mapFig = figure; imshow(imgMap(:,:,iFrameInterest),[]), hold on
    [~, htrackCircles] = drawClassifiedTracks(allDataClass,tracksNA,iFrameInterest,[],10);
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
    %% Saving
    if strcmp(useDefinedClassifier,'n')
        print(mapFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified.eps']);
        savefig(mapFig,outFilePaths{3,iChan})
        print(mapFig,'-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'FluorescenceChannelWithIdsClassified.tif'])
        save(outFilePaths{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3')
%         tableTracksNA = struct2table(tracksNA);
%         save(outFilePaths{5,iChan},'tracksNA','tableTracksNA','-v7.3') Doesn't need to store tracksNA in Step 8 because it's not changed 
    else
        print(mapFig,'-depsc2', '-r300', [p.OutputDirectory filesep 'eps' filesep 'FluorescenceChannelWithIdsClassified_otherClassifier.eps']);
        savefig(mapFig,[p.OutputDirectory filesep 'figs' filesep 'FluorescenceChannelWithIdsClassified_otherClassifier.fig'])
        print(mapFig,'-dtiff', '-loose', '-r300', [p.OutputDirectory filesep 'tif' filesep 'FluorescenceChannelWithIdsClassified_otherClassifier.tif'])
        save(outFilePaths{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3')
%         tableTracksNA = struct2table(tracksNA);
%         save(outFilePaths{5,iChan},'tracksNA','tableTracksNA','-v7.3') Doesn't need to store tracksNA in Step 8 because it's not changed 
    end
end
% saving presence
s = struct2table(tracksNA);
presenceAll = s.presence;
save(outFilePaths{5,iChan},'presenceAll')

disp('Classification Process Done!')
end
% function  [validationAccuracy,C,order] = validateClassifier(trainedClassifier,datasetTable) %This function is now separate.
% % Extract predictors and response
% predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
% predictors = datasetTable(:,predictorNames);
% predictors = table2array(varfun(@double, predictors));
% response = datasetTable.Group;
% % Perform cross-validation
% % figure; imagesc(predictors);
% % predictedLabels = predict(trainedClassifier,predictors);
% % [predictedLabels,NegLoss,PBScore] = trainedClassifier.predict(predictors);
% predictedLabels = trainedClassifier.predict(predictors);
% results = nan(1,numel(predictedLabels));
% for i = 1 : numel(predictedLabels)
%     results(i) = strcmp(predictedLabels{i},response{i});
% end
% validationAccuracy=sum(results)/length(results);
% % confusion matrix
% [C,order] = confusionmat(response,predictedLabels);
% end
%% classificationLearner
% SVMmodel = fitcecoc(meas,species);
% isLoss = resubLoss(SVMmodel)
% CVmodel = crossval(SVMmodel)
% oosLoss = kfoldLoss(CVmodel)
% % %% Here let's filter out group 4 and 5
% % % Look at each class - first is expected to have very negative
% % % advance
% % idGroup4 = [idGroup4Selected];
% % nonGroup135=[idGroup1Selected idGroup3Selected idGroup5Selected];
% % meas = [advanceDistNAs(idGroup4) edgeAdvanceDistNAs(idGroup4);
% %     advanceDistNAs(nonGroup135) edgeAdvanceDistNAs(nonGroup135)];
% % % meas = [advanceDistNAs(idGroup4Selected) edgeAdvanceDistNAs(idGroup4Selected);
% % %     advanceDistNAs(nonGroup24Selected) edgeAdvanceDistNAs(nonGroup24Selected)];
% % nG4=length(idGroup4);
% % nGn4 = length(nonGroup135);
% % nTotalG = nG4+nGn4;
% % species = cell(nTotalG,1);
% % for ii=1:nTotalG
% %     if ii<=nG4
% %         species{ii} = 'Group4';
% %     else 
% %         species{ii} = 'NonGroup4';
% %     end
% % %     elseif ii<=nG4+nGn4
% % %         species{ii} = 'NonGroup4';
% % %     else 
% % %         species{ii} = 'Group3';
% % %     end
% % end
% %% some manual group reassignment
% % species{1}= 'NonGroup4';
% % species{37}= 'Group4';
% % species{36}= 'Group4';
% % species{32}= 'Group4';
% 
% % %% plotting classes
% % mea1 = meas(:,1);
% % mea2 = meas(:,2);
% % figure
% % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % h1(1).LineWidth = 2;
% % h1(2).LineWidth = 2;
% % % h1(3).LineWidth = 2;
% % legend('Group 4','Non-Group 4','Location','best')
% % hold on
% % cls = fitcdiscr(meas,species);
% % % % Plot the classification boundaries.
% % K = cls.Coeffs(1,2).Const; % First retrieve the coefficients for the linear
% % L = cls.Coeffs(1,2).Linear;% boundary between the second and third classes
% %                            % (versicolor and virginica).
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2 = ezplot(f,[min(mea1(1:nG4+nGn4)) max(mea1(1:nG4+nGn4)) min(mea2(1:nG4+nGn4)) max(mea2(1:nG4+nGn4))]); 
% % h2.Color = 'r';
% % h2.LineWidth = 2;
% % xlabel('advanceDistNAs')
% % ylabel('edgeAdvanceDistNAs')
% 
% % showAdhesionTracks(pathForColocalization,70,'tracksNA',tracksNA);
% %% filter out group 4 from the population
% % % allData = [advanceDistNAs edgeAdvanceDistNAs];
% % % allDataClass = predict(cls,allData);
% % % idxGroup4 = strcmp(allDataClass,'Group4');
% % % idxGroup4(idGroup4Selected) = true;
% % % idGroup4 = find(idxGroup4);
% % % %% Working on group 5 - long life time, steadyly low distToEdgeChange expected...
% % % % group 6 - just noise...
% % % % Just get rid of them
% % % % disp('to do here')
% % % idGroup56 = [idGroup5Selected idGroup6Selected];
% % % 
% % % %% Group 2 - maturing
% % % idxMatureNAs = idGroup2Selected;
% % % % maxIntensityNAs(idGroup2Selected)
% % % % maxIntensityNAs(idxMatureNAs)
% % % % figure, histogram(maxIntensityNAs,50)
% % % % % Then, let's filter out group 2 using state information
% % % % idxMatureNAs =  arrayfun(@(x) x.maturing==1,tracksNA);
% % % % % showAdhesionTracks(pathForColocalization,idGroup6Selected,'tracksNA',tracksNA);
% % % % 
% % % % % Here, it is possible that wrong FA segmentation affected wrong group 2
% % % % % definition. We need to rescue tracks whose edge protrude and adhesion
% % % % % also moves along to the edge. In other words, distToEdge relatively
% % % % % constant AND edge protruded.
% % % % distToEdgeLastMatureNAs = arrayfun(@(x) x.distToEdge(x.endingFrame),tracksNA(idxMatureNAs));
% % % % % edgeAdvanceDistMatureNAs = arrayfun(@(x) x.edgeAdvanceDist,tracksNA(idxMatureNAs));
% % % % adhAdvanceDistMatureNAs = arrayfun(@(x) x.advanceDist,tracksNA(idxMatureNAs));
% % % % % figure, histogram(distToEdgeLastMatureNAs,50)
% % % % % figure, histogram(edgeAdvanceDistMatureNAs,50)
% % % % % figure, histogram(adhAdvanceDistMatureNAs,50)
% % % % % figure, plot(distToEdgeChangeMatureNAs,edgeAdvanceDistMatureNAs,'k.')
% % % % % figure, plot(adhAdvanceDistMatureNAs,edgeAdvanceDistMatureNAs,'k.')
% % % % idxRetractingAdhs = adhAdvanceDistMatureNAs<0 & distToEdgeLastMatureNAs>10;
% % % % % matureTracks = tracksNA(idxMatureNAs);
% % % % % reallyMatureTracks = matureTracks(idxRetractingAdhs);
% % % % p=0;
% % % % for ii=find(idxMatureNAs)'
% % % %    p= p+1;
% % % %    idxMatureNAs(ii)=idxRetractingAdhs(p);
% % % % end
% % % 
% % % %% other processes
% % % % meas = [advanceDistNAs(idGroup1Selected);
% % % %     advanceDistNAs(idGroup3Selected)];
% % % % meas = [MSDrate(idGroup1Selected) ;
% % % %     MSDrate(idGroup3Selected) ];
% % % % meas = [relMovWRTEdge(idGroup1Selected) MSDrate(idGroup1Selected) ;
% % % %     relMovWRTEdge(idGroup3Selected) MSDrate(idGroup3Selected) ];
% % % % meas = [MSDrate(idGroup1Selected) relMovWRTEdge(idGroup1Selected);
% % % %     MSDrate(idGroup3Selected) relMovWRTEdge(idGroup3Selected)];
% % % % meas = [MSDrate(idGroup1Selected) advanceDistNAs(idGroup1Selected);
% % % %     MSDrate(idGroup3Selected) advanceDistNAs(idGroup3Selected)];
% % % % meas = [advanceDistNAs(idGroup1Selected), asymTracks(idGroup1Selected);
% % % %     advanceDistNAs(idGroup3Selected),asymTracks(idGroup3Selected)];
% % % % meas = [advanceDistNAs(idGroup1), distToEdgeLastNAs(idGroup1);
% % % %     advanceDistNAs(idGroup3),distToEdgeLastNAs(idGroup3)];
% % % % meas = [asymTracks(idGroup1), MSDall(idGroup1);
% % % %     asymTracks(idGroup3),MSDall(idGroup3)];
% % % 
% % % % nG1=length(idGroup1Selected);
% % % % nG3 = length(idGroup3Selected);
% % % % nTotalG = nG1+nG3;
% % % % species = cell(nTotalG,1);
% % % % for ii=1:nTotalG
% % % %     if ii<=nG1
% % % %         species{ii} = 'group1';
% % % %     else
% % % %         species{ii} = 'group3';
% % % %     end
% % % % end
% % % 
% % % % In case for more than two features:
% % % % linclass = fitcdiscr(meas,species);
% % % % mea1 = meas(:,1);
% % % % mea2 = meas(:,2);
% % % % figure
% % % % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % % % h1(1).LineWidth = 2;
% % % % h1(2).LineWidth = 2;
% % % % legend('Group 1','Group 3','Location','best')
% % % % hold on
% % % % % Plot the classification boundaries.
% % % % K = linclass.Coeffs(1,2).Const; % First retrieve the coefficients for the linear
% % % % L = linclass.Coeffs(1,2).Linear;% boundary between the second and third classes
% % % %                            % (versicolor and virginica).
% % % % % Plot the curve K + [x,y]*L  = 0.
% % % % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % % % h2 = ezplot(f,[0 max(mea1) 0 max(mea2)]);
% % % % h2.Color = 'r';
% % % % h2.LineWidth = 2;
% % % % ylabel('MSD rate')
% % % % xlabel('Adhesion movement relative to edge rate')
% % % 
% % % %% Now separate group1 and 3 using MSD rate
% % % % In case of one feature, I 'll use mean of the interfacing tales of the
% % % % two distributions:
% % % sortedMSDrateG1 = sort(MSDrate(idGroup1Selected));
% % % sortedMSDrateG3 = sort(MSDrate(idGroup3Selected));
% % % thresMSDrate = mean([sortedMSDrateG1(end-1:end); sortedMSDrateG3(1:2)]);
% % % % additional filtering of maturing adhesions and adhesions in the
% % % % retracting edge ....
% % % arrayfun(@(x) plot(x.xCoord,x.yCoord,'b'),tracksNA(idGroup1Selected));
% % % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'bo'),tracksNA(idGroup1Selected));
% % % arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idGroup3Selected));
% % % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idGroup3Selected));
% % % showAdhesionTracks(pathForColocalization,idGroup2Selected,'tracksNA',tracksNA);
% % % % predict now
% % % allData = [MSDrate ];
% % % % allData = [relMovWRTEdge MSDrate ];
% % % % allData = [MSDrate relMovWRTEdge];
% % % % allData = [advanceDistNAs MSDall];
% % % % allData = [advanceDistNAs, distToEdgeLastNAs];
% % % % allDataClass = predict(linclass,allData);
% % % % idxStaticNAs = strcmp(allDataClass,'group1');
% % % % idxAdvancingNAs= strcmp(allDataClass,'group3');
% % % idxStaticNAs =allData<thresMSDrate;
% % % idxAdvancingNAs= allData>=thresMSDrate;
% % % idStaticNAs = find(idxStaticNAs);
% % % idStaticNAs = setdiff(idStaticNAs, [idxMatureNAs'; idGroup4; idGroup56']);
% % % idAdvancingNAs = find(idxAdvancingNAs);
% % % idAdvancingNAs = setdiff(idAdvancingNAs, [idxMatureNAs'; idGroup4; idGroup56']);
% % % % For advancing NAs, there is sub-group where adhesions are in stalling
% % % % edge. These should be separated. These should have low edgeAdvanceDist,
% % % % distToEdge, and low distToEdgeChange...
% % % % Look at MSDrate first
% % % figure, histogram(MSDrate(idGroup3Selected),0:1:15)
% % % % For those who have low MSDrate, check edgeAdvanceDist compared to those with high
% % % % MSDrate
% % % MSDrateG3 = MSDrate(idGroup3Selected);
% % % edgeAdvanceG3 = edgeAdvanceDistNAs(idGroup3Selected);
% % % distToEdgeChangeG3 =distToEdgeChangeNAs(idGroup3Selected);
% % % opts = statset('maxIter', 200);
% % % for n = 1:3
% % %     objMSDG3{n} = gmdistribution.fit(MSDrateG3, n, 'Options', opts);
% % % end
% % % [~,idx] = min(cellfun(@(i) i.BIC, objMSDG3));
% % % objMSDG3 = objMSDG3{idx};
% % % [mu,idx] = sort(objMSDG3.mu);
% % % svec = sqrt(squeeze(objMSDG3.Sigma(:,:,idx)));
% % % thresMSDrate= mu(idx(1))+2*svec(idx(1));
% % % idxLowMSDrateG3 = MSDrateG3<thresMSDrate;
% % % MSDrateG3(idxLowMSDrateG3)
% % % MSDrateG3(~idxLowMSDrateG3)
% % % edgeAdvanceG3(idxLowMSDrateG3)
% % % edgeAdvanceG3(~idxLowMSDrateG3)
% % % figure, scatter(MSDrate,edgeAdvanceDistNAs)
% % % hold on
% % % scatter(MSDrate(idGroup4),edgeAdvanceDistNAs(idGroup4),'ro')
% % % scatter(MSDrate(idGroup5Selected),edgeAdvanceDistNAs(idGroup5Selected),'co')
% % % scatter(MSDrate(idGroup6Selected),edgeAdvanceDistNAs(idGroup6Selected),'go')
% % % scatter(MSDrate(idGroup3Selected),edgeAdvanceDistNAs(idGroup3Selected),'ko')
% % % scatter(MSDrate(idGroup3Selected(idxLowMSDrateG3)),edgeAdvanceDistNAs(idGroup3Selected(idxLowMSDrateG3)),'go')
% % % scatter(MSDrate(idGroup1Selected),edgeAdvanceDistNAs(idGroup1Selected),'mo')
% % % 
% % % idStallingNAs;
% % confirm
% % figure, histogram(allData(idxStaticNAs),0:0.1:15); hold on; histogram(allData(idxAdvancingNAs),0:0.1:15)
% % figure, histogram(MSDrate(idGroup1Selected),0:0.1:15); hold on, histogram(MSDrate(idGroup3Selected),0:0.1:15)
% % figure
% % gscatter(allData(:,1),allData(:,2),allDataClass,'krb','ov^',[],'off');
% % legend('Group 1','Group 3','Location','best')
% % hold on
% % h1 = gscatter(mea1,mea2,species,'krb','ov^',[],'off');
% % h1(1).LineWidth = 2;
% % h1(2).LineWidth = 2;
% % % Plot the classification boundaries.
% % h2 = ezplot(f,[min(allData(:,1)) max(allData(:,1)) min(allData(:,2)) max(allData(:,2))]);
% % h2.Color = 'r';
% % h2.LineWidth = 2;
% % xlabel('MSD rate')
% % ylabel('Adhesion movement relative to edge rate')
% 
% figure, imshow(imgMap(:,:,end),[])
% hold on
% htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'g'),tracksNA(idStaticNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'go'),tracksNA(idStaticNAs));
% htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'r'),tracksNA(idAdvancingNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'ro'),tracksNA(idAdvancingNAs));
% htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'c'),tracksNA(idxMatureNAs),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'co'),tracksNA(idxMatureNAs));
% htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'y'),tracksNA(idxGroup4),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'yo'),tracksNA(idxGroup4));
% htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'m'),tracksNA(idGroup5Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'mo'),tracksNA(idGroup5Selected));
% htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup6Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup6Selected));
% htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup7Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup7Selected));
% htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'w'),tracksNA(idGroup8Selected),'UniformOutput',false);
% arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'wo'),tracksNA(idGroup8Selected));
% % arrayfun(@(x) plot(x.xCoord,x.yCoord,'y'),reallyMatureTracks);
% % arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'yo'),reallyMatureTracks);
% % idGroup1 = idNAs(idxStaticNAs);
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1}],{'G1:turn-over','G2:maturing','G3:moving along edge',...
%     'G4:retracting','G5:stable','G6:noise or very transient','G7:adhesions at stalling edge','G8:uncartegorized'},'TextColor','w','Location','best')
% legend('boxoff')
% 
% idGroup1 = idxValidNAs(idStaticNAs);
% idGroup2 = idxValidNAs(idxMatureNAs);
% % idGroup3 = idNAs(idxAdvancingNAs);
% idGroup3 = idxValidNAs(idAdvancingNAs);
% idGroup4 = idxValidNAs(idxGroup4);
% idGroup5 = idxValidNAs(idGroup5Selected);
% idGroup6 = idxValidNAs(idGroup6Selected);
% idGroup7 = idxValidNAs(idStallingNAs);
% disp('classified.')
% save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','tracksNA')
% % title('{\bf Linear Classification with Fisher Training Data}')
% 
% % return ids of 
% 
% %% Filtering is necessary
% % Filtering FAs out: tracks whose amplitude goes down during the first
% % minute or start 
% % initialAmpSlope = arrayfun(@(x) x.earlyAmpSlope,tracksNA);
% % thresInitSlope = multithresh(initialAmpSlope,4);
% % minSlope = thresInitSlope(find(thresInitSlope<0,1,'last'));
% % idxSlope = initialAmpSlope>minSlope;
% % %backup old tracks
% % 
% % initialAmpTotal=arrayfun(@(x) x.ampTotal(x.startingFrame),tracksNA);
% % thresAmp = multithresh(initialAmpTotal,2);
% % idxInitAmp = initialAmpTotal<thresAmp(1);
% % 
% % meanDist=arrayfun(@(x) x.advanceDist,tracksNA);
% % thresDist = multithresh(meanDist,10);
% % idxDist = meanDist>thresDist(1);
% % 
% % isNA = idxSlope & idxInitAmp & idxDist;
% % for ii=find(isNA)'
% %     tracksNA(ii).idxNA = true;
% % end
% % for ii=find(~isNA)'
% %     tracksNA(ii).idxNA = false;
% % end
% % newTracksNA = tracksNA(isNA);
% % % figure,hist(initialAmpTotal,200)
% % 
% % % % starting with median, find a edge disconnected with two consequtive
% % % % zeros.
% % % medC = median(pstruct.c);
% % % idxAfterMedC=find(edges>medC);
% % % qq=idxAfterMedC(1);
% % % while N(qq)>0 || N(qq+1)>0
% % %     qq=qq+1;
% % %     if qq>=length(edges)-1
% % %         break
% % %     end
% % % end
% % % idx = pstruct.c<edges(qq);
% % 
% % toc
% %% Make it to feature matrix and save
% % disp('Making it to feature matrix and saving...')
% % numFeatures = 5; % I have five features
% % naFeatMat = zeros(numFeatures,numel(newTracksNA)); 
% % for k=1:numel(newTracksNA)
% %     naFeatMat(:,k) = [newTracksNA(k).earlyAmpSlope; 
% %                                 newTracksNA(k).lateAmpSlope;
% %                                 newTracksNA(k).distToEdgeSlope;
% %                                 newTracksNA(k).edgeAdvanceDist;
% %                                 newTracksNA(k).advanceDist];
% % end
% % %throwing away observation containing NaN
% % idxNan=any(isnan(naFeatMat),1);
% % naFeatMat(:,idxNan)=[];
% % origNAFeatMat = naFeatMat;
% %% scaling of the X to -1 and 1 for slopes and -2 to 2 for distance
% % scaleMat = zeros(size(naFeatMat,1),1);
% % for jj=1:size(naFeatMat,1)
% %     if jj==1 
% %         % Find 2 percentile and 98 percentile
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Find which one is larger in magnitude
% %         biggerMaxSlop=max(abs(minRow), abs(maxRow));
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     elseif jj==2
% %         % use the same scale as column 1
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     else
% %         % Find 2 percentile and 98 percentile
% %         minRow = quantile(naFeatMat(jj,:),.02);
% %         maxRow = quantile(naFeatMat(jj,:),.98);
% %         % Find which one is larger in magnitude
% %         biggerMaxSlop = 0.5*max(abs(minRow), abs(maxRow));
% %         % Scale with biggerMaxSlop
% %         naFeatMat(jj,:) = naFeatMat(jj,:)/biggerMaxSlop;
% %         relMinRow = minRow/biggerMaxSlop;
% %         relMaxRow = maxRow/biggerMaxSlop;
% %         %Make data below 2 percentile relMinRow
% %         naFeatMat(jj,naFeatMat(jj,:)<relMinRow) = relMinRow;
% %         %Make data above 98 percentile one
% %         naFeatMat(jj,naFeatMat(jj,:)>relMaxRow) = relMaxRow;
% %         scaleMat(jj)=biggerMaxSlop;
% %     end
% % end
% % % save
% % save([pathForColocalization filesep 'data' filesep 'naFeatMat.mat'],'origNAFeatMat','naFeatMat','newTracksNA','isNA','tracksNA','scaleMat')
% disp('Done!')
% % Let's put a hard threshold first to see this criteria work. But later I
% % would like to do some type of machine learning 
% % group 1: distance increases AND ampTotal fall down with a negative slope
% 

% case study:
% % group 1: distance increases AND ampTotal fall down to noise level
% ii=166; fRange=tracksNA(ii).startingFrame:tracksNA(ii).endingFrame; figure, plot(fRange,tracksNA(ii).distToEdge(fRange))
% fRange=tracksNA(ii).startingFrameExtra:tracksNA(ii).endingFrameExtra; figure, plot(fRange,tracksNA(ii).ampTotal(fRange),'r');
% % group 2: distance increases AND ampTotal keep increasing
% ii=449; fRange=tracksNA(ii).startingFrame:tracksNA(ii).endingFrame; figure, plot(fRange,tracksNA(ii).distToEdge(fRange))
% fRange=tracksNA(ii).startingFrameExtra:tracksNA(ii).endingFrameExtra; figure, plot(fRange,tracksNA(ii).ampTotal(fRange),'r');
% % group 3: distance does not increase AND ampTotal does not increase


% see histogram
% figure,hist(arrayfun(@(x) x.ampSlope,tracksNA),200)

