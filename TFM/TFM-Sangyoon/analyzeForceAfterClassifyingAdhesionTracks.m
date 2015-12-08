function []=analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath, uptoColo)
%analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath)
%extracts an average behavior of traction forces from tracks after
%supervised classification.
% e.g.:
% pathForTheMovieDataFile = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/2015-04-29/cell2/Int4sec';
% outputPath = 'analysis2';
% analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath)
if nargin<3
    uptoColo = false;
end
%% Run colocalizationAdhesionsWithTFM if it's not run
% showAllTracks = false;
% plotEachTrack = false;
% tmaxEach = [];
tmax=2500;
[pathStr,name,ext]= fileparts(pathForTheMovieDataFile);
if strcmp([name ext],'movieData.mat')
    pathForTheMovieDataFile = pathStr;
end
outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
outputFile=strcat(dataPath,filesep,'tracksNA.mat');
orgTrackFile=strcat(dataPath,filesep,'tracksNAOrg.mat');
load([pathForTheMovieDataFile filesep 'movieData.mat'])
load([pathForTheMovieDataFile filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat'],'forceField')
forceSpacing = forceField(1).pos(2,2)-forceField(1).pos(1,2);
band = 2*forceSpacing;
% See if you can use existing tracks
if ~exist('tracksNA','var')
    if exist(outputFile,'file') && ~uptoColo
        disp('tracksNA file is found. Using it ... If you do not want to reuse this, please backup the file and rename it to other name than tracksNA.mat.')
        tic
        tracksNAFile = load(outputFile,'tracksNA');
        tracksNA = tracksNAFile.tracksNA;
        toc
        disp('Done loading tracks.')
    else
        % run analyzeAdhesionMaturation for obtaining tracks from paxillin channel
        tracksNA=colocalizationAdhesionsWithTFM(pathForTheMovieDataFile,band,tmax,false,false,[],...
            'onlyEdge',false,'outputPath',outputPath,'saveAnalysis',true,'matchWithFA',true,'minLifetime',10);
    end
end
%% classification after colocalization
if ~uptoColo
    %% filter out matured adhesion
    %back up tracksNA
    save(orgTrackFile,'tracksNA','-v7.3');
    pathForColocalization = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
    % load([pathForColocalization filesep 'data'  filesep 'tracksNA.mat'])
%     [tracksNA, ~]=separateMatureAdhesionTracks(tracksNA);
%     tracksNA = filterOutNonEmergingNA(tracksNA);
    if ~isfield(tracksNA,'SDC_applied')
        disp('Applying stage drift correction ...')
        iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
        if ~isempty(iSDCProc)
            SDCProc=MD.processes_{iSDCProc};
            if ~SDCProc.checkChannelOutput(1)
                error(['The channel must have been corrected ! ' ...
                    'Please apply stage drift correction to all needed channels before '...
                    'running displacement field calclation tracking!'])
            end
            if length(SDCProc.funParams_.ChannelIndex)>1
                iChan = 2;
            elseif length(SDCProc.funParams_.ChannelIndex) == 1
                iChan = SDCProc.funParams_.ChannelIndex;
            else
                error('No channel associated with SDC process!')
            end
            if iChan==2
                iBeadChan=1;
            else
                iBeadChan = SDCProc.funParams_.ChannelIndex(1);
            end
            s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
            T = s.T;
        end
        tracksNA = applyDriftToTracks(tracksNA, T); % need some other function....formatNATracks(tracksNAorg,detectedNAs,nFrames,T); 
    else
        disp('Stage drift correction was already applied to tracksNA.')
    end
    [~,idEmerging] = filterOutNonEmergingNA(tracksNA);
    [~,idFiltered] = filterTracksNAwithMask(pathForTheMovieDataFile,outputPath,tracksNA,band);
    tracksNA = recalculateLifeTimeTracks(tracksNA);
    %% Classify with ampTotal and distToEdge
    [idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]= classifyNascentAdhesionTracks(pathForColocalization,'tracksNA',tracksNA);
    %% filter with idFiltered
    idGroup1filtered = idGroup1 & idFiltered & idEmerging;
    idGroup2filtered = idGroup2 & idFiltered;
    idGroup3filtered = idGroup3 & idFiltered & idEmerging;
    idGroup4filtered = idGroup4 & idFiltered;
    idGroup5filtered = idGroup5 & idFiltered;
    idGroup6filtered = idGroup6 & idFiltered;
    idGroup7filtered = idGroup7 & idFiltered & idEmerging;
    idGroup8filtered = idGroup8 & idFiltered;
    idGroup9filtered = idGroup9 & idFiltered & idEmerging;
    %% turn all zeros to NaNs in ampTotal
    for ii=1:numel(tracksNA)
        tracksNA(ii).ampTotal(tracksNA(ii).ampTotal==0)=NaN;
        tracksNA(ii).forceMag(tracksNA(ii).forceMag==0)=NaN;
    end
    %% Filtering for group1
    % lateAmpSlopeG1 = arrayfun(@(x) x.lateAmpSlope,tracksNA(idGroup1));
    % idxLateAmpSlopeG1 = lateAmpSlopeG1<0;
    % We filter out tracks whose ampTotal is too high, bigger than mean value
    % of ampTotal maxima
    meanAmpMaximum = mean(arrayfun(@(x) nanmax(x.ampTotal),tracksNA(idGroup1filtered)));
    ampSlopeG1 = arrayfun(@(x) x.ampSlope,tracksNA(idGroup1filtered));
    lateAmpTotalG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNA(idGroup1filtered));
    initForceG1 = arrayfun(@(x) x.forceMag(x.startingFrame),tracksNA(idGroup1filtered));
    idxIncreasingAmpG1 = ampSlopeG1>0;
    idxLowInitForceG1= initForceG1<500;
    idxLateAmpLow = lateAmpTotalG1<meanAmpMaximum;
    idGroup1nl = find(idGroup1filtered); %converting to non-logical index
    idGroup1f = idGroup1nl(idxLateAmpLow & idxIncreasingAmpG1 & idxLowInitForceG1);

    %% drawing group1
    epsPath=[pathForColocalization filesep 'eps'];
    fileStore = [epsPath filesep 'ampForcePlotG1.eps'];
    plotIntensityForce(tracksNA(idGroup1f),fileStore,false,false)
    %% Filtering for group2
    ampSlopeG2 = arrayfun(@(x) x.ampSlope,tracksNA(idGroup2filtered));
    initForceG2 = arrayfun(@(x) x.forceMag(x.startingFrame),tracksNA(idGroup2filtered));
    lifeTimeG2 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup2filtered));
    ampEndingG2 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNA(idGroup2filtered));
    ampStartingG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),tracksNA(idGroup2filtered));
    idxIncreasingAmpG2 = ampSlopeG2>=0 & ampEndingG2>ampStartingG2;
    idxLowInitForceG2= initForceG2<500;
    idxLongLifeTimeG2=lifeTimeG2>40;
    idGroup2nl = find(idGroup2filtered); %converting to non-logical index
    idGroup2f = idGroup2nl(idxIncreasingAmpG2 & idxLowInitForceG2 & idxLongLifeTimeG2);
    %% group 2
    fileStoreG2 = [epsPath filesep 'ampForcePlotG2.eps'];
    plotIntensityForce(tracksNA(idGroup2f),fileStoreG2,false,true)

    %% group 3
    % lifeTimeG3 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup3));
    % [~,longestID]=max(lifeTimeG3);
    % % longestID = idGroup3(longestID);
    % idxGroup3=find(idGroup3);
    % ii=idxGroup3(longestID);
    % figure
    % plot(1:(tracksNA(ii).lifeTime),tracksNA(ii).amp(tracksNA(ii).presence))
    %% group 3 plotting

    fileStoreG3 = [epsPath filesep 'ampForcePlotG3.eps'];
    plotIntensityForce(tracksNA(idGroup3filtered),fileStoreG3,false,false)

    %% group4 plotting
    fileStoreG4 = [epsPath filesep 'ampForcePlotG4.eps'];
    plotIntensityForce(tracksNA(idGroup4filtered),fileStoreG4,false,false)
    %% group5 plotting
    fileStoreG5 = [epsPath filesep 'ampForcePlotG5.eps'];
    plotIntensityForce(tracksNA(idGroup5filtered),fileStoreG5,false,false)
    %% group6 plotting
    fileStoreG6 = [epsPath filesep 'ampForcePlotG6.eps'];
    plotIntensityForce(tracksNA(idGroup6filtered),fileStoreG6,false,false)
    %% group7 plotting
    fileStoreG7 = [epsPath filesep 'ampForcePlotG7.eps'];
    plotIntensityForce(tracksNA(idGroup7filtered),fileStoreG7,false,false)
    %% group8 plotting
    fileStoreG8 = [epsPath filesep 'ampForcePlotG8.eps'];
    plotIntensityForce(tracksNA(idGroup8filtered),fileStoreG8,false,false)
    %% group9 plotting
    fileStoreG9 = [epsPath filesep 'ampForcePlotG9.eps'];
    plotIntensityForce(tracksNA(idGroup9filtered),fileStoreG9,false,false)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% For group 1,3,7, look at the cross-correlation score, and collect those with high scores and see if they represent adhesion/edge protrusion  behavior
    % group 1
%     ccG1 = arrayfun(@(x) x.CCscore,tracksNA(idGroup1f),'UniformOutput',false);
    ccG1 = arrayfun(@(x) x.CCscoreMax,tracksNA(idGroup1f));
    % if we pick up adhesions wth high correlation, there should be an
    % increasing phase in the force too...
    idxHighCCG1=ccG1>0.5;
    % delete very short life time
    LTG1 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup1f));
    idxLongLTG1=LTG1>20;
    fileStoreG1_hiCC = [epsPath filesep 'ampForcePlotG1_hiCC.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCC,false,true) %)),fileStoreG1_hiCC,false,true)% 
    ccLagG1Hi = arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)));
    mean(ccLagG1Hi)
    %% those with intermediate cc
    idxIntmedCCG1=ccG1>0.1 & ccG1<=0.5;
    % delete very short life time
    fileStoreG1_middleCC = [epsPath filesep 'ampForcePlotG1_middleCC.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_middleCC,false,true)
    %% those with low and negative cc
    idxLowCCG1=ccG1<0.1;
    % delete very short life time
    fileStoreG1_lowCC = [epsPath filesep 'ampForcePlotG1_lowCC.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_lowCC,false,true)
    %% Which group has tracks with high cc-score?
    ccAll = arrayfun(@(x) x.CCscoreMax,tracksNA);
    hiCCNAs = ccAll>0.5;
    idGroups = {idGroup1filtered,idGroup2filtered,idGroup3filtered,idGroup4filtered,idGroup5filtered,idGroup6filtered,idGroup7filtered,idGroup8filtered,idGroup9filtered};
    numhiCC = zeros(9,1);
    ratiohiCC = zeros(9,1);
    numTracksPerGroups = zeros(9,1);
    for ii=1:9
        numTracksPerGroups(ii) = sum(idGroups{ii});
        numhiCC(ii) = sum(idGroups{ii} & hiCCNAs);
        ratiohiCC(ii) = numhiCC(ii) /  numTracksPerGroups(ii);
    end
    [~,indMaxCC]=max(ratiohiCC);
    disp(['A group containing the most ratio of high CC score: group ' num2str(indMaxCC) ' with ratio of ' num2str(ratiohiCC(indMaxCC)) '.'])
    iGroups = (1:9)';

    ratioT = table(iGroups,numhiCC,numTracksPerGroups,ratiohiCC);
    disp(ratioT)
    writetable(ratioT,[pathForColocalization filesep 'data' filesep 'ratioTracksHiCCscore.csv'])
    %% Event alignment for Group1
    % since there is variance in peaks of vinculin, the rise-and-fall behavior
    % is not picked up in the average. So I'll align these event with the
    % maximum of the intensity timeseries
    fileStoreG1_hiCCshifted = [epsPath filesep 'ampForcePlotG1_hiCC_shifted.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCCshifted,true,true)

    fileStoreG1_middleCCshifted = [epsPath filesep 'ampForcePlotG1_middleCC_shifted.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_middleCCshifted,true,true)

    fileStoreG1_lowCCshifted = [epsPath filesep 'ampForcePlotG1_lowCC_shifted.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_lowCCshifted,true,true)

    fileStoreG1_Allshifted = [epsPath filesep 'ampForcePlotG1_All_shifted.eps'];
    plotIntensityForce(tracksNA(idGroup1f(idxLongLTG1)),fileStoreG1_Allshifted,true,true)
    %% Event alignment for Group3
%     % since there is variance in peaks of vinculin, the rise-and-fall behavior
%     % is not picked up in the average. So I'll align these event with the
%     % maximum of the intensity timeseries
%     ccG3 = arrayfun(@(x) x.CCscore,tracksNA(idGroup3));
%     LTG3 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup3));
%     meanForceG3 = arrayfun(@(x) mean(x.forceMag),tracksNA(idGroup3));
%     idxLongLTG3=LTG3>20;
%     idxHighCCG3=ccG3>0.5;
%     idxIntmedCCG3=ccG3>0.1 & ccG3<=0.5;
%     idxLowCCG3=ccG3<0.1;
%     lowForceG3 = meanForceG3<20;
%     idxGroup3 = find(idGroup3);
%     fileStoreG3_hiCCshifted = [epsPath filesep 'ampForcePlotG3_hiCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idxGroup3(idxHighCCG3 & idxLongLTG3)),fileStoreG3_hiCCshifted,1)
% 
%     fileStoreG3_hiCCshiftedLowForce = [epsPath filesep 'ampForcePlotG3_hiCC_shifted_lowForce.eps'];
%     plotIntensityForce(tracksNA(idxGroup3(idxHighCCG3 & idxLongLTG3 & lowForceG3)),fileStoreG3_hiCCshiftedLowForce,1)
% 
%     fileStoreG3_middleCCshifted = [epsPath filesep 'ampForcePlotG3_middleCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idxGroup3(idxIntmedCCG3 & idxLongLTG3)),fileStoreG3_middleCCshifted,1)
% 
%     fileStoreG3_lowCCshifted = [epsPath filesep 'ampForcePlotG3_lowCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idxGroup3(idxLowCCG3 & idxLongLTG3)),fileStoreG3_lowCCshifted,1)
% 
%     fileStoreG3_Allshifted = [epsPath filesep 'ampForcePlotG3_All_shifted.eps'];
%     plotIntensityForce(tracksNA(idxGroup3(idxLongLTG3)),fileStoreG3_Allshifted,1)


    %% Where are they?
    % Now we have to make a function that shows in the last vinculin image, the
    % tracks of interest with labels 
    imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
    imgMap = imgMap.paxImgStack;
    colors = distinguishable_colors(3,'k');

    figure, imshow(imgMap(:,:,end),[])
    hold on
    htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),'UniformOutput',false);
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)));
    htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),'UniformOutput',false);
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)));
    htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),'UniformOutput',false);
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)));
    try
        legend([htrackG1{1} htrackG2{1} htrackG3{1}],{'G1 with high CC','G1 with intermediate CC','G1 with low CC'},'TextColor','w','Location','best')
    catch
        try
            legend([htrackG1{1} htrackG2{1}],{'G1 with high CC','G1 with intermediate CC'},'TextColor','w','Location','best')
        catch
            legend([htrackG2{1}],{'G1 with intermediate CC'},'TextColor','w','Location','best')
        end
    end
    legend('boxoff')
    print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithTracksInGroup1.eps']);
    %% See if there is any time lag
%     arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLongLTG1)));
%     idGroup1 = idGroup1 | idGroup7;
    ccLagG1 = arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1filtered));
    % see the profiles of tracks that have too much lags
%     idxTooMuchLag=find(ccLagG1>20);
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchLag(1))),[],false,true)
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchLag(2))),[],false,true)
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchLag(3))),[],false,true)
%     idxTooMuchNegLag=find(ccLagG1<-20);
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchNegLag(1))),[],false,true)
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchNegLag(2))),[],false,true)
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchNegLag(3))),[],false,true)
%     plotIntensityForce(tracksNA(idGroup1f(idxTooMuchNegLag(4))),[],false,true)
    meanTimelagG1=mean(ccLagG1);
    stdTimelagG1=std(ccLagG1);
    avgTimeLag = mean(ccLagG1);
    avgTimeLagHiCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1))));
    stdTimeLagHiCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1))));
    avgTimeLagMidCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1))));
    stdTimeLagMidCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1))));
    avgTimeLagLoCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1))));
    stdTimeLagLoCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1))));
    %% initial time estimation
    % noise level is estimated in the previous time point, and the time
    % when the signal is above the background for the first time
%     curTrack = tracksNA(287);
%     plotIntensityForce(curTrack,[],false,true)
    % Check the noise level in the previous frame
    if ~exist('imgMap','var')
        disp('Loading image stacks ...')
        tic
        imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
        imgMap = imgMap.paxImgStack;
        toc
    end
    if ~exist('tMap','var')
        disp('Loading traction stacks ...')
        tic
        tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
        tMap = tMap.tMap;
        toc
    end
    %% read intensity again
    tracksNA = readIntensityFromTracks(tracksNA,imgMap,1,'reTrack',false,'extraLength',30);
    tracksNA = readIntensityFromTracks(tracksNA,tMap,2,'reTrack',true);
    
%     % look at 95 % confidence interval of before-detection
%     meanIntenBeforeLevel = mean(curTrack.ampTotal(curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra));
%     stdIntenBeforeLevel = std(curTrack.ampTotal(curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra));
    %% Force confiddence
    % First we need to filter out force-inconfident tracks
    if ~isfield(tracksNA,'forceUncertainty')
        tracksNA = readForceUncertaintyFromTracks(pathForTheMovieDataFile,'tracksNA',tracksNA);
    end
    % Collect fCfd
    fCfdG1 = (arrayfun(@(x) mean(x.forceUncertainty(x.startingFrame:x.endingFrame)),tracksNA(idGroup1f)));
    figure, histogram(fCfdG1)
    % Maybe I can throw away tracks with below 0.3 of force confidence.
    % Then jump directly to figure out whether each track increases its
    % force during its presence compared to the earlier time points
    % clear tMap imgMap
    %% Initial force-increase time quantification!
    tInterval = MD.timeInterval_;
    curIndices = find(idGroup1filtered | idGroup3filtered | idGroup7filtered | idGroup9filtered)';
%     curIndices = idGroup1f';
    for ii=curIndices
        curTrack = tracksNA(ii);
        sFEE = curTrack.startingFrameExtraExtra;
%         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
        % See how many frames you have before the startingFrameExtra
        numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
        numPreFrames = max(1,floor(1/5*numFramesBefore));
        numPreSigStart = min(20,numFramesBefore);
        sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
        bkgMaxInt = max(curTrack.ampTotal(sFEE:sF5before));
        firstIncreseTimeInt = find(curTrack.ampTotal>bkgMaxInt,1);
%         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstIncreseTimeInt)
            bkgMaxForce = max(curTrack.forceMag(sFEE:sF5before));
            firstIncreseTimeForce = find(curTrack.forceMag>bkgMaxForce,1);
            if isempty(firstIncreseTimeForce) || firstIncreseTimeForce>curTrack.endingFrameExtra
                tracksNA(ii).forceTransmitting = false;
                tracksNA(ii).firstIncreseTimeInt = [];
                tracksNA(ii).firstIncreseTimeForce = [];
                tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            else
                tracksNA(ii).forceTransmitting = true;
                tracksNA(ii).firstIncreseTimeInt = firstIncreseTimeInt*tInterval; % in sec
                tracksNA(ii).firstIncreseTimeForce = firstIncreseTimeForce*tInterval;
                tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeInt*tInterval - firstIncreseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            end
        else
            tracksNA(ii).forceTransmitting = false;
            tracksNA(ii).firstIncreseTimeInt = [];
            tracksNA(ii).firstIncreseTimeForce = [];
            tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
            tracksNA(ii).bkgMaxInt = bkgMaxInt;
            tracksNA(ii).bkgMaxForce = [];
        end
    end
    % statistics about firstIncreseTimeIntAgainstForce
    firstIncreseTimeIntAgainstForceAllIdx = arrayfun(@(x) ~isempty(x.forceTransmitting) && x.forceTransmitting,tracksNA(curIndices));
    firstIncreseTimeIntAgainstForceAll = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdx)));
    numForceTransmitting = length(firstIncreseTimeIntAgainstForceAll);
    numNonTransmitting = sum(~firstIncreseTimeIntAgainstForceAllIdx);
    figure, histogram(firstIncreseTimeIntAgainstForceAll)
    median(firstIncreseTimeIntAgainstForceAll)
    %% See ccLag
    figure, histogram(ccLagG1)
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'timeInitLagsG1.mat'],'ccLagG1','firstIncreseTimeIntAgainstForceAll','numForceTransmitting','numNonTransmitting')
    %% Now it's time to quantify peak intensity vs. peak force in Group 1
    %   perform spline filter for ampTotal
    splineParam=0.01;
    for ii=curIndices
        d = tracksNA(ii).ampTotal;
        nTime = length(d);
        tRange = 1:nTime;
        numNan = find(isnan(d),1,'last');
        tRange(isnan(d)) = [];
        d(isnan(d)) = [];
        sd_spline= csaps(tRange,d,splineParam);
        sd=ppval(sd_spline,tRange);
        d = [NaN(1,numNan) d];
%         tRange = [NaN(1,numNan) tRange];
        sd = [NaN(1,numNan) sd];
%         sd(isnan(d)) = NaN;
        % Find the maximum
        [~,curFrameMaxAmp]=max(sd);
        curFrameLocMaxes=locmax1d(sd,7);
        
        if ismember(curFrameMaxAmp,curFrameLocMaxes) && curFrameMaxAmp>tRange(1) && curFrameMaxAmp<tRange(end)
            tracksNA(ii).intenPeakness = true;
            tracksNA(ii).intenPeakFrame = curFrameMaxAmp;
            curForce = tracksNA(ii).forceMag;
            curForce(isnan(curForce)) = [];
            sCurForce_spline= csaps(tRange,curForce,splineParam);
            sCurForce=ppval(sCurForce_spline,tRange);
            sCurForce = [NaN(1,numNan) sCurForce];
            curForceLocMaxes=locmax1d(sCurForce,7);
            % delete the first and last frame from loc maxes
            curForceLocMaxes = setdiff(curForceLocMaxes, [tRange(1) tRange(end)]);
            if ~isempty(curForceLocMaxes)
                [forceMacMax,indMaxForceAmongLMs] = max(sCurForce(curForceLocMaxes));
                tracksNA(ii).forcePeakness = true;
                tracksNA(ii).forcePeakFrame = curForceLocMaxes(indMaxForceAmongLMs);
                tracksNA(ii).forcePeakMag = forceMacMax; 
                tracksNA(ii).forcePeakMagRel = forceMacMax - min(sCurForce); % this is a relative difference
                tracksNA(ii).peakIntTimeIntAgainstForcePeak = -tracksNA(ii).forcePeakFrame*tInterval + tracksNA(ii).intenPeakFrame*tInterval; % - means intensity comes first; + means force peak comes first
            else
                tracksNA(ii).forcePeakness = false;
                tracksNA(ii).forcePeakFrame = NaN;
                tracksNA(ii).forcePeakMag = NaN; 
                tracksNA(ii).forcePeakMagRel = NaN; % this is a relative difference
                tracksNA(ii).peakIntTimeIntAgainstForcePeak = NaN; %
            end
        else
            tracksNA(ii).intenPeakness = false;
            tracksNA(ii).intenPeakFrame = [];
            tracksNA(ii).forcePeakness = false;
            tracksNA(ii).forcePeakFrame = NaN;
            tracksNA(ii).forcePeakMag = NaN; 
            tracksNA(ii).forcePeakMagRel = NaN; % this is a relative difference
            tracksNA(ii).peakIntTimeIntAgainstForcePeak = NaN; %
        end
    end
    %% statistics about peakTimeIntAgainstForceAll
    peakTimeIntAgainstForceAllIdx = arrayfun(@(x) ~isempty(x.intenPeakness) && x.intenPeakness && ~isempty(x.forcePeakness) && x.forcePeakness,tracksNA(curIndices));
    peakTimeIntAgainstForceAll = arrayfun(@(x) x.peakIntTimeIntAgainstForcePeak, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    peakForceAll = arrayfun(@(x) x.forcePeakMag, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    peakForceRelAll = arrayfun(@(x) x.forcePeakMagRel, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    figure, histogram(peakTimeIntAgainstForceAll)
    figure, histogram(peakForceAll)
    figure, histogram(peakForceRelAll)
    median(peakTimeIntAgainstForceAll)
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'timePeaksG1.mat'],'peakTimeIntAgainstForceAll','peakForceAll','peakForceRelAll')
    
    %% Correlation between firstIncreseTimeIntAgainstForce vs. peakTimeIntAgainstForceAll
    bothTiTpIdx = arrayfun(@(x) ~isempty(x.intenPeakness) && x.intenPeakness && ...
        ~isempty(x.forcePeakness) && x.forcePeakness && ~isempty(x.forceTransmitting) && x.forceTransmitting,tracksNA(curIndices));
    firstIncreseTimeIntAgainstForce_both = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndices(bothTiTpIdx)));
    peakTimeIntAgainstForceAll_both = arrayfun(@(x) x.peakIntTimeIntAgainstForcePeak, tracksNA(curIndices(bothTiTpIdx)));
    figure, plot(firstIncreseTimeIntAgainstForce_both,peakTimeIntAgainstForceAll_both,'ko'),xlabel('t_{init}lag'),ylabel('t_{peak}lag')
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'corrTinitTpeak.mat'],'firstIncreseTimeIntAgainstForce_both','peakTimeIntAgainstForceAll_both')
    %% Lifetime vs. force-trasmissivity
    % Hypothesis: short-lived NAs will have low chance of force-increase.
    % recalculate lifeTime
    lifeTimesG1 = arrayfun(@(x) x.endingFrameExtra-x.startingFrameExtra, tracksNA(curIndices));
    forceTransmissivityG1 = arrayfun(@(x) x.forceTransmitting, tracksNA(curIndices));
    % Accumulate force-transmissivities to lifetime bins
    [NLT,edgesLT]=histcounts(lifeTimesG1);
    % figure,plot((edgesLT(1:end-1)+edgesLT(2:end))/2,NLT,'o-')
    lifeTimeBins = [edgesLT(1) edgesLT(2) edgesLT(3) edgesLT(4) edgesLT(5) edgesLT(6) edgesLT(8) edgesLT(11) edgesLT(end)];
    lifeTimeBinPoints = (lifeTimeBins(1:end-1)+lifeTimeBins(2:end))/2;
    nLTBins = length(lifeTimeBinPoints);
    cumFTrans = zeros(size(lifeTimeBinPoints));
    cumFTransAll = zeros(size(lifeTimeBinPoints));
    nG1 = length(curIndices);
    for ii=1:nG1
        curLT = lifeTimesG1(ii);
        for jj=1:nLTBins
            if curLT>=lifeTimeBins(jj) && curLT<lifeTimeBins(jj+1)
                cumFTrans(jj)=cumFTrans(jj)+forceTransmissivityG1(ii);
                cumFTransAll(jj)=cumFTransAll(jj)+1;
                break
            end
        end
    end
    % make them to normalized probability
    probFTrans = cumFTrans./cumFTransAll;
    figure, plot(lifeTimeBinPoints,probFTrans,'s-')
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'corrLifetimeForceTransmissivity.mat'],'probFTrans','cumFTrans','cumFTransAll','lifeTimeBinPoints',...
        'NLT','lifeTimesG1','lifeTimeBins','nLTBins','nG1')
    %% figure about extended time series - this is used for showing quantification representation
%     firstIncreseTimeIntAgainstForceAllIdxIdx = find(firstIncreseTimeIntAgainstForceAllIdx);
%     curID = idGroup1f(firstIncreseTimeIntAgainstForceAllIdxIdx(157));
%     idLongG1 = curIndices(lifeTimesG1>20);
%     curID=idLongG1(39);
% force-non-transmitting
%     curIndices = find(idGroup1)';
%     curID = 1395;
%% For non-transmitting
    nonTransmittingIdx = ~firstIncreseTimeIntAgainstForceAllIdx;
    nonTransmittingIdx = find(nonTransmittingIdx);
    curID = nonTransmittingIdx(25);
    curTrack = tracksNA(curIndices(curID));%4210);
    plotIntensityForceSingle(curTrack,curIndices(curID));
    %% for transmitting
    firstIncreseTimeIntAgainstForceAllIdxIDs = find(firstIncreseTimeIntAgainstForceAllIdx);
    %% Actual plotting
    curID = firstIncreseTimeIntAgainstForceAllIdxIDs(56);
    curTrack = tracksNA(curIndices(curID));%4210);
%     curTrack = readIntensityFromTracks(curTrack,imgMap,1,'reTrack',true,'extraLength',30);
%     curTrack = readIntensityFromTracks(curTrack,tMap,2,'reTrack',false);
%     clear imgMap tMap
%     showAdhesionTracks(pathForColocalization,'all','tracksNA',curTrack)
%  showing
    plotIntensityForceSingle(curTrack,curIndices(curID));
    % tracksNA = readIntensityFromTracks(tracksNA,imgMap,1,'trackOnlyDetected',true);
    %% Edge protrusion for force-transmitting vs. no-transmitting adhesions
    
    edgeMaxProtAll = arrayfun(@(x) x.maxEdgeAdvanceDistChange, tracksNA(curIndices));
    edgeMaxProtForceTrans = edgeMaxProtAll(firstIncreseTimeIntAgainstForceAllIdx);
    edgeMaxProtNonTrans = edgeMaxProtAll(~firstIncreseTimeIntAgainstForceAllIdx);
    meanEdgeMaxProtForceTrans=mean(edgeMaxProtForceTrans)
    meanEdgeMaxProtNonTrans = mean(edgeMaxProtNonTrans)
    [hMaxP,pMaxP]=ttest2(edgeMaxProtForceTrans,edgeMaxProtNonTrans)
    edgeMeanProtAll = arrayfun(@(x) mean(x.edgeAdvanceDist), tracksNA(curIndices));
    edgeMeanProtForceTrans = edgeMeanProtAll(firstIncreseTimeIntAgainstForceAllIdx);
    edgeMeanProtNonTrans = edgeMeanProtAll(~firstIncreseTimeIntAgainstForceAllIdx);
    meanEdgeMeanProtForceTrans=mean(edgeMeanProtForceTrans)
    meanEdgeMeanProtNonTrans=mean(edgeMeanProtNonTrans)
    [hMeanP,pMeanP]=ttest2(edgeMeanProtForceTrans,edgeMeanProtNonTrans)
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'edgeProtrusion.mat'],'edgeMaxProtForceTrans','edgeMaxProtNonTrans','meanEdgeMaxProtForceTrans',...
        'meanEdgeMaxProtNonTrans','hMaxP','pMaxP',...
        'edgeMeanProtForceTrans','edgeMeanProtNonTrans','meanEdgeMeanProtForceTrans','meanEdgeMeanProtNonTrans','hMeanP','pMeanP')
%     save([pathForColocalization filesep 'data' filesep 'tempAllData.mat'])
%% Maturing adhesions' t_init quantification
    curIndicesG2 = find(idGroup2filtered)';
%     curIndicesG2 = idGroup1f';
    for ii=curIndicesG2
        curTrack = tracksNA(ii);
        sFEE = curTrack.startingFrameExtraExtra;
        sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
        bkgMaxInt = max(curTrack.ampTotal(sFEE:sF5before));
        firstIncreseTimeInt = find(curTrack.ampTotal>bkgMaxInt,1);
%         firstIncreseTimeInt = curTrack.startingFrameExtra;
        if ~isempty(firstIncreseTimeInt)
            bkgMaxForce = max(curTrack.forceMag(sFEE:sF5before));
            firstIncreseTimeForce = find(curTrack.forceMag>bkgMaxForce,1);
            if isempty(firstIncreseTimeForce) || firstIncreseTimeForce>curTrack.endingFrameExtra
                tracksNA(ii).forceTransmitting = false;
                tracksNA(ii).firstIncreseTimeInt = [];
                tracksNA(ii).firstIncreseTimeForce = [];
                tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            else
                tracksNA(ii).forceTransmitting = true;
                tracksNA(ii).firstIncreseTimeInt = firstIncreseTimeInt*tInterval;
                tracksNA(ii).firstIncreseTimeForce = firstIncreseTimeForce*tInterval;
                tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeForce*tInterval-firstIncreseTimeInt*tInterval; % +:intensity comes first; -: force comes first.
                tracksNA(ii).bkgMaxInt = bkgMaxInt;
                tracksNA(ii).bkgMaxForce = bkgMaxForce;
            end
        end
    end
    % statistics about firstIncreseTimeIntAgainstForce
    firstIncreseTimeIntAgainstForceAllIdxG2 = arrayfun(@(x) ~isempty(x.forceTransmitting) && ~isempty(x.firstIncreseTimeInt) && x.forceTransmitting,tracksNA(curIndicesG2));
    firstIncreseTimeIntAgainstForceAllG2 = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndicesG2(firstIncreseTimeIntAgainstForceAllIdxG2)));
    figure, histogram(firstIncreseTimeIntAgainstForceAllG2)
    ratioForceTrasG2 = length(firstIncreseTimeIntAgainstForceAllG2)/length(curIndicesG2)
    median(firstIncreseTimeIntAgainstForceAllG2)
    %% save
    save([pathForColocalization filesep 'data' filesep 't_initForMaturingAdhesion.mat'],'firstIncreseTimeIntAgainstForceAllG2','ratioForceTrasG2')

%% Plot edge protrusion distance
%     figure, hold all
%     fileStoreG1_hiCC_edge = [epsPath filesep 'edgeDistPlotG1_hiCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCC_edge,true,'Source', 'edgeAdvanceDist')
%     fileStoreG1_midCC_edge = [epsPath filesep 'edgeDistPlotG1_midCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_midCC_edge,true,'Source', 'edgeAdvanceDist')
%     fileStoreG1_loCC_edge = [epsPath filesep 'edgeDistPlotG1_loCC_shifted.eps'];
%     plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_loCC_edge,true,'Source', 'edgeAdvanceDist')
%     fileStoreG1_loCC_edgeNoshift = [epsPath filesep 'edgeDistPlotG1_loCC.eps'];
%     plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_loCC_edgeNoshift,false,'Source', 'edgeAdvanceDist')
    % ii=idGroup1f(4);
    % fileStoreG1_1 = [epsPath filesep 'ampForcePlotG1_' num2str(ii) '.eps'];
    % plotIntensityForce(tracksNA(ii),fileStoreG1_1)
    %% Backup original tracksNA and save filtered one
    movefile([pathForColocalization filesep 'data' filesep 'tracksNA.mat'], [pathForColocalization filesep 'data' filesep 'tracksNA_org.mat'])
    save([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA','-v7.3')
    %% Backup original tracksNA and save filtered one
    movefile([pathForColocalization filesep 'data' filesep 'idsClassified.mat'], [pathForColocalization filesep 'data' filesep 'idsClassified_org.mat'])
    save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9',...
        'idGroup1filtered','idGroup2filtered','idGroup3filtered','idGroup4filtered','idGroup5filtered','idGroup6filtered','idGroup7filtered','idGroup8filtered','idGroup9filtered','tracksNA','-v7.3')
    %% Show some example tracks
    % get one of median first time lag
    medFirstShift=median(firstIncreseTimeIntAgainstForceAll);
    medFirstShiftIDsAllIdx=(firstIncreseTimeIntAgainstForceAll>=medFirstShift-tInterval) & (firstIncreseTimeIntAgainstForceAll <= medFirstShift+tInterval);
    lifeTimeMedFirstG1 = arrayfun(@(x) x.lifeTime,tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdx)));
    longFirstG1idx = lifeTimeMedFirstG1>20;
    longMedFirstG1idx = longFirstG1idx & medFirstShiftIDsAllIdx;
    longMedFirstG1 = find(longMedFirstG1idx);
    %% loading
    if ~exist('imgMap','var')
        tic
        imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
        imgMap = imgMap.paxImgStack;
        tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
        tMap = tMap.tMap;
        toc
    end
    %% initialization
    curID = 10;
    %% actual showing
    close all
    gPath = [pathForColocalization filesep 'eps'  filesep 'representativeMedian'];
    if ~exist(gPath,'dir')
        mkdir(gPath)
    end
%     medID = longMedFirstG1(51);
%     medID = medID+1;
    additionalName = '';
    curID = curID+1;
    medID = longMedFirstG1(curID);
%     curID = firstIncreseTimeIntAgainstForceAllIdxIDs(56);
    firstIncreseTimeIntAgainstForceAllIdxIDs=find(firstIncreseTimeIntAgainstForceAllIdx);
    curTrack = tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs(medID)));%4210);
    h=showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,curIndices(medID),gPath,additionalName);
    %% Save all workspace
    save([pathForColocalization filesep 'data' filesep 'allDataAfterClassification.mat'],'-v7.3')
    %% Make movie
    close all
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,1500,100)
end
end
