function []=analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath, uptoColo, varargin)
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
ip =inputParser;
ip.addRequired('pathForTheMovieDataFile',@ischar)
ip.addRequired('outputPath',@ischar)
ip.addOptional('uptoColo',false,@islogical)
ip.addParamValue('labeledData',[],@iscell); % This is the master channle index.
ip.parse(pathForTheMovieDataFile,outputPath, uptoColo, varargin{:});
samFolders=ip.Results.labeledData;

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
try
    load([pathForTheMovieDataFile filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat'],'forceField')
catch
    TFMpackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
    forceProc =TFMpackage.processes_{4};
    forceField = forceProc.loadChannelOutput;
end
forceSpacing = forceField(1).pos(2,2)-forceField(1).pos(1,2);
band = 2*forceSpacing;
idClassifiedFile=strcat(dataPath,filesep,'idsClassified.mat');
timeInterval = MD.timeInterval_/60; % in min
% See if you can use existing tracks
if ~exist('tracksNA','var')
    if exist(outputFile,'file') && ~uptoColo  
        disp('tracksNA file is found. Using it ... If you do not want to reuse this, please backup the file and rename it to other name than tracksNA.mat.')
        tic
        if exist(idClassifiedFile,'file') 
            disp('idClassified.mat is also found. Using tracksNA in there...')
            load(idClassifiedFile)
        end
        if ~exist('tracksNA','var')
            tracksNAFile = load(outputFile,'tracksNA');
            tracksNA = tracksNAFile.tracksNA;
        end
        toc
        disp('Done loading tracks.')
    else
        % run analyzeAdhesionMaturation for obtaining tracks from paxillin channel
        tracksNA=colocalizationAdhesionsWithTFM(pathForTheMovieDataFile,band,tmax,false,false,[],...
            'onlyEdge',false,'outputPath',outputPath,'saveAnalysis',true,'matchWithFA',true,'minLifetime',7);
    end
end
%% classification after colocalization
if ~uptoColo
    %% filter out matured adhesion
    %back up tracksNA
%     save(orgTrackFile,'tracksNA','-v7.3');
    pathForColocalization = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
    % load([pathForColocalization filesep 'data'  filesep 'tracksNA.mat'])
%     [tracksNA, ~]=separateMatureAdhesionTracks(tracksNA);
%     tracksNA = filterOutNonEmergingNA(tracksNA);
    iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isfield(tracksNA,'SDC_applied') && ~isempty(iSDCProc)
        disp('Applying stage drift correction ...')
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
        if isa(SDCProc,'EfficientSubpixelRegistrationProcess')
            tracksNA = applyDriftToTracks(tracksNA, T, 1); % need some other function....formatNATracks(tracksNAorg,detectedNAs,nFrames,T); 
        else
            tracksNA = applyDriftToTracks(tracksNA, T, 0);
        end
%         tracksNA = applyDriftToTracks(tracksNA, T); % need some other function....formatNATracks(tracksNAorg,detectedNAs,nFrames,T); 
    else
        disp('Stage drift correction was already applied to tracksNA.')
    end
    % try to pull idGroups - I've already done this before.
%     if exist(strcat(dataPath,filesep,'idsClassified.mat'),'file')
%         idStruct=load(strcat(dataPath,filesep,'idsClassified.mat'));
%         idGroup1=idStruct.idGroup1;
%         idGroup2=idStruct.idGroup2;
%         idGroup3=idStruct.idGroup3;
%         idGroup4=idStruct.idGroup4;
%         idGroup5=idStruct.idGroup5;
%         idGroup6=idStruct.idGroup6;
%         idGroup7=idStruct.idGroup7;
%         idGroup8=idStruct.idGroup8;
%         idGroup9=idStruct.idGroup9;
%         tracksNA=idStruct.tracksNA;
%     else
        [~,idEmerging,idFCstarting] = filterOutNonEmergingNA(tracksNA);
        tracksNA=tracksNA(idEmerging | idFCstarting); % This is to exclude NAs that appear from the very first frame but not FC or FA that can be included in g5, g4, g8
        [~,idFiltered] = filterTracksNAwithMask(pathForTheMovieDataFile,outputPath,tracksNA,band);
        tracksNA = recalculateLifeTimeTracks(tracksNA);
%     end
    %% Classify with ampTotal and distToEdge
%     if ~exist('idGroup1','var')
%     if ~exist(idClassifiedFile,'file') 
        [idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]= ...
            classifyNascentAdhesionTracks(pathForColocalization,'tracksNA',tracksNA,'labeledData',samFolders);
%     end
%     end
    %% filter with idFiltered
    try
        idGroup1filtered = idGroup1 & idFiltered; % & idEmerging;
        idGroup2filtered = idGroup2 & idFiltered;
        idGroup3filtered = idGroup3 & idFiltered; % & idEmerging;
        idGroup4filtered = idGroup4 & idFiltered;
        idGroup5filtered = idGroup5 & idFiltered;
        idGroup6filtered = idGroup6 & idFiltered;
        idGroup7filtered = idGroup7 & idFiltered; % & idEmerging;
        idGroup8filtered = idGroup8 & idFiltered;
        idGroup9filtered = idGroup9 & idFiltered; % & idEmerging;
    catch
        load(idClassifiedFile)
        [~,idFiltered] = filterTracksNAwithMask(pathForTheMovieDataFile,outputPath,tracksNA,band);
        idGroup1filtered = idGroup1 & idFiltered & idEmerging;
        idGroup2filtered = idGroup2 & idFiltered;
        idGroup3filtered = idGroup3 & idFiltered & idEmerging;
        idGroup4filtered = idGroup4 & idFiltered;
        idGroup5filtered = idGroup5 & idFiltered;
        idGroup6filtered = idGroup6 & idFiltered;
        idGroup7filtered = idGroup7 & idFiltered & idEmerging;
        idGroup8filtered = idGroup8 & idFiltered;
        idGroup9filtered = idGroup9 & idFiltered & idEmerging;
    end
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
    dataPath=[pathForColocalization filesep 'data'];
    epsPath=[pathForColocalization filesep 'eps'];
    figPath=[pathForColocalization filesep 'figs'];
    fileStore = [epsPath filesep 'ampForcePlotG1.eps'];
%     plotIntensityForce(tracksNA(idGroup1f),fileStore,false,false)
    plotIntensityForce(tracksNA(idGroup1),fileStore,false,false)
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
    close all
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
    nanmean(ccLagG1Hi)
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
            try
                legend([htrackG2{1}],{'G1 with intermediate CC'},'TextColor','w','Location','best')
            catch
                disp('Nothing to show for G1')
            end
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
%     tracksNA = readIntensityFromTracks(tracksNA,imgMap,1,'reTrack',false,'extraLength',30);
%     tracksNA = readIntensityFromTracks(tracksNA,tMap,2,'reTrack',true);
    
%     % look at 95 % confidence interval of before-detection
%     meanIntenBeforeLevel = mean(curTrack.ampTotal(curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra));
%     stdIntenBeforeLevel = std(curTrack.ampTotal(curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra));
    %% Force confiddence
    % First we need to filter out force-inconfident tracks
%     if ~isfield(tracksNA,'forceUncertainty')
%         tracksNA = readForceUncertaintyFromTracks(pathForTheMovieDataFile,'tracksNA',tracksNA);
%     end
%     % Collect fCfd
%     fCfdG1 = (arrayfun(@(x) mean(x.forceUncertainty(x.startingFrame:x.endingFrame)),tracksNA(idGroup1f)));
%     figure, histogram(fCfdG1)
    % Maybe I can throw away tracks with below 0.3 of force confidence.
    % Then jump directly to figure out whether each track increases its
    % force during its presence compared to the earlier time points
    % clear tMap imgMap
    %% tracks having long-enough pre-detection period will be considered
%     preDetectionPeriod = 20;
%     idGroupLongPreDetecPeriod = arrayfun(@(x) x.startingFrameExtra-x.startingFrameExtraExtra>preDetectionPeriod,tracksNA);
%     preDetecFactor=3/5; %for paxillin because diffuse signal before point-source-like signal
%     preDetecFactor=1/5; %for vinculin 
    preDetecFactor=1/10; %for talin
    %% Initial force-increase time quantification! - time
    tInterval = MD.timeInterval_;
    %% Initial force-increase time quantification! - time
%     curIndices = find(idGroup1filtered | idGroup3filtered | idGroup7filtered | idGroup9filtered)';
%     curIndices = find(idGroup1filtered | idGroup3filtered | idGroup7filtered)';
%     curIndices = find((idGroup1filtered | idGroup3filtered | idGroup7filtered) & idGroupLongPreDetecPeriod)';
%     curIndices = find(idGroup1filtered | idGroup7filtered)';
%     curIndices = find(idGroup7filtered)';
    curIndices = find(idGroup1filtered)';
% %     curIndices = idGroup1f';
%     useSmoothing=true;
    splineParamInit=0.99;
    tempTracks = calculateFirstIncreaseTimeTracks(tracksNA(curIndices),splineParamInit,preDetecFactor,tInterval);
    [tracksNA(curIndices).forceTransmitting] = tempTracks.forceTransmitting;
    [tracksNA(curIndices).firstIncreseTimeInt] = tempTracks.firstIncreseTimeInt;
    [tracksNA(curIndices).firstIncreseTimeForce] = tempTracks.firstIncreseTimeForce;
    [tracksNA(curIndices).firstIncreseTimeIntAgainstForce] = tempTracks.firstIncreseTimeIntAgainstForce;
    [tracksNA(curIndices).bkgMaxInt] = tempTracks.bkgMaxInt;
    [tracksNA(curIndices).bkgMaxForce] = tempTracks.bkgMaxForce;
%     showSingleAdhesionTrackSummary(MD,testTrack,imgMap,tMap)
%     for ii=curIndices
%         curTrack = tracksNA(ii);
%         sFEE = curTrack.startingFrameExtraExtra;
% %         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
%         % See how many frames you have before the startingFrameExtra
%         numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
%         numPreFrames = max(1,floor(preDetecFactor*numFramesBefore));
%         numPreSigStart = min(20,numFramesBefore);
%         sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
%         if useSmoothing
%             d = tracksNA(ii).ampTotal;
%             nTime = length(d);
%             tRange = 1:nTime;
%             numNan = find(isnan(d),1,'last');
%             if isempty(numNan)
%                 numNan=0;
%             end
%             tRange(isnan(d)) = [];
%             d(isnan(d)) = [];
%             sd_spline= csaps(tRange,d,splineParamInit);
%             sd=ppval(sd_spline,tRange);
%             d = [NaN(1,numNan) d];
%     %         tRange = [NaN(1,numNan) tRange];
%             sd = [NaN(1,numNan) sd];
%             bkgMaxInt = max(sd(sFEE:sF5before));
%             firstIncreseTimeInt = find(sd>bkgMaxInt,1);
%         else
%             bkgMaxInt = max(curTrack.ampTotal(sFEE:sF5before));
%             firstIncreseTimeInt = find(curTrack.ampTotal>bkgMaxInt,1);
%         end
% %         firstIncreseTimeInt = curTrack.startingFrameExtra;
%         if ~isempty(firstIncreseTimeInt)
%             if useSmoothing
%                 curForce = tracksNA(ii).forceMag;
%                 sCurForce_spline= csaps(tRange,curForce,splineParamInit);
%                 sCurForce=ppval(sCurForce_spline,tRange);
%                 sCurForce = [NaN(1,numNan) sCurForce];
%                 bkgMaxForce = max(sCurForce(sFEE:sF5before));
%                 firstIncreseTimeForce = find(sCurForce>bkgMaxForce,1);
%             else
%                 bkgMaxForce = max(curTrack.forceMag(sFEE:sF5before));
%                 firstIncreseTimeForce = find(curTrack.forceMag>bkgMaxForce,1);
%             end
%             if isempty(firstIncreseTimeForce) || firstIncreseTimeForce>curTrack.endingFrameExtra
%                 tracksNA(ii).forceTransmitting = false;
%                 tracksNA(ii).firstIncreseTimeInt = [];
%                 tracksNA(ii).firstIncreseTimeForce = [];
%                 tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
%                 tracksNA(ii).bkgMaxInt = bkgMaxInt;
%                 tracksNA(ii).bkgMaxForce = bkgMaxForce;
%             else
%                 tracksNA(ii).forceTransmitting = true;
%                 tracksNA(ii).firstIncreseTimeInt = firstIncreseTimeInt*tInterval; % in sec
%                 tracksNA(ii).firstIncreseTimeForce = firstIncreseTimeForce*tInterval;
%                 tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeInt*tInterval - firstIncreseTimeForce*tInterval; % -:intensity comes first; +: force comes first. in sec
%                 tracksNA(ii).bkgMaxInt = bkgMaxInt;
%                 tracksNA(ii).bkgMaxForce = bkgMaxForce;
%             end
%         else
%             tracksNA(ii).forceTransmitting = false;
%             tracksNA(ii).firstIncreseTimeInt = [];
%             tracksNA(ii).firstIncreseTimeForce = [];
%             tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
%             tracksNA(ii).bkgMaxInt = bkgMaxInt;
%             tracksNA(ii).bkgMaxForce = [];
%         end
%     end
    % statistics about firstIncreseTimeIntAgainstForce
    firstIncreseTimeIntAgainstForceAllIdx = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNA(curIndices));
    firstIncreseTimeIntAgainstForceAll = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdx)));
    numForceTransmitting = length(firstIncreseTimeIntAgainstForceAll);
    numNonTransmitting = sum(~firstIncreseTimeIntAgainstForceAllIdx);
%     figure, histogram(firstIncreseTimeIntAgainstForceAll,-50:2:50)
    figure, histogram(firstIncreseTimeIntAgainstForceAll)
    disp(['Median of firstIncreseTimeIntAgainstForceAll = ' num2str(median(firstIncreseTimeIntAgainstForceAll))])

    %% See ccLag
    close all
    figure, histogram(ccLagG1*tInterval)
    disp(['Median of ccLagG1 = ' num2str(median(ccLagG1*tInterval))])
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'timeInitLagsG1.mat'],'ccLagG1','firstIncreseTimeIntAgainstForceAll','numForceTransmitting','numNonTransmitting')
    %% Now it's time to quantify peak intensity vs. peak force in Group 1
    %   perform spline filter for ampTotal
    splineParam=0.1;
%     for ii=curIndices
    for ii=curIndices
        d = tracksNA(ii).ampTotal;
        tRange = tracksNA(ii).iFrame;
        warning('off','SPLINES:CHCKXYWP:NaNs')
        d(d==0)=NaN;
        try
            sd_spline= csaps(tRange,d,splineParam);
        catch
            d = tracksNA(ii).amp;
            d(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra) = ...
                tracksNA(ii).ampTotal(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra);
            sd_spline= csaps(tRange,d,splineParamInit);
        end
        sd=ppval(sd_spline,tRange);
        sd(isnan(d))=NaN;
        [~,curFrameMaxAmp]=nanmax(sd);
%         timeToMaxInten(ii) = curFrameMaxAmp-tracksNA(ii).startingFrameExtra;

%         d = tracksNA(ii).ampTotal;
%         nTime = length(d);
%         tRange = 1:nTime;
%         numNan = find(isnan(d),1,'last');
%         if isempty(numNan)
%             numNan=0;
%         end
%         tRange(isnan(d)) = [];
%         d(isnan(d)) = [];
%         framesBeforeFII = 0;
%         if ~isempty(tracksNA(ii).forceTransmitting) && tracksNA(ii).forceTransmitting
%             firstIncreseIntFrame = round(tracksNA(ii).firstIncreseTimeInt/tInterval);
%             if firstIncreseIntFrame>tRange(1)
%                 framesBeforeFII = firstIncreseIntFrame-tRange(1);
%             end
%         end
%         sd_spline= csaps(tRange,d,splineParam);
%         sd=ppval(sd_spline,tRange(1+framesBeforeFII:end));
%         d = [NaN(1,numNan) d];
% %         tRange = [NaN(1,numNan) tRange];
%         sd = [NaN(1,numNan+framesBeforeFII) sd];
% %         sd(isnan(d)) = NaN;
%         % Find the maximum
%         [~,curFrameMaxAmp]=max(sd);
        if tracksNA(ii).lifeTime>31
            movingWindowSize=31;
        else
            movingWindowSize=5;
        end
        
        curFrameLocMaxes=locmax1d(sd,movingWindowSize);
        
        if ismember(curFrameMaxAmp,curFrameLocMaxes) && curFrameMaxAmp>tRange(1) && curFrameMaxAmp<tRange(end) && ~isempty(tracksNA(ii).forceTransmitting) && tracksNA(ii).forceTransmitting
            tracksNA(ii).intenPeakness = true;
            tracksNA(ii).intenPeakFrame = curFrameMaxAmp;

            curForce=d;
            curForce(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra) = tracksNA(ii).forceMag(tracksNA(ii).startingFrameExtraExtra:tracksNA(ii).endingFrameExtraExtra);
            
%             curForce = tracksNA(ii).forceMag;
%             curForce(isnan(curForce)) = [];
            sCurForce_spline= csaps(tRange,curForce,splineParam);
%             
%             framesBeforeFTI = 0;
%             if ~isempty(tracksNA(ii).forceTransmitting) && tracksNA(ii).forceTransmitting
%                 firstIncreseTractionFrame = round(tracksNA(ii).firstIncreseTimeForce/tInterval);
%                 if firstIncreseTractionFrame>tRange(1)
%                     framesBeforeFTI = firstIncreseTractionFrame-tRange(1);
%                 end
%             end

            sCurForce_sd=ppval(sCurForce_spline,tRange);
            sCurForce_sd(isnan(curForce))=NaN;
            
%             sCurForce=ppval(sCurForce_spline,tRange(1+framesBeforeFTI:end));
%             sCurForce = [NaN(1,numNan+framesBeforeFTI) sCurForce];
            curForceLocMaxes=locmax1d(sCurForce_sd,movingWindowSize);
            % delete the first and last frame from loc maxes
            curForceLocMaxes = setdiff(curForceLocMaxes, [tracksNA(ii).startingFrameExtra tracksNA(ii).endingFrameExtra]);
            if ~isempty(curForceLocMaxes)
%                 [forceMacMax,indMaxForceAmongLMs] = max(sCurForce(curForceLocMaxes));
                forceMagLMCand=sCurForce_sd(curForceLocMaxes)';
                weightForceMag = (forceMagLMCand - min(sCurForce_sd))/(max(forceMagLMCand) - min(sCurForce_sd));
                [forceMacMax,indMaxForceAmongLMs] = min((abs(curForceLocMaxes-curFrameMaxAmp)).^0.5/length(tracksNA(ii).lifeTime)./weightForceMag);
%                 [forceMacMax,indMaxForceAmongLMs] = min(abs(curForceLocMaxes-curFrameMaxAmp)./weightForceMag);
                tracksNA(ii).forcePeakness = true;
                tracksNA(ii).forcePeakFrame = curForceLocMaxes(indMaxForceAmongLMs);
                tracksNA(ii).forcePeakMag = forceMacMax; 
                tracksNA(ii).forcePeakMagRel = forceMacMax - min(sCurForce_sd); % this is a relative difference
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
    % statistics about peakTimeIntAgainstForceAll
    peakTimeIntAgainstForceAllIdx = arrayfun(@(x) ~isempty(x.intenPeakness) && x.intenPeakness && ~isempty(x.forcePeakness) && x.forcePeakness,tracksNA(curIndices));
    peakTimeIntAgainstForceAll = arrayfun(@(x) x.peakIntTimeIntAgainstForcePeak, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    peakForceAll = arrayfun(@(x) x.forcePeakMag, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    peakForceRelAll = arrayfun(@(x) x.forcePeakMagRel, tracksNA(curIndices(peakTimeIntAgainstForceAllIdx)));
    figure, histogram(peakTimeIntAgainstForceAll)
%     figure, histogram(peakForceAll)
%     figure, histogram(peakForceRelAll)
    disp(['Median of peakTimeIntAgainstForceAll = ' num2str(median(peakTimeIntAgainstForceAll))])
%     median(peakTimeIntAgainstForceAll)
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
    try
        lifeTimeBins = [edgesLT(1) edgesLT(2) edgesLT(3) edgesLT(4) edgesLT(5) edgesLT(6) edgesLT(8) edgesLT(11) edgesLT(end)];
    catch
        lifeTimeBins = edgesLT;
    end
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
%     curID = nonTransmittingIdx(1);       
%     curTrack = tracksNA(curIndices(curID));%4210);
%     plotIntensityForceSingle(curTrack,curIndices(curID));
    %% for transmitting
    firstIncreseTimeIntAgainstForceAllIdxIDs = find(firstIncreseTimeIntAgainstForceAllIdx);
    FTID=ceil(length(firstIncreseTimeIntAgainstForceAllIdxIDs)/2);
    %% Actual plotting
    curID = firstIncreseTimeIntAgainstForceAllIdxIDs(FTID);
    curTrack = tracksNA(curIndices(curID));%4210);
%     curTrack = readIntensityFromTracks(curTrack,imgMap,1,'reTrack',true,'extraLength',30);
%     curTrack = readIntensityFromTracks(curTrack,tMap,2,'reTrack',false);
%     clear imgMap tMap
%     showAdhesionTracks(pathForColocalization,'all','tracksNA',curTrack)
%  showing
%     plotIntensityForceSingle(curTrack,curIndices(curID));
    close all
    figure,set(gcf,'Position',[200,400,400,200]),plotIntensityForceSingle(curTrack,curIndices(curID),0,0.01,splineParamInit,tInterval)
    % tracksNA = readIntensityFromTracks(tracksNA,imgMap,1,'trackOnlyDetected',true);
    FTID=FTID+1;
    %% Edge protrusion for force-transmitting vs. no-transmitting adhesions
    
    edgeMaxProtAll = arrayfun(@(x) x.maxEdgeAdvanceDistChange, tracksNA(curIndices));
    edgeMaxProtForceTrans = edgeMaxProtAll(firstIncreseTimeIntAgainstForceAllIdx);
    edgeMaxProtNonTrans = edgeMaxProtAll(~firstIncreseTimeIntAgainstForceAllIdx);
    edgeMaxProtCell={edgeMaxProtForceTrans,edgeMaxProtNonTrans};
    namesList={'Force-Transmitting','Non-transmitting'};
    h1=figure;
    barPlotCellArray(edgeMaxProtCell,namesList); title('Max Protrusion'); ylabel('Max protrusion (um)')
    hgexport(h1,strcat(epsPath,'/MaxProtrusionForceTran'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/MaxProtrusionForceTran'),'-v7.3')
    
    meanEdgeMaxProtForceTrans=mean(edgeMaxProtForceTrans);
    meanEdgeMaxProtNonTrans = mean(edgeMaxProtNonTrans);
    [hMaxP,pMaxP]=ttest2(edgeMaxProtForceTrans,edgeMaxProtNonTrans)
    edgeMeanProtAll = arrayfun(@(x) mean(x.edgeAdvanceDist), tracksNA(curIndices));
    edgeMeanProtForceTrans = edgeMeanProtAll(firstIncreseTimeIntAgainstForceAllIdx);
    edgeMeanProtNonTrans = edgeMeanProtAll(~firstIncreseTimeIntAgainstForceAllIdx);
    meanEdgeMeanProtForceTrans=mean(edgeMeanProtForceTrans);
    meanEdgeMeanProtNonTrans=mean(edgeMeanProtNonTrans);
    [hMeanP,pMeanP]=ttest2(edgeMeanProtForceTrans,edgeMeanProtNonTrans);
    edgeMaxProtCell={edgeMeanProtForceTrans,edgeMeanProtNonTrans};
    namesList={'Force-Transmitting','Non-transmitting'};
    h1=figure;
    barPlotCellArray(edgeMaxProtCell,namesList); title('Mean Protrusion'); ylabel('Mean protrusion (um)')
    hgexport(h1,strcat(epsPath,'/MeanProtrusionForceTran'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/MeanProtrusionForceTran'),'-v7.3')
    %% Save output
    save([pathForColocalization filesep 'data' filesep 'edgeProtrusion.mat'],'edgeMaxProtForceTrans','edgeMaxProtNonTrans','meanEdgeMaxProtForceTrans',...
        'meanEdgeMaxProtNonTrans','hMaxP','pMaxP',...
        'edgeMeanProtForceTrans','edgeMeanProtNonTrans','meanEdgeMeanProtForceTrans','meanEdgeMeanProtNonTrans','hMeanP','pMeanP')
%     save([pathForColocalization filesep 'data' filesep 'tempAllData.mat'])
%% Maturing adhesions' t_init quantification
%     curIndicesG2 = find(idGroup2filtered)';
    curIndicesG2 = idGroup2f;
    tempTracks2 = calculateFirstIncreaseTimeTracks(tracksNA(curIndicesG2),splineParamInit,preDetecFactor,tInterval);
    [tracksNA(curIndicesG2).forceTransmitting] = tempTracks2.forceTransmitting;
    [tracksNA(curIndicesG2).firstIncreseTimeInt] = tempTracks2.firstIncreseTimeInt;
    [tracksNA(curIndicesG2).firstIncreseTimeForce] = tempTracks2.firstIncreseTimeForce;
    [tracksNA(curIndicesG2).firstIncreseTimeIntAgainstForce] = tempTracks2.firstIncreseTimeIntAgainstForce;
    [tracksNA(curIndicesG2).bkgMaxInt] = tempTracks2.bkgMaxInt;
    [tracksNA(curIndicesG2).bkgMaxForce] = tempTracks2.bkgMaxForce;
%     for ii=curIndicesG2
%         curTrack = tracksNA(ii);
%         sFEE = curTrack.startingFrameExtraExtra;
%         numFramesBefore = curTrack.startingFrameExtra - curTrack.startingFrameExtraExtra;
%         numPreFrames = max(1,floor(preDetecFactor*numFramesBefore)); % for maturing adhesions, we expect there can be ealier growth that's not detected. (10/29/15)
%         numPreSigStart = min(20,numFramesBefore);
%         sF5before = max(curTrack.startingFrameExtra-numPreSigStart,curTrack.startingFrameExtra-numPreFrames);
% %         sF5before = max(curTrack.startingFrameExtraExtra,curTrack.startingFrameExtra-5);
%         bkgMaxInt = max(curTrack.ampTotal(sFEE:sF5before));
%         firstIncreseTimeInt = find(curTrack.ampTotal>bkgMaxInt,1);
% %         firstIncreseTimeInt = curTrack.startingFrameExtra;
%         if ~isempty(firstIncreseTimeInt)
%             bkgMaxForce = max(curTrack.forceMag(sFEE:sF5before));
%             firstIncreseTimeForce = find(curTrack.forceMag>bkgMaxForce,1);
%             if isempty(firstIncreseTimeForce) || firstIncreseTimeForce>curTrack.endingFrameExtra
%                 tracksNA(ii).forceTransmitting = false;
%                 tracksNA(ii).firstIncreseTimeInt = [];
%                 tracksNA(ii).firstIncreseTimeForce = [];
%                 tracksNA(ii).firstIncreseTimeIntAgainstForce = []; 
%                 tracksNA(ii).bkgMaxInt = bkgMaxInt;
%                 tracksNA(ii).bkgMaxForce = bkgMaxForce;
%             else
%                 tracksNA(ii).forceTransmitting = true;
%                 tracksNA(ii).firstIncreseTimeInt = firstIncreseTimeInt*tInterval;
%                 tracksNA(ii).firstIncreseTimeForce = firstIncreseTimeForce*tInterval;
%                 tracksNA(ii).firstIncreseTimeIntAgainstForce = firstIncreseTimeForce*tInterval-firstIncreseTimeInt*tInterval; % +:intensity comes first; -: force comes first.
%                 tracksNA(ii).bkgMaxInt = bkgMaxInt;
%                 tracksNA(ii).bkgMaxForce = bkgMaxForce;
%             end
%         end
%     end
    % statistics about firstIncreseTimeIntAgainstForce
    firstIncreseTimeIntAgainstForceAllIdxG2 = arrayfun(@(x) ~isempty(x.forceTransmitting) && ~isempty(x.firstIncreseTimeInt) && x.forceTransmitting,tracksNA(curIndicesG2));
    firstIncreseTimeIntAgainstForceAllG2 = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndicesG2(firstIncreseTimeIntAgainstForceAllIdxG2)));
    figure, histogram(firstIncreseTimeIntAgainstForceAllG2,20)
    ratioForceTrasG2 = length(firstIncreseTimeIntAgainstForceAllG2)/length(curIndicesG2)
    disp(['Median of firstIncreseTimeIntAgainstForceAllG2 = ' num2str(median(firstIncreseTimeIntAgainstForceAllG2))])
    FTID2 = 1;
    %% save
    save([pathForColocalization filesep 'data' filesep 't_initForMaturingAdhesion.mat'],'firstIncreseTimeIntAgainstForceAllG2','ratioForceTrasG2')
    %% Example for Group2
%     curID2 = firstIncreseTimeIntAgainstForceAllIdxG2(FTID2);
    try
        curTrack2 = tracksNA(curIndicesG2(FTID2));%4210);
    %     curTrack = readIntensityFromTracks(curTrack,imgMap,1,'reTrack',true,'extraLength',30);
    %     curTrack = readIntensityFromTracks(curTrack,tMap,2,'reTrack',false);
    %     clear imgMap tMap
    %     showAdhesionTracks(pathForColocalization,'all','tracksNA',curTrack)
    %  showing
    %     plotIntensityForceSingle(curTrack,curIndices(curID));
        close all
        figure,set(gcf,'Position',[200,400,400,200]),plotIntensityForceSingle(curTrack2,curIndicesG2(FTID2),0,0.01,splineParamInit,tInterval)
        % tracksNA = readIntensityFromTracks(tracksNA,imgMap,1,'trackOnlyDetected',true);
        FTID2=FTID2+1;
    catch
        disp(['No adhesion in G2 in ' num2str(FTID2) 'th. Skipping...'])
    end
    %% average time from t_init to t_peak
    firstIncreseTimeForce_both = arrayfun(@(x) x.firstIncreseTimeForce, tracksNA(curIndices(bothTiTpIdx)));
    peakTimeForceAll_both = arrayfun(@(x) x.forcePeakFrame*tInterval, tracksNA(curIndices(bothTiTpIdx)));
    timeToPeakForceFT = peakTimeForceAll_both-firstIncreseTimeForce_both;
%     figure, histogram(timeToPeakForceFT)
%     figure, plot(firstIncreseTimeForce_both,peakTimeForceAll_both,'.')
    disp(['Median of timeToPeakForceFT = ' num2str(median(timeToPeakForceFT))])
%     median(timeToPeakForceFT)

    firstIncreseTimeInt_both = arrayfun(@(x) x.firstIncreseTimeInt, tracksNA(curIndices(bothTiTpIdx)));
    peakTimeIntAll_both = arrayfun(@(x) x.intenPeakFrame*tInterval, tracksNA(curIndices(bothTiTpIdx)));
    timeToPeakIntFT = peakTimeIntAll_both-firstIncreseTimeInt_both;
%     figure, histogram(timeToPeakIntFT)
%     figure, plot(firstIncreseTimeInt_both,peakTimeIntAll_both,'.')
    disp(['Median of timeToPeakIntFT = ' num2str(median(timeToPeakIntFT))])
%     median(timeToPeakIntFT)

    save([pathForColocalization filesep 'data' filesep 'timeToPeaks.mat'],'timeToPeakForceFT','timeToPeakIntFT')
    
    %% What is the initial time lags that have peaks? - just to check
    tLagInitial_both = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, tracksNA(curIndices(bothTiTpIdx)));
    disp(['Median of tLagInitial_both = ' num2str(median(tLagInitial_both))])
%     median(tLagInitial_both)
    tLagPeak_both = arrayfun(@(x) x.peakIntTimeIntAgainstForcePeak, tracksNA(curIndices(bothTiTpIdx)));
    disp(['Median of tLagPeak_both = ' num2str(median(tLagPeak_both))])
%     median(tLagPeak_both)
    % Save output
%     save([pathForColocalization filesep 'data' filesep 'timeInitLagsG1.mat'],'ccLagG1','firstIncreseTimeIntAgainstForceAll','numForceTransmitting','numNonTransmitting')
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
    %% Show some example tracks with median init time lag and peak time lag
    % get one of median first time lag
    medFirstShift=median(firstIncreseTimeIntAgainstForce_both);
    medPeakShift=median(peakTimeIntAgainstForceAll_both);
%     medFirstShift = medFirstShift;
%     medFirstShift = medFirstShift-10*tInterval;
%     medFirstShiftIDsAllIdx=(firstIncreseTimeIntAgainstForceAll>medFirstShift-tInterval) & (firstIncreseTimeIntAgainstForceAll < medFirstShift+tInterval);
%     medInitAndPeakAllIdx=(firstIncreseTimeIntAgainstForce_both>medFirstShift-3*tInterval) & ...
%         (firstIncreseTimeIntAgainstForce_both < medFirstShift+2*tInterval) & ...
%         (peakTimeIntAgainstForceAll_both>medPeakShift-20*tInterval) & ...
%         (peakTimeIntAgainstForceAll_both < medPeakShift+5*tInterval);
%     iMedInitPeak = find(medInitAndPeakAllIdx);
    medInitAndPeakAllIdx=(firstIncreseTimeIntAgainstForce_both>medFirstShift-3*tInterval) & ...
        (firstIncreseTimeIntAgainstForce_both < medFirstShift+2*tInterval);
    iMedInitPeak = find(medInitAndPeakAllIdx);
%     medFirstShiftIDsAllIdx=(firstIncreseTimeIntAgainstForceAll==medFirstShift);
%     lifeTimeMedFirstG1 = arrayfun(@(x) x.lifeTime,tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdx)));
%     longFirstG1idx = lifeTimeMedFirstG1>20;
%     longMedFirstG1idx = longFirstG1idx & medFirstShiftIDsAllIdx;
%     longMedFirstG1 = find(longMedFirstG1idx);
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
    curID = 0;
    %% actual showing
    close all
    gPath = [pathForColocalization filesep 'eps'  filesep 'representativeMedian'];
    if ~exist(gPath,'dir')
        mkdir(gPath)
    end
%     medID = longMedFirstG1(51);
%     medID = medID+1;
    additionalName = '';
%     additionalName = 'forcePeakAfterFluorLT';
%     additionalName = 'forcePeakSame';
%     additionalName = 'forcePeakLater';
%     additionalName = 'forcePeakEarly';
    if isempty(additionalName)
        curID = curID+1;
    end
    medID = iMedInitPeak(curID);
%     curID = firstIncreseTimeIntAgainstForceAllIdxIDs(56);
%     firstIncreseTimeIntAgainstForceAllIdxIDs=find(firstIncreseTimeIntAgainstForceAllIdx);
%     curTrack = tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs(medID)));%4210);
%     h=showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,curIndices(medID),gPath,additionalName);
    bothTiTpIdxIDs=find(bothTiTpIdx);
    curTrack = tracksNA(curIndices(bothTiTpIdxIDs(medID)));%4210);
    showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,curIndices(bothTiTpIdxIDs(medID)),gPath,additionalName);
    if ~isempty(additionalName)
        save([gPath filesep 'track' num2str(curIndices(medID)) additionalName '.mat'],'curTrack','medID','curIndices','gPath','additionalName')
    end
    %% actual showing for G2 -initialization
    curID2=0;
    %% actual showing for G2
    close all
    gPath = [pathForColocalization filesep 'eps'  filesep 'representativeMedian'];
    if ~exist(gPath,'dir')
        mkdir(gPath)
    end
%     medID = longMedFirstG1(51);
%     medID = medID+1;
    additionalName = '';
%     additionalName = 'forcePeakAfterFluorLT';
%     additionalName = 'forcePeakSame';
%     additionalName = 'forcePeakLater';
%     additionalName = 'forcePeakEarly';
    if isempty(additionalName)
        curID2 = curID2+1;
    end
    try
        medID2 = curIndicesG2(curID2);
%     curID = firstIncreseTimeIntAgainstForceAllIdxIDs(56);
%     firstIncreseTimeIntAgainstForceAllIdxIDs=find(firstIncreseTimeIntAgainstForceAllIdx);
%     curTrack = tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs(medID)));%4210);
%     h=showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,curIndices(medID),gPath,additionalName);
        curTrack2 = tracksNA(medID2);%4210);
        showSingleAdhesionTrackSummary(MD,curTrack2,imgMap,tMap,medID2,gPath,additionalName);
        if ~isempty(additionalName)
            save([gPath filesep 'track' num2str(medID2) additionalName '.mat'],'curTrack2','medID2','curIndices2','gPath','additionalName')
        end
    catch
        disp('no G2 example...')
    end
    %% Look at feature difference per each group
    pixSize=MD.pixelSize_/1000; % in um
    distToEdge{1} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup1filtered));
    distToEdge{2} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup2filtered));
    distToEdge{3} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup3filtered));
    distToEdge{4} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup4filtered));
    distToEdge{5} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup5filtered));
    distToEdge{6} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup6filtered));
    distToEdge{7} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup7filtered));
    distToEdge{8} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup8filtered));
    distToEdge{9} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup9filtered));
    groupNameList={'g1','g2','g3','g4','g5','g6','g7','g8','g9'};
    figure;
    boxPlotCellArray(distToEdge,groupNameList,pixSize)
%     [lengthLongest]=max(cellfun(@(x) length(x),distToEdge));
%     matrixDistToEdge = NaN(lengthLongest,9);
%     for ii=1:9
%         matrixDistToEdge(1:length(distToEdge{ii}),ii) = distToEdge{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{['G1' '(N=' num2str(length(distToEdge{1})) ')'],['G2 (N=' num2str(length(distToEdge{2})) ')'],...
%         ['G3 (N=' num2str(length(distToEdge{3})) ')'],['G4 (N=' num2str(length(distToEdge{4})) ')'],...
%         ['G5 (N=' num2str(length(distToEdge{5})) ')'],['G6 (N=' num2str(length(distToEdge{6})) ')'],...
%         ['G7 (N=' num2str(length(distToEdge{7})) ')'],['G8 (N=' num2str(length(distToEdge{8})) ')'],...
%         ['G9 (N=' num2str(length(distToEdge{9})) ')']},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{sprintf('G%d\n(N=%d)',1, length(distToEdge{1})),sprintf('G%d\n(N=%d)',2, length(distToEdge{2})),...
%         sprintf('G%d\n(N=%d)',3, length(distToEdge{3})),sprintf('G%d\n(N=%d)',4, length(distToEdge{4})),...
%         sprintf('G%d\n(N=%d)',5, length(distToEdge{5})),sprintf('G%d\n(N=%d)',6, length(distToEdge{6})),...
%         sprintf('G%d\n(N=%d)',7, length(distToEdge{7})),sprintf('G%d\n(N=%d)',8, length(distToEdge{8})),...
%         sprintf('G%d\n(N=%d)',9, length(distToEdge{9}))},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('Distance to edge')
    ylabel('Distance to edge (um)')
    save([pathForColocalization filesep 'data' filesep 'distToEdge.mat'],'distToEdge','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'distToEdgeForAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/distToEdgeForAllGroups'),'-v7.3')
    
    %% Look at feature difference per each group - advanceDist
    advanceDist{1} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup1filtered));
    advanceDist{2} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup2filtered));
    advanceDist{3} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup3filtered));
    advanceDist{4} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup4filtered));
    advanceDist{5} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup5filtered));
    advanceDist{6} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup6filtered));
    advanceDist{7} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup7filtered));
    advanceDist{8} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup8filtered));
    advanceDist{9} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup9filtered));

    figure;
    boxPlotCellArray(advanceDist,groupNameList)

%     [lengthLongestAdvDist]=max(cellfun(@(x) length(x),advanceDist));
% 
%     matrixAdvanceDist = NaN(lengthLongestAdvDist,9);
%     for ii=1:9
%         matrixAdvanceDist(1:length(advanceDist{ii}),ii) = advanceDist{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixAdvanceDist,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('Adhesion advancement forward')
    ylabel('Adhesion advancement (um)')
    save([pathForColocalization filesep 'data' filesep 'advanceDist.mat'],'advanceDist','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'advanceDistAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/advanceDistAllGroups'),'-v7.3')
    %% Look at feature difference per each group - ampTotal
    ampTotal{1} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup1filtered));
    ampTotal{2} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup2filtered));
    ampTotal{3} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup3filtered));
    ampTotal{4} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup4filtered));
    ampTotal{5} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup5filtered));
    ampTotal{6} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup6filtered));
    ampTotal{7} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup7filtered));
    ampTotal{8} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup8filtered));
    ampTotal{9} =arrayfun(@(x) nanmean(x.ampTotal),tracksNA(idGroup9filtered));
    figure;
    boxPlotCellArray(ampTotal,groupNameList)
    
%     [lengthLongestAmpTotal]=max(cellfun(@(x) length(x),ampTotal));
% 
%     matrixAmpTotal = NaN(lengthLongestAmpTotal,9);
%     for ii=1:9
%         matrixAmpTotal(1:length(ampTotal{ii}),ii) = ampTotal{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixAmpTotal,'orientation','vertical','whisker',0.5,'notch','off',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
% %     ylim([-2 50])
    title('ampTotal')
    ylabel('Fluorescence intensity (A.U.)')
    save([pathForColocalization filesep 'data' filesep 'ampTotal.mat'],'ampTotal','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'ampTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/ampTotalAllGroups'),'-v7.3')
    %% Look at feature difference per each group - starting ampTotal
    startingAmpTotal{1} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup1filtered));
    startingAmpTotal{2} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup2filtered));
    startingAmpTotal{3} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup3filtered));
    startingAmpTotal{4} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup4filtered));
    startingAmpTotal{5} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup5filtered));
    startingAmpTotal{6} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup6filtered));
    startingAmpTotal{7} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup7filtered));
    startingAmpTotal{8} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup8filtered));
    startingAmpTotal{9} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup9filtered));
    figure;
    boxPlotCellArray(startingAmpTotal,groupNameList)
%     [lengthLongestStartingAmpTotal]=max(cellfun(@(x) length(x),startingAmpTotal));
% 
%     matrixStartingAmpTotal = NaN(lengthLongestStartingAmpTotal,9);
%     for ii=1:9
%         matrixStartingAmpTotal(1:length(startingAmpTotal{ii}),ii) = startingAmpTotal{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixStartingAmpTotal,'orientation','vertical','whisker',0.5,'notch','off',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
% %     ylim([-2 50])
    title('Starting Amplitude')
    ylabel('Fluorescence intensity (A.U.)')
    save([pathForColocalization filesep 'data' filesep 'startingAmpTotal.mat'],'startingAmpTotal','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'startingAmpTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/startingAmpTotalAllGroups'),'-v7.3')
    %% Look at feature difference per each group - starting edgeAdvanceDistChange
    edgeAdvanceDistChange{1} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup1filtered));
    edgeAdvanceDistChange{2} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup2filtered));
    edgeAdvanceDistChange{3} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup3filtered));
    edgeAdvanceDistChange{4} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup4filtered));
    edgeAdvanceDistChange{5} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup5filtered));
    edgeAdvanceDistChange{6} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup6filtered));
    edgeAdvanceDistChange{7} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup7filtered));
    edgeAdvanceDistChange{8} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup8filtered));
    edgeAdvanceDistChange{9} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup9filtered));
%     [lengthLongestStartingAmpTotal]=max(cellfun(@(x) length(x),edgeAdvanceDistChange));
% 
%     matrixEdgeAdvanceDistChange = NaN(lengthLongestStartingAmpTotal,9);
%     for ii=1:9
%         matrixEdgeAdvanceDistChange(1:length(edgeAdvanceDistChange{ii}),ii) = edgeAdvanceDistChange{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixEdgeAdvanceDistChange,'orientation','vertical','whisker',0.5,'notch','off',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
% %     ylim([-2 50])
    figure;
    boxPlotCellArray(edgeAdvanceDistChange,groupNameList)
    title('edgeAdvanceDistChange at the end of tracks (to see g7 has nearly zero edge advance)')
    ylabel('edgeAdvanceDistChange (um)')
    save([pathForColocalization filesep 'data' filesep 'edgeAdvanceDistChange.mat'],'edgeAdvanceDistChange','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'edgeAdvanceDistChangeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/edgeAdvanceDistChangeAllGroups'),'-v7.3')
      %% Look at feature difference per each group - starting forceMag
    startingForceMag{1} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup1filtered));
    startingForceMag{2} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup2filtered));
    startingForceMag{3} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup3filtered));
    startingForceMag{4} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup4filtered));
    startingForceMag{5} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup5filtered));
    startingForceMag{6} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup6filtered));
    startingForceMag{7} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup7filtered));
    startingForceMag{8} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup8filtered));
    startingForceMag{9} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup9filtered));
%     [lengthLongestStartingForceMag]=max(cellfun(@(x) length(x),startingForceMag));
% 
%     matrixStartingForceMag = NaN(lengthLongestStartingForceMag,9);
%     for ii=1:9
%         matrixStartingForceMag(1:length(startingForceMag{ii}),ii) = startingForceMag{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixStartingForceMag,'orientation','vertical','whisker',0.5,'notch','off',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    figure;
    boxPlotCellArray(startingForceMag,groupNameList)
    title('startingForceMag (to see g1,2,3,7 start nearly at similar force)')
    ylabel('startingForceMag (Pa)')
    save([pathForColocalization filesep 'data' filesep 'startingForceMag.mat'],'startingForceMag','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'startingForceMagAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    hgsave(strcat(figPath,'/startingForceMagAllGroups'),'-v7.3')
    %% recalculate force slope 
    tracksNA = calculateTrackSlopes(tracksNA,tInterval);
    
%     tIntervalMin = tInterval/60;
%     for k=1:numel(tracksG1)
%         sF=tracksG1(k).startingFrameExtra;
%         eF=tracksG1(k).endingFrameExtra;
%         curLT = tracksG1(k).lifeTime;
%         halfLT = ceil(curLT/2);
%         earlyPeriod = min(halfLT,floor(60/tInterval)); % frames per a minute or half life time
% 
%         lastFrame = min(sum(~isnan(tracksG1(k).ampTotal)),sF+earlyPeriod-1);
%         lastFrameFromOne = lastFrame - sF+1;
%     %     [curR,curM] = regression(timeInterval*(1:lastFrameFromOne),tracksNA(k).amp(tracksNA(k).startingFrame:lastFrame));
%     %     tracksNA(k).ampSlope = curM; % in a.u./min
%     %     tracksNA(k).ampSlopeR = curR; % Pearson's correlation coefficient
%         [curForceR,curForceM] = regression(tIntervalMin*(1:lastFrameFromOne),tracksG1(k).forceMag(sF:lastFrame));
% %         figure, plot(tIntervalMin*(1:lastFrameFromOne),tracksG1(k).forceMag(sF:lastFrame))
%         tracksG1(k).forceSlope = curForceM; % in Pa/min
%         tracksG1(k).forceSlopeR = curForceR; % Pearson's correlation coefficient
%     end
%     plotIntensityForceSingle(tracksG1(k),1,0,0.01,0.3,tInterval)
        %% assembly rate and force growth rate in curIndices
    forceSlopeG1 =arrayfun(@(x) (x.forceSlope),tracksNA(curIndices));
    earlyAmpSlopeG1 =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(curIndices));
   
    save([pathForColocalization filesep 'data' filesep 'assemblyRateForceSlopes.mat'],'forceSlopeG1','earlyAmpSlopeG1')
    
    %% non force transmitting adhesion statistics
    % average fluorescence intensity
    avgFluoInten_nonTransmitting = arrayfun(@(x) nanmean(x.ampTotal(logical(x.presence))),tracksNA(curIndices(nonTransmittingIdx)));
    nanmean(avgFluoInten_nonTransmitting);
    avgFluoInten_forceTransmitting = arrayfun(@(x) nanmean(x.ampTotal(logical(x.presence))),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    nanmean(avgFluoInten_forceTransmitting);
    [hFI,pFI]=ttest2(avgFluoInten_nonTransmitting,avgFluoInten_forceTransmitting);
    avgFluoIntenCell={avgFluoInten_forceTransmitting,avgFluoInten_nonTransmitting};
%     namesList={'Force-Transmitting','Non-transmitting'};
    h1=figure;
    barPlotCellArray(avgFluoIntenCell,namesList); title('avgFluoIntenCell'); ylabel('avgFluoIntenCell (a.u.)')
    hgexport(h1,strcat(epsPath,'/avgFluoInten'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/avgFluoInten'),'-v7.3')
   
    % only for pre-detection period
    avgPreDetecFluoInten_nonTransmitting = arrayfun(@(x) nanmean(x.ampTotal(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(nonTransmittingIdx)));
    nanmean(avgPreDetecFluoInten_nonTransmitting)
    avgPreDetecFluoInten_forceTransmitting = arrayfun(@(x) nanmean(x.ampTotal(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    nanmean(avgPreDetecFluoInten_forceTransmitting)
    [hPreDetecFI,pPreDetecFI]=ttest2(avgPreDetecFluoInten_forceTransmitting,avgPreDetecFluoInten_nonTransmitting);
    
    avgPreDetecFluoIntenCell={avgPreDetecFluoInten_forceTransmitting,avgPreDetecFluoInten_nonTransmitting};
    h1=figure;
    barPlotCellArray(avgPreDetecFluoIntenCell,namesList); title('avgPreDetecFluoIntenCell'); ylabel('avgPreDetecFluoIntenCell (a.u.)')
    hgexport(h1,strcat(epsPath,'/avgPreDetecFluoInten'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/avgPreDetecFluoInten'),'-v7.3')
    % earlyAmpSlope
    avgEarlyAmpSlope_nonTransmitting = arrayfun(@(x) x.earlyAmpSlope,tracksNA(curIndices(nonTransmittingIdx)));
    nanmean(avgEarlyAmpSlope_nonTransmitting)
    avgEarlyAmpSlope_forceTransmitting = arrayfun(@(x) x.earlyAmpSlope,tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    nanmean(avgEarlyAmpSlope_forceTransmitting)
    [hEarlyAmpSlope,pEarlyAmpSlope]=ttest2(avgEarlyAmpSlope_nonTransmitting,avgEarlyAmpSlope_forceTransmitting)
    avgEarlyAmpSlopeCell={avgEarlyAmpSlope_forceTransmitting,avgEarlyAmpSlope_nonTransmitting};
    h1=figure;
    barPlotCellArray(avgEarlyAmpSlopeCell,namesList); title('avgEarlyAmpSlopeCell'); ylabel('avgEarlyAmpSlopeCell (a.u./min)')
    hgexport(h1,strcat(epsPath,'/avgEarlyAmpSlope'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/avgEarlyAmpSlope'),'-v7.3')
    % forceMag for pre-detection period
    avgPreDetecForceMag_nonTransmitting = arrayfun(@(x) nanmean(x.forceMag(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(nonTransmittingIdx)));
    nanmean(avgPreDetecForceMag_nonTransmitting)
    avgPreDetecForceMag_forceTransmitting = arrayfun(@(x) nanmean(x.forceMag(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    nanmean(avgPreDetecForceMag_forceTransmitting)
    [hPreDetecForce,pPreDetecForce]=ttest2(avgPreDetecForceMag_nonTransmitting,avgPreDetecForceMag_forceTransmitting)
    avgPreDetecForceMagCell={avgPreDetecForceMag_forceTransmitting,avgPreDetecForceMag_nonTransmitting};
    h1=figure;
    barPlotCellArray(avgPreDetecForceMagCell,namesList); title('avgPreDetecForceMagCell'); ylabel('avgPreDetecForceMagCell (Pa)')
    hgexport(h1,strcat(epsPath,'/avgPreDetecForceMag'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/avgPreDetecForceMag'),'-v7.3')
    % forceSlope
    avgForceSlope_nonTransmitting = arrayfun(@(x) x.forceSlope,tracksNA(curIndices(nonTransmittingIdx)));
    nanmean(avgForceSlope_nonTransmitting)
    avgForceSlope_forceTransmitting = arrayfun(@(x) x.forceSlope,tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    nanmean(avgForceSlope_forceTransmitting)
    [hForceSlope,pForceSlope]=ttest2(avgForceSlope_nonTransmitting,avgForceSlope_forceTransmitting)
    avgForceSlopeCell={avgForceSlope_forceTransmitting,avgForceSlope_nonTransmitting};
    h1=figure;
    boxPlotCellArray(avgForceSlopeCell,namesList); title('avgForceSlopeCell'); ylabel('avgForceSlopeCell (Pa/min)')
    hgexport(h1,strcat(epsPath,'/avgForceSlope'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/avgForceSlope'),'-v7.3')
    %% Save these parameters
    save([pathForColocalization filesep 'data' filesep 'nonTransmittingVsForceTransmitting.mat'],'avgFluoInten_nonTransmitting',...
        'avgFluoInten_forceTransmitting','hFI','pFI','avgPreDetecFluoInten_nonTransmitting','avgPreDetecFluoInten_forceTransmitting',...
        'hPreDetecFI','pPreDetecFI','avgEarlyAmpSlope_nonTransmitting','avgEarlyAmpSlope_forceTransmitting',...
        'hEarlyAmpSlope','pEarlyAmpSlope','avgPreDetecForceMag_nonTransmitting','avgPreDetecForceMag_forceTransmitting',...
        'hPreDetecForce','pPreDetecForce','avgForceSlope_forceTransmitting','avgForceSlope_nonTransmitting','hForceSlope','pForceSlope','-v7.3')
%% Look at feature difference per each group - earlyAmpSlope
    earlyAmpSlope{1} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup1filtered));
    earlyAmpSlope{2} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup2filtered));
    earlyAmpSlope{3} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup3filtered));
    earlyAmpSlope{4} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup4filtered));
    earlyAmpSlope{5} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup5filtered));
    earlyAmpSlope{6} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup6filtered));
    earlyAmpSlope{7} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup7filtered));
    earlyAmpSlope{8} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup8filtered));
    earlyAmpSlope{9} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup9filtered));
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
    figure;
    boxPlotCellArray(earlyAmpSlope,groupNameList)
    title('earlyAmpSlope (this shows that g7 has the same early slope as g3).')
    ylabel('earlyAmpSlope (A.U./min)')
    hgsave(strcat(figPath,'/earlyAmpSlopeAllGroups'),'-v7.3')
    save([pathForColocalization filesep 'data' filesep 'earlyAmpSlopeAllGroups.mat'],'earlyAmpSlope','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'earlyAmpSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - force slope
    forceSlope{1} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup1filtered));
    forceSlope{2} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup2filtered));
    forceSlope{3} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup3filtered));
    forceSlope{4} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup4filtered));
    forceSlope{5} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup5filtered));
    forceSlope{6} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup6filtered));
    forceSlope{7} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup7filtered));
    forceSlope{8} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup8filtered));
    forceSlope{9} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup9filtered));
%     [lengthForceSlope]=max(cellfun(@(x) length(x),forceSlope));
% 
%     matrixForceSlope = NaN(lengthForceSlope,9);
%     for ii=1:9
%         matrixForceSlope(1:length(forceSlope{ii}),ii) = forceSlope{ii};
%     end
%     boxWidth=0.5;
%     figure
%     boxplot(matrixForceSlope,'orientation','vertical','whisker',1,'notch','off',...
%         'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'colors','k')
%     set(findobj(gca,'LineStyle','--'),'LineStyle','-')
%     set(findobj(gca,'tag','Median'),'LineWidth',2)
% %     ylim([-2 50])
    figure;
    boxPlotCellArray(forceSlope,groupNameList)
    title('forceSlope (all nascent adhesions (g1,2,3,7) show the same force slopes).')
    ylabel('forceSlope (Pa/min)')
    hgsave(strcat(figPath,'/forceSlopeAllGroups'),'-v7.3')
    save([pathForColocalization filesep 'data' filesep 'forceSlopeAllGroups.mat'],'forceSlope','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'forceSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% export tracksG1, G2, G3 and G7 separately
    tracksG1 = tracksNA(idGroup1filtered);
    tracksG2 = tracksNA(idGroup2filtered);
    tracksG3 = tracksNA(idGroup3filtered);
    tracksG7 = tracksNA(idGroup7filtered);
    save([pathForColocalization filesep 'data' filesep 'tracksG1.mat'],'tracksG1','-v7.3')
    save([pathForColocalization filesep 'data' filesep 'tracksG2.mat'],'tracksG2','-v7.3')
    save([pathForColocalization filesep 'data' filesep 'tracksG3.mat'],'tracksG3','-v7.3')
    save([pathForColocalization filesep 'data' filesep 'tracksG7.mat'],'tracksG7','-v7.3')
    %% Filtering G7 for those that increase edge advance 
    timeInterval = MD.timeInterval_/60; % in min
    close all
    numG7=numel(tracksG7);
    risingAndStoppingG7=[];
    edgeRisingStoppingG7{1}=[];
    adhAdvanceRisingStoppingG7{1}=[];
    divPointRS=[];
    tractionRisingStoppingG7{1}=[];
    risingSlopeG7=[];
    forceIncrease1Half=[];
    forceIncrease2Half=[];
    edgeIncrease1Half=[];
    edgeIncrease2Half=[];
    adhIncrease1Half=[];
    adhIncrease2Half=[];
    p=0;
    q=0;
    for kk=1:numG7
        curG7=tracksG7(kk);
        if curG7.lifeTime<5
            continue
        end
        tRange = 1:length(curG7.edgeAdvanceDist);
        if curG7.startingFrameExtra>1
            curG7.edgeAdvanceDist(1:curG7.startingFrameExtra-1)=NaN;
            curG7.forceMag(1:curG7.startingFrameExtra-1)=NaN;
            curG7.advanceDist(1:curG7.startingFrameExtra-1)=NaN;
        end
        % fit with partial sloped line + flat line in the later
        x=(0:(sum(~isnan(curG7.edgeAdvanceDist))-1))*timeInterval;
        y=curG7.edgeAdvanceDist(~isnan(curG7.edgeAdvanceDist));
        fun=@(z) [(z(1)*x(x<=z(3))+z(2)) (z(1)*z(3)+z(2))*ones(1,sum(x>z(3)))]-y;
        z0=[0.1,0,20*timeInterval];
    %     [x0out,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0i);
        z=lsqnonlin(fun,z0);
        % Since the purpose is to separate out those that increase about a
        % half and stay the rest of half, I'll filter out a<=0,
        % xt/lifetime>0.6.
        if z(1)>0.1 && z(3)/x(end)<0.8 && z(3)/x(end)>0.2
            p=p+1;
            q=q+1;
            risingAndStoppingG7=[risingAndStoppingG7 kk];
            curG7Advance=curG7.advanceDist(~isnan(curG7.advanceDist));
            figure, subplot(2,1,1),plot(x,curG7Advance,'o--','Color',[1,140/255,0]), hold on,
            plot(x,y,'ko--'), 
            plot(x,[(z(1)*x(x<=z(3))+z(2)) (z(1)*z(3)+z(2))*ones(1,sum(x>z(3)))],'k-','Linewidth',2)
            title(['Edge Protrusion, ID:' num2str(kk)]); ylabel('Protrusion distance (um)')
            xt = round(z(3)/timeInterval);
            edgeRisingStoppingG7{p}=y;
            adhAdvanceRisingStoppingG7{p}=curG7Advance;
            divPointRS=[divPointRS xt];
            risingSlopeG7=[risingSlopeG7 z(1)];

            try
                sForce= csaps(tRange,curG7.forceMag,splineParam);
            catch
                curG7.forceMag(curG7.endingFrameExtra+1:end)=[];
                sForce= csaps(tRange,curG7.forceMag,splineParam);
            end
            sdForce=ppval(sForce,tRange);
            sdForce(isnan(curG7.forceMag))=NaN;
            curG7Force = sdForce(~isnan(curG7.edgeAdvanceDist));
            curG7ForceRaw = curG7.forceMag(~isnan(curG7.edgeAdvanceDist));
            tractionRisingStoppingG7{p}=curG7ForceRaw;
        %     [~,m_first]=regression(x(1:xt),y(1:xt))
        %     [~,m_second]=regression(x(xt:end),y(xt:end))
        %     [~,m_first]=regression(x(1:xt),curG7Force(1:xt))
        %     [~,m_second]=regression(x(xt:end),curG7Force(xt:end))
            forceIncrease1Half(p)=(max(curG7ForceRaw(xt-3:xt))-min(curG7ForceRaw(1:3)))/z(3);
            forceIncrease2Half(q)=(max(curG7ForceRaw(end-5:end))-min(curG7ForceRaw(xt-1:xt+1)))/(x(end)-z(3));
            edgeIncrease1Half(p)=(max(y(xt-3:xt))-min(y(1:3)))/z(3);
            edgeIncrease2Half(q)=(max(y(end-5:end))-min(y(xt-1:xt+1)))/(x(end)-z(3));
            adhIncrease1Half(p)=(max(curG7Advance(xt-3:xt))-min(curG7Advance(1:3)))/z(3);
            adhIncrease2Half(q)=(max(curG7Advance(end-5:end))-min(curG7Advance(xt-1:xt+1)))/(x(end)-z(3));
            subplot(2,1,2),plot(x,curG7ForceRaw,'ko-'),hold on,plot(x,curG7Force,'k-','LineWidth',2)
            title('Traction'); ylabel('Traction Magnitude (Pa)'); xlabel('Time (min)');
            hgsave(strcat(figPath,'/EdgeAndForceInG7',num2str(kk)),'-v7.3')
            print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'EdgeAndForceInG7' num2str(kk) '.eps']);% histogramPeakLagVinVsTal -transparent
        else
            q=q+1;
            curG7ForceRaw = curG7.forceMag(~isnan(curG7.edgeAdvanceDist));
            forceIncrease2Half(q)=(max(curG7ForceRaw(end-5:end))-min(curG7ForceRaw(1:3)))/(x(end)-x(1));
        end
    end
    save(strcat(dataPath,'/EdgeAndForceInG7selected.mat'),'risingAndStoppingG7','risingSlopeG7',...
        'forceIncrease1Half','forceIncrease2Half','edgeIncrease1Half','edgeIncrease2Half','adhIncrease1Half','adhIncrease2Half',...
        'edgeRisingStoppingG7','adhAdvanceRisingStoppingG7','divPointRS','tractionRisingStoppingG7','-v7.3')
    %% Filtering G3 for those that increase edge advance 
    close all
    numG3=numel(tracksG3);
    risingAndStoppingG3=[];
    edgeRisingStoppingG3{1}=[];
    adhAdvanceRisingStoppingG3{1}=[];
    forceIncreaseG3=[];
    edgeIncreaseG3=[];
    adhIncreaseG3=[];
    tractionRisingStoppingG3{1}=[];
    risingSlopeG3=[];
    p=0;
    for kk=1:numG3
        curG3=tracksG3(kk);
        if curG3.lifeTime<5
            continue
        end
        tRange = 1:length(curG3.edgeAdvanceDist);
        if curG3.startingFrameExtra>1
            curG3.edgeAdvanceDist(1:curG3.startingFrameExtra-1)=NaN;
            curG3.forceMag(1:curG3.startingFrameExtra-1)=NaN;
            curG3.advanceDist(1:curG3.startingFrameExtra-1)=NaN;
        end
        % fit with partial sloped line + flat line in the later
        x=(0:(sum(~isnan(curG3.edgeAdvanceDist))-1))*timeInterval;
        y=curG3.edgeAdvanceDist(~isnan(curG3.edgeAdvanceDist));
        fun=@(z) (z(1)*x+z(2))-y;
        z0=[0.1,0];
    %     [x0out,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0i);
        z=lsqnonlin(fun,z0);
        % Since the purpose is to separate out those that increase about a
        % half and stay the rest of half, I'll filter out a<=0,
        % xt/lifetime>0.6.
        if z(1)>0.1
            p=p+1;
            risingAndStoppingG3=[risingAndStoppingG3 kk];
            curG3Advance=curG3.advanceDist(~isnan(curG3.advanceDist));
            figure, subplot(2,1,1),plot(x,curG3Advance,'o--','Color',[1,140/255,0]), hold on,
            plot(x,y,'ko--'), 
            plot(x,(z(1)*x+z(2)),'k-','Linewidth',2)
            title(['Edge Protrusion, ID:' num2str(kk)]); ylabel('Protrusion distance (um)')
            edgeRisingStoppingG3{p}=y;
            adhAdvanceRisingStoppingG3{p}=curG3Advance;
            risingSlopeG3=[risingSlopeG3 z(1)];

            try
                sForce= csaps(tRange,curG3.forceMag,splineParam);
            catch
                curG3.forceMag(curG3.endingFrameExtra+1:end)=[];
                sForce= csaps(tRange,curG3.forceMag,splineParam);
            end
            sdForce=ppval(sForce,tRange);
            sdForce(isnan(curG3.forceMag))=NaN;
            curG3Force = sdForce(~isnan(curG3.edgeAdvanceDist));
            curG3ForceRaw = curG3.forceMag(~isnan(curG3.edgeAdvanceDist));
            tractionRisingStoppingG3{p}=curG3ForceRaw;

            forceIncreaseG3(p)=(max(curG3Force(end-2:end))-min(curG3Force(1:3)))/x(end);
            edgeIncreaseG3(p)=(max(y(end-2:end))-min(y(1:3)))/x(end);
            adhIncreaseG3(p)=(max(curG3Advance(end-2:end))-min(curG3Advance(1:3)))/x(end);
            subplot(2,1,2),plot(x,curG3ForceRaw,'ko-'),hold on,plot(x,curG3Force,'k-','LineWidth',2)
            title('Traction'); ylabel('Traction Magnitude (Pa)'); xlabel('Time (min)');
            hgsave(strcat(figPath,'/EdgeAndForceInG3',num2str(kk)),'-v7.3')
            print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'EdgeAndForceInG3' num2str(kk) '.eps']);% histogramPeakLagVinVsTal -transparent
        end
    end
    save(strcat(dataPath,'/EdgeAndForceInG3selected.mat'),'risingAndStoppingG3','edgeRisingStoppingG3',...
        'adhAdvanceRisingStoppingG3','forceIncreaseG3','edgeIncreaseG3','adhIncreaseG3','tractionRisingStoppingG3','-v7.3')
    %% Last, compare force increase among forceIncreaseG3,forceIncrease1Half and forceIncrease2Half
    close all
    forceIncreaseG3G7Cell={ forceIncreaseG3,forceIncrease1Half, forceIncrease2Half};
    nameG3G7={'G3', 'G7-protruding phase', 'G7-stalling phase'};
    figure;
    barPlotCellArray(forceIncreaseG3G7Cell,nameG3G7)
    title('Increase in force in G3 and two separate phases in G7')
    ylabel('Change in force (Pa/min)')
    hgsave(strcat(figPath,'/forceIncreaseG3G7'),'-v7.3')
    save([pathForColocalization filesep 'data' filesep 'forceIncreaseG3G7.mat'],'forceIncreaseG3G7Cell','nameG3G7','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'forceIncreaseG3G7.eps']);% histogramPeakLagVinVsTal -transparent

    edgeAdhIncreaseG3G7Cell={ edgeIncreaseG3,adhIncreaseG3,edgeIncrease1Half,adhIncrease1Half,...
        edgeIncrease2Half,adhIncrease2Half};
    nameG3G7edgeAdh={'G3-edge','Ge-adh', 'G7-edge-protruding phase','G7-adh-protruding phase',...
        'G7-edge-stalling phase','G7-adh-stalling phase'};
    figure;
    barPlotCellArray(edgeAdhIncreaseG3G7Cell,nameG3G7edgeAdh,pixSize)
    title('Advance in edge and adhesion in G3 and two separate phases in G7')
    ylabel('Change in edge (um/min)')
    hgsave(strcat(figPath,'/edgeAdhIncreaseG3G7'),'-v7.3')
    save([pathForColocalization filesep 'data' filesep 'edgeAdhIncreaseG3G7.mat'],'edgeAdhIncreaseG3G7Cell','nameG3G7edgeAdh','-v7.3')
    print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'edgeAdhIncreaseG3G7.eps']);% histogramPeakLagVinVsTal -transparent

    
    %% Look at G7 in detail and see force characteristics
    % First let's look at edge advance
%     curG7ID=0;
    %% individual G7 cells: [5 10 13 21 24 28 71 91 113 118 145] 
%     close
%     curG7ID=curG7ID+1;
%     curG7=tracksG7(curG7ID);
%     if curG7.startingFrameExtra>1
%         curG7.edgeAdvanceDist(1:curG7.startingFrameExtra-1)=NaN;
%     end
%     figure, plot(curG7.edgeAdvanceDist), title(num2str(curG7ID))
    %% Look at force in selected G7
%     idsG7=[5 10 13 21 24 28 71 91 113 118 145]; 
%     curG7ID=0;
    %% Actual looking
%     curG7ID=curG7ID+1;
%     curG7=tracksG7(idsG7(curG7ID));
%     if curG7.startingFrameExtra>1
%         curG7.edgeAdvanceDist(1:curG7.startingFrameExtra-1)=NaN;
%         curG7.forceMag(1:curG7.startingFrameExtra-1)=NaN;
%     end
%     figure, subplot(2,1,1), plot(curG7.edgeAdvanceDist), title(num2str(curG7ID)), ylabel('Edge advance (um)')
%     tRange = 1:length(curG7.edgeAdvanceDist);
%     sEdgeAdv= csaps(tRange,curG7.edgeAdvanceDist,splineParam);
%     sdEdgeAdv=ppval(sEdgeAdv,tRange);
%     sdEdgeAdv(isnan(curG7.edgeAdvanceDist))=NaN;
%     hold on, plot(sdEdgeAdv)
%     
%     subplot(2,1,2), plot(curG7.forceMag), ylabel('Force (Pa)')
%     sForce= csaps(tRange,curG7.forceMag,splineParam);
%     sdForce=ppval(sForce,tRange);
%     sdForce(isnan(curG7.forceMag))=NaN;
%     hold on, plot(sdForce)
% 
%     hgsave(strcat(figPath,'/EdgeAndForce',num2str(curG7ID)),'-v7.3')
%     print('-depsc','-loose',[pathForColocalization filesep 'eps' filesep 'EdgeAndForce' num2str(curG7ID) '.eps']);% histogramPeakLagVinVsTal -transparent
    %% Now I'll smartly separate G7 (G3) by least-square difference
    % z=[a,b,xt]
%     x=1:sum(~isnan(curG7.edgeAdvanceDist));
%     y=curG7.edgeAdvanceDist(~isnan(curG7.edgeAdvanceDist));
%     fun=@(z) [(z(1)*x(x<=z(3))+z(2)) (z(1)*z(3)+z(2))*ones(1,sum(x>z(3)))]-y;
%     x0i=[0.1,0,20];
% %     [x0out,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0i);
%     x0out=lsqnonlin(fun,x0i);
%     figure, plot(x,y,'ko',x,[(x0out(1)*x(x<=x0out(3))+x0out(2)) (x0out(1)*x0out(3)+x0out(2))*ones(1,sum(x>x0out(3)))],'b-')
%     % Force increase during both periods
%     %First phase
%     xt=round(x0out(3));
% %     curG7Force = curG7.forceMag(~isnan(curG7.forceMag));
%     curG7Force = sdForce(~isnan(curG7.forceMag));
%     
% %     [~,m_first]=regression(x(1:xt),y(1:xt))
% %     [~,m_second]=regression(x(xt:end),y(xt:end))
% %     [~,m_first]=regression(x(1:xt),curG7Force(1:xt))
% %     [~,m_second]=regression(x(xt:end),curG7Force(xt:end))
%     forceIncrease1Half=max(curG7Force(xt-3:xt))-min(curG7Force(1:3));
%     forceIncrease2Half=max(curG7Force(end-5:end))-min(curG7Force(xt-1:xt+1));
    %% Backup original tracksNA and save filtered one
    movefile([pathForColocalization filesep 'data' filesep 'tracksNA.mat'], [pathForColocalization filesep 'data' filesep 'tracksNA_org.mat'])
    save([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA','-v7.3')
    %% Backup original tracksNA and save filtered one
    movefile([pathForColocalization filesep 'data' filesep 'idsClassified.mat'], [pathForColocalization filesep 'data' filesep 'idsClassified_org.mat'])
    save([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9',...
        'idGroup1filtered','idGroup2filtered','idGroup3filtered','idGroup4filtered','idGroup5filtered','idGroup6filtered','idGroup7filtered','idGroup8filtered','idGroup9filtered','tracksNA','-v7.3')
    %% Save all workspace
    save([pathForColocalization filesep 'data' filesep 'allDataAfterClassification.mat'],'-v7.3')
    %% Make movie
    close all
%     makeMoviesAdhesionTracks(tracksNA,pathForColocalization,1500,100)
%     makeMoviesAdhesionTracks(tracksNA,pathForColocalization,400,100) % for 0724paxillin1
%     makeMoviesAdhesionTracks(tracksNA,pathForColocalization,2000,100) % for 0128Fibro3t3pax
%     makeMoviesAdhesionTracks(tracksNA,pathForColocalization,200,100) % for 0724paxillin2
%     makeMoviesAdhesionTracks(tracksNA,pathForColocalization,200,100) % for 0724talin
    tmaxG1 = max(avgPreDetecForceMag_forceTransmitting);
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,10) % automatic
    %% Make movie - without classification
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,0) % automatic
    %% Make movie - only for class 1
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,1) % automatic
    %% Make movie - only for class 2
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,2) % automatic
    %% Make movie - only for class 3
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,3) % automatic
    %% Make movie - only for class 7
    makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmaxG1,16,7) % automatic
end
end
