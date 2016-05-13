function [tracksAD] = analyzeNAMaturationAndFADecay(pathForTheMovieDataFile,showAllTracks,varargin)
% [tracksAD,indNAFail,indNAMature,indFAgrowing,indFAdecaying] = analyzeNAMaturationAndFADecay(pathForTheMovieDataFile,outputPath,showAllTracks,plotEachTrack)
% performs point sources vs. segmented FA colocalization analysis and
% determine growing NAs, failing NAs, growing FAs and decaying FAs based on
% the size.

% input:    pathForTheMovieDataFile:    path to the movieData file (FA
%                   segmentation and NA tracking package should be run beforehand)
%           outputPath                  outputPath
%           showAllTracks             true if you want to show all tracks
%                                               for whole cell
%           plotEachTrack              true if you want to plot each individual track
%                                               with small roi
% output:   images will be stored in pathForTheMovieDataFile/trackFrames
%           tracksNAfailing,          tracks of failing NAs that go on to turn-over
%           tracksNAmaturing,          tracks of failing NAs that matures to FAs
%           lifeTimeNAfailing,            lifetime of all NAs that turn over 
%           lifeTimeNAmaturing,            lifetime of all maturing NAs until their final turn-over
%           maturingRatio,            ratio of maturing NAs w.r.t. all NA tracks 
%           NADensity                   density of nascent adhesions, unit: number/um2
%           FADensity                   density of focal adhesions , unit: number/um2

% status of each track
%           BA,          Before Adhesion
%           NA,          Nascent Adhesion
%           FC,            Focal Contact
%           FA,            Focal Adhesion
%           ANA,           After Nascent Adhesion
%           Out_of_Band,            Out of band from the cell edge

% Sangyoon Han Feb 2016

%% Inputs
ip =inputParser;
ip.addRequired('pathForTheMovieDataFile',@ischar)
ip.addOptional('showAllTracks',false,@(x)islogical(x)||isempty(x))
ip.addOptional('plotEachTrack',false,@(x)islogical(x)||isempty(x))
ip.addParamValue('onlyEdge',false,@islogical); % collect NA tracks that ever close to cell edge
ip.addParamValue('outputPath','AdhesionMaturationAndDecaying',@ischar)
ip.addParamValue('saveAnalysis',true,@islogical)
ip.addParamValue('matchWithFA',true,@islogical) %For cells with only NAs, we turn this off.
ip.addParamValue('minLifetime',5,@isscalar) %For cells with only NAs, we turn this off.
ip.addParamValue('skipImages',1,@isscalar) %For cells with only NAs, we turn this off.

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
ip.parse(pathForTheMovieDataFile,showAllTracks,varargin{:});
outputPath=ip.Results.outputPath;
saveAnalysis=ip.Results.saveAnalysis;
matchWithFA=ip.Results.matchWithFA;
minLifetime=ip.Results.minLifetime;
% showAllTracks=ip.Results.showAllTracks;
% plotEachTrack=ip.Results.plotEachTrack;
onlyEdge=ip.Results.onlyEdge;
skipImages = ip.Results.skipImages;
% Load the Paxillin channel

%% Data Set up
% Load the MovieData
% movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
[pathStr,name,ext]= fileparts(pathForTheMovieDataFile);
if strcmp([name ext],'movieData.mat')
    movieDataPath = pathForTheMovieDataFile;
    pathForTheMovieDataFile = pathStr;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
end
MD = MovieData.load(movieDataPath);
% MD = MovieData.load(pathForTheMovieDataFile);
% Get whole frame number
nFrames = MD.nFrames_;
% Get U-Track package (actually NA package)
try
    NAPackage = MD.getPackage(MD.getPackageIndex('UTrackPackage'));
    % Load tracks
    iTracking = 2;
    trackNAProc = NAPackage.processes_{iTracking};
    detectedNAProc = NAPackage.processes_{1};
catch % this is to accoutn for the analysis that is done in the old version of colocalizationAdhesionWithTFM
    trackNAProc = MD.getProcess(MD.getProcessIndex('TrackingProcess'));
    detectedNAProc = MD.getProcess(MD.getProcessIndex('AnisoGaussianDetectionProcess'));
end
% Get FA segmentation package 
FASegPackage = MD.getPackage(MD.getPackageIndex('FocalAdhesionSegmentationPackage'));

% Set up the output file path
outputFilePath = [pathForTheMovieDataFile filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
paxPath = [outputFilePath filesep 'pax'];
paxtifPath = [outputFilePath filesep 'paxtifs'];
epsPath = [outputFilePath filesep 'eps'];
figPath = [outputFilePath filesep 'figs'];
if ~exist(paxtifPath,'dir') || ~exist(paxPath,'dir') || ~exist(figPath,'dir') || ~exist(epsPath,'dir') || ~exist(dataPath,'dir') 
    mkdir(paxPath);
    mkdir(paxtifPath);
    mkdir(figPath);
    mkdir(epsPath);
    mkdir(dataPath);
end
iiformat = ['%.' '3' 'd'];
%     paxLevel = zeros(nFrames,1);
% SegmentationPackage = MD.getPackage(MD.getPackageIndex('SegmentationPackage'));
% minSize = round((500/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=.5um*.5um
minLifetime = min(nFrames,minLifetime);
markerSize = 4;
% tracks
iPaxChannel = 1; % this should be intentionally done in the analysis level
if ~trackNAProc.checkChannelOutput(iPaxChannel)
    iPaxChannel = 2;
end
% filter out tracks that have lifetime less than 2 frames
disp('loading NA tracks...')
tic
tracksNAorg = trackNAProc.loadChannelOutput(iPaxChannel);
toc
% if ii>minLifetime
SEL = getTrackSEL(tracksNAorg); %SEL: StartEndLifetime
% Remove any less than 3-frame long track.
isValid = SEL(:,3) >= minLifetime;
tracksNAorg = tracksNAorg(isValid);
% end
detectedNAs = detectedNAProc.loadChannelOutput(iPaxChannel);

% re-express tracksNA so that each track has information for every frame
disp('reformating NA tracks...')
tic
tracksAD = formatTracks(tracksNAorg,detectedNAs,nFrames); 
toc

% disp('loading segmented FAs...')
iFASeg = 6;
FASegProc = FASegPackage.processes_{iFASeg};
bandArea = zeros(nFrames,1);
NADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
FADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
numFAs = zeros(nFrames,1);
nChannels = numel(MD.channels_);
iChan = 0;
% Finding which channel has a cell mask information
maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
for k=1:nChannels
    if maskProc.checkChannelOutput(k)
        iChan = k;
    end
end
% Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered
if onlyEdge
    disp(['Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered'])
    trackIdx = false(numel(tracksAD),1);
    bandwidthNA = 7; %um 
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    for ii=1:nFrames
        % Cell Boundary Mask 
        mask = maskProc.loadChannelOutput(iChan,ii);
        % mask for band from edge
        iMask = imcomplement(mask);
        distFromEdge = bwdist(iMask);
        bandMask = distFromEdge <= bandwidthNA_pix;

        maskOnlyBand = bandMask & mask;
        bandArea(ii) = sum(maskOnlyBand(:)); % in pixel
        
        % collect index of tracks which first appear at frame ii
        idxFirstAppear = arrayfun(@(x) x.startingFrame==ii,tracksAD);
        % now see if these tracks ever in the maskOnlyBand
        for k=find(idxFirstAppear)'
            if maskOnlyBand(round(tracksAD(k).yCoord(ii)),round(tracksAD(k).xCoord(ii)))
                trackIdx(k) = true;
            end
        end
    end    
else
    disp('Entire adhesion tracks are considered.')
    trackIdx = true(numel(tracksAD),1);
    mask = maskProc.loadChannelOutput(iChan,1);
end
% get rid of tracks that have out of bands...
tracksAD = tracksAD(trackIdx);
%% add intensity of tracks including before and after NA status
% get the movie stack
% [h,w]=size(mask);
% imgStack = zeros(h,w,nFrames);
% for ii=1:nFrames
%     imgStack(:,:,ii)=MD.channels_(iAdh).loadImage(ii); 
% end
% get the intensity
% disp('Reading intensities without additional tracking...')
% tic
% tracksAD = readIntensityFromTracks(tracksAD,imgStack,1,'reTrack',false); % 1 means intensity collection from pax image
% toc
firstMask=maskProc.loadChannelOutput(iChan,1);
cropMaskStack = false(size(firstMask,1),size(firstMask,2),nFrames);
numTracks=numel(tracksAD);
progressText(0,'Matching with segmented adhesions:');
%% Matching with segmented adhesions
% skipImages=6; %This is for Tristan's.
for ii=1:nFrames
    % Cell Boundary Mask 
    mask = maskProc.loadChannelOutput(iChan,ii);
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');
    cropMaskStack(:,:,ii) = maskProc.loadChannelOutput(iChan,ii);
    if matchWithFA
        % filter tracks with naMasks
%         disp(['Processing ' num2str(ii) 'th frame out of ' num2str(nFrames) ' total frames...'])
        % focal contact (FC) analysis
        % Get the mask for FAs
        maskFAs = FASegProc.loadChannelOutput(iPaxChannel,ii);
        maskAdhesion = maskFAs>0;
        Adhs = regionprops(maskAdhesion,'Area','Eccentricity','PixelIdxList','PixelList' );
    %     propFAs = regionprops(maskFAs,'Area','Eccentricity','PixelIdxList','PixelList' );
        numFAs(ii) = numel(Adhs);
        minFASize = round((1000/MD.pixelSize_)*(1000/MD.pixelSize_)); %adhesion limit=1um*1um
        minFCSize = round((600/MD.pixelSize_)*(400/MD.pixelSize_)); %adhesion limit=0.6um*0.4um

        fcIdx = arrayfun(@(x) x.Area<minFASize & x.Area>minFCSize, Adhs);
        FCIdx = find(fcIdx);
        adhBound = bwboundaries(maskAdhesion,'noholes');    

        % for larger adhesions
        faIdx = arrayfun(@(x) x.Area>=minFASize, Adhs);
        FAIdx =  find(faIdx);
        neighPix = 2;


        % Deciding each adhesion maturation status
        for k=1:numel(tracksAD)
            if tracksAD(k).presence(ii)
                if ~strcmp(tracksAD(k).state{ii} , 'NA') && ii>1
                    tracksAD(k).state{ii} = tracksAD(k).state{ii-1};
                end
                % decide if each track is associated with FC or FA
                p = 0;
                for jj=FCIdx'
                    p=p+1;
                    if any(round(tracksAD(k).xCoord(ii))==Adhs(jj).PixelList(:,1) & round(tracksAD(k).yCoord(ii))==Adhs(jj).PixelList(:,2))
                        tracksAD(k).state{ii} = 'FC';
                        tracksAD(k).area(ii) = Adhs(jj).Area;% in pixel
                        tracksAD(k).FApixelList{ii} = Adhs(jj).PixelList;
                        tracksAD(k).adhBoundary{ii} = adhBound{jj};
                        tracksAD(k).faID(ii) = maskFAs(round(tracksAD(k).yCoord(ii)),round(tracksAD(k).xCoord(ii)));
                    end
                end
                p = 0;
                for jj=FAIdx'
                    p=p+1;
                    if any(round(tracksAD(k).xCoord(ii))==Adhs(jj).PixelList(:,1) & round(tracksAD(k).yCoord(ii))==Adhs(jj).PixelList(:,2))
                        tracksAD(k).state{ii} = 'FA';
                        tracksAD(k).area(ii) = Adhs(jj).Area;% in pixel
                        tracksAD(k).FApixelList{ii} = Adhs(jj).PixelList;
                        tracksAD(k).adhBoundary{ii} = adhBound{jj};
                        tracksAD(k).faID(ii) = maskFAs(round(tracksAD(k).yCoord(ii)),round(tracksAD(k).xCoord(ii)));
                    end
                end
            elseif ii>tracksAD(k).endingFrame && (strcmp(tracksAD(k).state{tracksAD(k).endingFrame},'FA')...
                    || strcmp(tracksAD(k).state{tracksAD(k).endingFrame},'FC'))
                % starting from indexed maskFAs, find out segmentation that is
                % closest to the last track point.
                subMaskFAs = maskFAs==tracksAD(k).faID(tracksAD(k).endingFrame);
                if max(subMaskFAs(:))==0
                    tracksAD(k).state{ii} = 'ANA';
                    tracksAD(k).FApixelList{ii} = NaN;
                    tracksAD(k).adhBoundary{ii} = NaN;
                    continue
                else
                    propSubMaskFAs = regionprops(subMaskFAs,'PixelList');
                    minDist = zeros(length(propSubMaskFAs),1);
                    for q=1:length(propSubMaskFAs)
                        minDist(q) = min(sqrt((propSubMaskFAs(q).PixelList(:,1)-(tracksAD(k).xCoord(tracksAD(k).endingFrame))).^2 +...
                            (propSubMaskFAs(q).PixelList(:,2)-(tracksAD(k).yCoord(tracksAD(k).endingFrame))).^2));
                    end
                    % find the closest segment
                    [~,subMaskFAsIdx] = min(minDist);
                    subAdhBound = bwboundaries(subMaskFAs,'noholes');    
                    [~,closestPixelID] = min(sqrt((propSubMaskFAs(subMaskFAsIdx).PixelList(:,1)-(tracksAD(k).xCoord(tracksAD(k).endingFrame))).^2 +...
                        (propSubMaskFAs(subMaskFAsIdx).PixelList(:,2)-(tracksAD(k).yCoord(tracksAD(k).endingFrame))).^2));

                    tracksAD(k).state{ii} = 'FC';
%                     tracksNA(k).xCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,1);
%                     tracksNA(k).yCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,2);
                    tracksAD(k).FApixelList{ii} = propSubMaskFAs(subMaskFAsIdx).PixelList;
                    tracksAD(k).adhBoundary{ii} = subAdhBound{subMaskFAsIdx};
                end
            end
        end
    else
        FCIdx = [];
        FAIdx = [];
    end
    % recording features
    % get the point on the boundary closest to the adhesion
    allBdPoints = [];
    for kk=1:nBD
        boundary = B{kk};
        allBdPoints = [allBdPoints; boundary(:,2), boundary(:,1)];
    end
    for k=1:numel(tracksAD)
        % distance to the cell edge
%         if tracksNA(k).presence(ii)
        if ii>=tracksAD(k).startingFrame && ii<=tracksAD(k).endingFrame
            xCropped = tracksAD(k).xCoord(ii);
            yCropped = tracksAD(k).yCoord(ii);
            distToAdh = sqrt(sum((allBdPoints- ...
                ones(size(allBdPoints,1),1)*[xCropped, yCropped]).^2,2));
            [minDistToBd,indMinBdPoint] = min(distToAdh);
            tracksAD(k).distToEdge(ii) = minDistToBd;
            tracksAD(k).closestBdPoint(ii,:) = allBdPoints(indMinBdPoint,:); % this is lab frame of reference. (not relative to adhesion position)
        end
    end
    
    if showAllTracks && ismember(ii,1:skipImages:nFrames)
        if ii==1
            h2=figure;
        else
            figure(h2)
        end
        %Scale bar 2 um
    %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
        paxImageCropped=MD.channels_(iPaxChannel).loadImage(ii); 
        % size of the region of interest
        imshow(paxImageCropped,[]), hold on
        line([10 10+round(2000/MD.pixelSize_)],[15 15],'LineWidth',2,'Color','w')
        
        for kk=1:nBD
            boundary = B{kk};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        end
    %     plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
        for k = FCIdx'
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
        end
        % for larger adhesions
        for k = FAIdx'
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
        end
        for k=1:numel(tracksAD)
            if tracksAD(k).presence(ii)
                if strcmp(tracksAD(k).state{ii} , 'NA')
                    % drawing tracks
                    plot(tracksAD(k).xCoord(1:ii),tracksAD(k).yCoord(1:ii),'r', 'LineWidth', 0.5)
                    plot(tracksAD(k).xCoord(ii),tracksAD(k).yCoord(ii),'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksAD(k).state{ii} , 'FC')
                    % drawing tracks
                    plot(tracksAD(k).xCoord(1:ii),tracksAD(k).yCoord(1:ii),'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                    plot(tracksAD(k).xCoord(ii),tracksAD(k).yCoord(ii),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksAD(k).state{ii} , 'FA')
                    % drawing tracks
                    plot(tracksAD(k).xCoord(1:ii),tracksAD(k).yCoord(1:ii),'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                    plot(tracksAD(k).xCoord(ii),tracksAD(k).yCoord(ii),'Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                end
            end
        end

        print(h2, '-depsc2', strcat(epsPath,'/pax',num2str(ii,iiformat),'.eps'));
        print(h2, '-dtiff', strcat(paxtifPath,'/pax',num2str(ii,iiformat),'.tif'));
    %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
        hold off
    end
%     if saveAnalysis
%         imwrite(paxImageCropped,strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
%     end
    progressText(ii/(nFrames-1),'Matching with segmented adhesions:');
end
%% disp('Intermediate saving before post analysis...')
disp('Intermediate saving before post analysis...')
save([dataPath filesep 'tracksAD_beforePostAnalysis.mat'], 'tracksAD','-v7.3')
% save([dataPath filesep 'intermediateWorkspace.mat'], '-v7.3')
%% protrusion/retraction information
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.
for ii=1:numTracks
    idxZeros = tracksAD(ii).closestBdPoint(:)==0;
    tracksAD(ii).closestBdPoint(idxZeros)=NaN;
end

disp('Post-analysis on adhesion movement...')
deltaT = MD.timeInterval_; % sampling rate (in seconds, every deltaT seconds)
tIntervalMin = deltaT/60; % in min
periodMin = 1;
periodFrames = floor(periodMin/tIntervalMin); % early period in frames
frames2min = floor(2/tIntervalMin); % early period in frames
progressText(0,'Post-analysis:');
for k=1:numTracks
    % cross-correlation scores
    presIdx = logical(tracksAD(k).presence);
    % get the instantaneous velocity
    % get the distance first
    curTrackVelLength=sum(presIdx)-1;
    distTrajec=zeros(curTrackVelLength,1);
    presIdxSeq = find(presIdx);
    for kk=1:curTrackVelLength
        real_kk = presIdxSeq(kk);
        distTrajec(kk) = sqrt(sum((tracksAD(k).closestBdPoint(real_kk+1,:)- ...
            tracksAD(k).closestBdPoint(real_kk,:)).^2,2));
        lastPointIntX = round(tracksAD(k).closestBdPoint(real_kk+1,1));
        lastPointIntY = round(tracksAD(k).closestBdPoint(real_kk+1,2));
        if cropMaskStack(lastPointIntY,lastPointIntX,real_kk) %if the last point is in the first mask, it is inward
            distTrajec(kk) = -distTrajec(kk);
        end
    end
    if any(distTrajec~=0)
        [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
        if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
            tracksAD(k).isProtrusion = true;
        else
            tracksAD(k).isProtrusion = false;
        end
        % average velocity (positive for protrusion)
        curProtVel = (Protrusion.Veloc); curProtVel(isnan(curProtVel))=0;
        curProtPersTime = (Protrusion.persTime); curProtPersTime(isnan(curProtPersTime))=0;
        curRetVel = (Retraction.Veloc); curRetVel(isnan(curRetVel))=0;
        curRetPersTime = (Retraction.persTime); curRetPersTime(isnan(curRetPersTime))=0;

        tracksAD(k).edgeVel = (mean(curProtVel.*curProtPersTime)-mean(curRetVel.*curRetPersTime))/mean([curProtPersTime;curRetPersTime]);
    else
        tracksAD(k).edgeVel = 0;
    end
    % lifetime information
    sF=tracksAD(k).startingFrame;
    eF=tracksAD(k).endingFrame;
    tracksAD(k).lifeTime = eF-sF;    
    % Inital intensity slope for one min
    timeInterval = deltaT/60; % in min
    earlyPeriod = floor(1/timeInterval); % frames per minute
    lastFrame = min(sum(~isnan(tracksAD(k).amp)),sF+earlyPeriod-1);
    lastFrameFromOne = lastFrame - sF+1;
%     lastFrameFromOne = sF;
%     lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
    [curR,curM] = regression(timeInterval*(1:lastFrameFromOne),tracksAD(k).amp(sF:lastFrame));
    tracksAD(k).ampSlope = curM; % in a.u./min
    tracksAD(k).ampSlopeR = curR; % Pearson's correlation coefficient
    
    curEndFrame = min(sF+periodFrames-1,eF);
    curEarlyPeriod = curEndFrame - sF+1;
    [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),tracksAD(k).amp(sF:curEndFrame));
    tracksAD(k).earlyAmpSlope = curM; % in a.u./min

    curStartFrame = max(tracksAD(k).startingFrame,eF-periodFrames+1);
    curLatePeriod = eF - curStartFrame+1;
    [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksAD(k).amp(curStartFrame:eF));
    tracksAD(k).lateAmpSlope = curMlate; % in a.u./min

    curEndFrame = min(sF+periodFrames-1,eF);
    curEarlyPeriod = curEndFrame - sF+1;
    [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),tracksAD(k).distToEdge(sF:curEndFrame));
    tracksAD(k).distToEdgeSlope = curMdist; % in a.u./min
    tracksAD(k).distToEdgeChange = (tracksAD(k).distToEdge(end)-tracksAD(k).distToEdge(tracksAD(k).startingFrame)); % in pixel
    % Determining protrusion/retraction based on closestBdPoint and [xCoord
    % yCoord]. If the edge at one time point does not cross the adhesion
    % tracks over entire life time, there is no problem. But that's not
    % always the case: adhesion tracks crosses the cell boundary at first
    % or last time point. And we don't know if the adhesion tracks are
    % moving in the direction of protrusion or retraction. Thus, what I
    % will do is to calculate the inner product of vectors from the
    % adhesion tracks from the closestBdPoint at the first frame. If the
    % product is positive, it means both adhesion points are in the same
    % side. And if distance from the boundary to the last point is larger
    % than the distance to the first point, it means the adhesion track is
    % retracting. In the opposite case, the track is protruding (but it
    % will happen less likely because usually the track would cross the
    % first frame boundary if it is in the protrusion direction). 
    
    % We need to find out boundary points projected on the line of adhesion track 
    fitobj = fit(tracksAD(k).xCoord(sF:eF)',tracksAD(k).yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
    x0=nanmedian(tracksAD(k).xCoord);
    y0=fitobj(x0);
    dx = 1;
    dy = fitobj.p1;
    trackLine = createLineGeom2d(x0,y0,dx,dy); % this is a geometric average linear line of the adhesion track
    
%     trackLine = edgeToLine(edge);
    firstBdPoint = [tracksAD(k).closestBdPoint(sF,1) tracksAD(k).closestBdPoint(sF,2)];
    firstBdPointProjected = projPointOnLine(firstBdPoint, trackLine); % this is an edge boundary point at the first time point projected on the average line of track.
    % try to record advanceDist and edgeAdvanceDist for every single time
    % point ...
    for ii=sF:eF
        curBdPoint = [tracksAD(k).closestBdPoint(ii,1) tracksAD(k).closestBdPoint(ii,2)];
        curBdPointProjected = projPointOnLine(curBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.

        fromFirstBdPointToFirstAdh = [tracksAD(k).xCoord(sF)-firstBdPointProjected(1), tracksAD(k).yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
        fromFirstBdPointToLastAdh = [tracksAD(k).xCoord(ii)-firstBdPointProjected(1), tracksAD(k).yCoord(ii)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
        fromCurBdPointToFirstAdh = [tracksAD(k).xCoord(sF)-curBdPointProjected(1), tracksAD(k).yCoord(sF)-curBdPointProjected(2)]; % a vector from the last edge point to the first track point
        fromCurBdPointToLastAdh = [tracksAD(k).xCoord(ii)-curBdPointProjected(1), tracksAD(k).yCoord(ii)-curBdPointProjected(2)]; % a vector from the last edge point to the last track point
        firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
        curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
        if firstBDproduct>0 && firstBDproduct>curBDproduct% both adhesion points are in the same side
            tracksAD(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                            (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
            tracksAD(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                            (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
        else
            if curBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
                tracksAD(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                tracksAD(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
            else % Neither products are positive. This means the track crossed both the first and last boundaries. These would show shear movement. Relative comparison is performed.
                % Using actual BD points instead of projected ones because
                % somehow the track might be tilted...
                fromFirstBdPointToFirstAdh = [tracksAD(k).xCoord(sF)-tracksAD(k).closestBdPoint(sF,1), tracksAD(k).yCoord(sF)-tracksAD(k).closestBdPoint(sF,2)];
                fromFirstBdPointToLastAdh = [tracksAD(k).xCoord(ii)-tracksAD(k).closestBdPoint(sF,1), tracksAD(k).yCoord(ii)-tracksAD(k).closestBdPoint(sF,2)];
                fromCurBdPointToFirstAdh = [tracksAD(k).xCoord(sF)-tracksAD(k).closestBdPoint(ii,1), tracksAD(k).yCoord(sF)-tracksAD(k).closestBdPoint(ii,2)];
                fromCurBdPointToLastAdh = [tracksAD(k).xCoord(ii)-tracksAD(k).closestBdPoint(ii,1), tracksAD(k).yCoord(ii)-tracksAD(k).closestBdPoint(ii,2)];
                firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
                if firstBDproduct>curBDproduct % First BD point is in more distant position from the two adhesion points than the current BD point is.
                    tracksAD(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                    tracksAD(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                else
                    tracksAD(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                    tracksAD(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                end
            end
        end
    end
    % Record average advanceDist and edgeAdvanceDist for entire lifetime,
    % last 2 min, and every 2 minutes, and maximum of those
    % protrusion/retraction distance to figure out if the adhesion
    % experienced any previous protrusion ...
    for ii=sF:eF
        i2minBefore = max(sF,ii-frames2min);
        tracksAD(k).advanceDistChange2min(ii) = tracksAD(k).advanceDist(ii)-tracksAD(k).advanceDist(i2minBefore);
        tracksAD(k).edgeAdvanceDistChange2min(ii) = tracksAD(k).edgeAdvanceDist(ii)-tracksAD(k).edgeAdvanceDist(i2minBefore);
    end
    % Get the maximum of them. 
    tracksAD(k).maxAdvanceDistChange = max(tracksAD(k).advanceDistChange2min(sF:eF-1));
    tracksAD(k).maxEdgeAdvanceDistChange = max(tracksAD(k).edgeAdvanceDistChange2min(sF:eF-1));
%     lastBdPoint = [tracksNA(k).closestBdPoint(eF,1) tracksNA(k).closestBdPoint(eF,2)];
%     lastBdPointProjected = projPointOnLine(lastBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.
%     
%     fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
%     fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(eF)-firstBdPointProjected(1), tracksNA(k).yCoord(eF)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
%     fromLastBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-lastBdPointProjected(1), tracksNA(k).yCoord(sF)-lastBdPointProjected(2)]; % a vector from the last edge point to the first track point
%     fromLastBdPointToLastAdh = [tracksNA(k).xCoord(eF)-lastBdPointProjected(1), tracksNA(k).yCoord(eF)-lastBdPointProjected(2)]; % a vector from the last edge point to the last track point
%     firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
%     lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
%     if firstBDproduct>0 && firstBDproduct>lastBDproduct% both adhesion points are in the same side
%         tracksNA(k).advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
%         tracksNA(k).edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
%                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
%     else
%         if lastBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
%             tracksNA(k).advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                             (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
%             tracksNA(k).edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                             (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
%             % this code is nice to check:
% %             figure, imshow(paxImgStack(:,:,tracksNA(k).endingFrame),[]), hold on, plot(tracksNA(k).xCoord,tracksNA(k).yCoord,'w'), plot(tracksNA(k).closestBdPoint(:,1),tracksNA(k).closestBdPoint(:,2),'r')
% %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
% %             plot(tracksNA(k).xCoord(tracksNA(k).startingFrame),tracksNA(k).yCoord(tracksNA(k).startingFrame),'yo'),plot(tracksNA(k).xCoord(tracksNA(k).endingFrame),tracksNA(k).yCoord(tracksNA(k).endingFrame),'mo')
%         else
% %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
%             fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(sF,2)];
%             fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(sF,2)];
%             fromLastBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(eF,2)];
%             fromLastBdPointToLastAdh = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(eF,2)];
%             firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
%             lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
%             if firstBDproduct>lastBDproduct
%                 tracksNA(k).advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
%                 tracksNA(k).edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
%                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
%             else
%                 tracksNA(k).advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                                 (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
%                 tracksNA(k).edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
%                                                                 (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
%             end
%         end
%     end
    %Mean Squared Displacement
    meanX = mean(tracksAD(k).xCoord(logical(tracksAD(k).presence)));
    meanY = mean(tracksAD(k).yCoord(logical(tracksAD(k).presence)));
    tracksAD(k).MSD=sum((tracksAD(k).xCoord(logical(tracksAD(k).presence))'-meanX).^2+...
        (tracksAD(k).yCoord(logical(tracksAD(k).presence))'-meanY).^2);
    tracksAD(k).MSDrate = tracksAD(k).MSD/tracksAD(k).lifeTime;
    progressText(k/(numTracks-1),'Post-analysis:');
end
%% saving
save([dataPath filesep 'tracksAD.mat'], 'tracksAD','-v7.3')
%% saving
if saveAnalysis
    % saving
    %% NA FA Density analysis
    numNAs = zeros(nFrames,1);
    for ii=1:nFrames
        numNAs(ii) = sum(arrayfun(@(x) (x.presence(ii)==true), tracksAD));
        NADensity(ii) = numNAs(ii)/(bandArea(ii)*MD.pixelSize_^2*1e6);  % unit: number/um2 
        FADensity(ii) = numFAs(ii)/(bandArea(ii)*MD.pixelSize_^2*1e6);  % unit: number/um2 
    end
    save([dataPath filesep 'NAFADensity.mat'], 'NADensity','FADensity')

    %% Lifetime analysis
    p=0;
    idx = false(numel(tracksAD),1);
    for k=1:numel(tracksAD)
        % look for tracks that had a state of 'BA' and become 'NA'
        firstNAidx = find(strcmp(tracksAD(k).state,'NA'),1,'first');
        % see if the state is 'BA' before 'NA' state
        if (~isempty(firstNAidx) && firstNAidx>1 && strcmp(tracksAD(k).state(firstNAidx-1),'BA')) || (~isempty(firstNAidx) &&firstNAidx==1)
            p=p+1;
            idx(k) = true;
            tracksAD(k).emerging = true;
            tracksAD(k).emergingFrame = firstNAidx;
        else
            tracksAD(k).emerging = false;
        end        
    end
    %% Analysis of those whose force was under noise level: how long does it take
    % Analysis shows that force is already developed somewhat compared to
    % background. 
    % Filter out any tracks that has 'Out_of_ROI' in their status (especially
    % after NA ...)
    trNAonly = tracksAD(idx);
    indNAMature = false(numel(trNAonly));
    indNAFail = false(numel(trNAonly));
    p=0; q=0;

    for k=1:numel(trNAonly)
        if trNAonly(k).emerging 
            % maturing NAs
            if (any(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'FC')) || ...
                    any(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'FA'))) && ...
                    sum(trNAonly(k).presence)>8

                trNAonly(k).maturing = true;
                indNAMature(k) = true;
                p=p+1;
                % lifetime until FC
                lifeTimeNAmaturing(p) = sum(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'NA'));
                % it might be beneficial to store amplitude time series. But
                % this can be done later from trackNAmature
            elseif sum(tracksAD(k).presence)<61 && sum(tracksAD(k).presence)>6
            % failing NAs
                tracksAD(k).maturing = false;
                indNAFail(k) = true;
                q=q+1;
                % lifetime until FC
                lifeTimeNAfailing(q) = sum(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'NA'));
            end
        end
    end
    maturingRatio = p/(p+q);
    tracksNAmaturing = trNAonly(indNAMature);
    tracksNAfailing = trNAonly(indNAFail);
    save([dataPath filesep 'allData.mat'], 'trNAonly', 'tracksNAfailing','tracksNAmaturing','maturingRatio','lifeTimeNAfailing','lifeTimeNAmaturing')
else
    trNAonly = tracksAD;
end
disp('Done!')
end
function newTracks = formatTracks(tracks,detectedNAs,nFrames)
% Format tracks structure into tracks with every frame

newTracks(numel(tracks),1) = struct('xCoord', [], 'yCoord', [],'state',[],'iFrame',[],'presence',[],'amp',[],'bkgAmp',[]);
% BA: before adhesion, NA: nascent adh, FC: focal complex, FA: focal adh,
% ANA: after NA (failed to be matured.
for i = 1:numel(tracks)
    % Get the x and y coordinate of all compound tracks
    startNA = true;
    endNA = true;
    for  jj = 1 : nFrames
        newTracks(i).iFrame(jj) = jj;
        if jj<tracks(i).seqOfEvents(1,1)
            newTracks(i).state{jj} = 'BA';
            newTracks(i).xCoord(jj) = NaN;
            newTracks(i).yCoord(jj) = NaN;
            newTracks(i).presence(jj) = false;
            newTracks(i).amp(jj) = NaN;
            newTracks(i).bkgAmp(jj) = NaN;
        elseif jj>tracks(i).seqOfEvents(end,1)
            newTracks(i).state{jj} = 'ANA';
            newTracks(i).xCoord(jj) = NaN;
            newTracks(i).yCoord(jj) = NaN;
            newTracks(i).amp(jj) = NaN;
            newTracks(i).bkgAmp(jj) = NaN;
            newTracks(i).presence(jj) = false;
            if endNA
                newTracks(i).endingFrame = jj-1;
                endNA = false;
            end
        elseif jj==tracks(i).seqOfEvents(2,1)
            newTracks(i).state{jj} = 'NA';
            newTracks(i).xCoord(jj) = tracks(i).tracksCoordAmpCG(1,1+8*(jj-tracks(i).seqOfEvents(1,1)));
            newTracks(i).yCoord(jj) = tracks(i).tracksCoordAmpCG(1,2+8*(jj-tracks(i).seqOfEvents(1,1)));
            newTracks(i).amp(jj) = tracks(i).tracksCoordAmpCG(1,4+8*(jj-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(jj) = NaN;
            else
                jFromBirth=jj-tracks(i).seqOfEvents(1,1)+1;
                newTracks(i).bkgAmp(jj) = detectedNAs(jj).bkg(tracks(i).tracksFeatIndxCG(jFromBirth));
                newTracks(i).sigma(jj) = detectedNAs(jj).sigmaX(tracks(i).tracksFeatIndxCG(jFromBirth));
            end
            newTracks(i).presence(jj) = true;
            if startNA
                newTracks(i).startingFrame = jj;
                startNA = false;
            end
            if endNA
                newTracks(i).endingFrame = jj;
                endNA = false;
            end
        else
            newTracks(i).state{jj} = 'NA';
            newTracks(i).xCoord(jj) = tracks(i).tracksCoordAmpCG(1,1+8*(jj-tracks(i).seqOfEvents(1,1)));
            newTracks(i).yCoord(jj) = tracks(i).tracksCoordAmpCG(1,2+8*(jj-tracks(i).seqOfEvents(1,1)));
            newTracks(i).amp(jj) = tracks(i).tracksCoordAmpCG(1,4+8*(jj-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(jj) = NaN;
            else % tracksFeatIndxCG: [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%           .tracksFeatIndxCG: Connectivity matrix of features between
%                              frames, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = number of frames
%                              the compound track spans. Zeros indicate
%                              frames where track segments do not exist
%                              (either because those frames are before the
%                              segment starts or after it ends, or because
%                              of losing parts of a segment.
%                 newTracks(i).bkgAmp(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).bkg(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
%                 newTracks(i).sigma(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).sigmaX(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).bkgAmp(jj) = detectedNAs(jj).bkg(tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(jj) = detectedNAs(jj).sigmaX(tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(jj) = true;
            if startNA
                newTracks(i).startingFrame = jj;
                startNA = false;
            end
        end
            
        if isfield(tracks, 'label'),
            newTracks(iTrack).label = tracks(i).label;
        end
    end
    % go through frames again and fill NaNs with numbers at the gap
    % position
    jj=newTracks(i).startingFrame;
    gap=0;
    while jj<newTracks(i).endingFrame
        if newTracks(i).presence(jj) && ~isnan(newTracks(i).xCoord(jj))
            % jump to the next broken block
            nNextConsecBlock = find(isnan(newTracks(i).xCoord) & newTracks(i).iFrame>jj,1,'first');
            if isempty(nNextConsecBlock)
                break % there is no gap afterward
            else
                % see if abscence (or gap) is until the end of frame
                % find the next presence after this gap
                nNextNextConsecBlock = find(~isnan(newTracks(i).xCoord) & newTracks(i).iFrame>nNextConsecBlock,1,'first');
                gap = nNextNextConsecBlock-nNextConsecBlock;
                jj=nNextConsecBlock;
            end
        else
            for kk=1:gap
                newTracks(i).xCoord(jj+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(jj-1)+kk*newTracks(i).xCoord(jj+gap))/(gap+1);
                newTracks(i).yCoord(jj+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(jj-1)+kk*newTracks(i).yCoord(jj+gap))/(gap+1);
                newTracks(i).amp(jj+kk-1) = ((gap+1-kk)*newTracks(i).amp(jj-1)+kk*newTracks(i).amp(jj+gap))/(gap+1);
                newTracks(i).bkgAmp(jj+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(jj-1)+kk*newTracks(i).bkgAmp(jj+gap))/(gap+1);
            end
            jj=jj+gap;
        end
    end
    if isempty(newTracks(i).startingFrame)
        warning(['startingFrame is empty for track #' num2str(i)])
    end
end
end

