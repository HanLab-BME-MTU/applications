function [] = analyzeAdhesionMaturationProcess(movieData,varargin)
% [tracksNA,lifeTimeNA] = analyzeAdhesionMaturation(pathForTheMovieDataFile,outputPath,showAllTracks,plotEachTrack)
% filter out NA tracks, obtain life time of each NA tracks

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

% Sangyoon Han April 2014

%% Inputs
% ip =inputParser;
% ip.addRequired('pathForTheMovieDataFile',@(x)ischar(x)||isa(x,MovieData))
% ip.addParamValue('showAllTracks',false,@(x)islogical(x)||isempty(x))
% ip.addParamValue('plotEachTrack',false,@(x)islogical(x)||isempty(x))
% ip.addParamValue('onlyEdge',false,@islogical); % collect NA tracks that ever close to cell edge
% ip.addParamValue('outputPath','AdhesionTracking',@ischar)
% ip.addParamValue('saveAnalysis',true,@islogical)
% ip.addParamValue('matchWithFA',true,@islogical) %For cells with only NAs, we turn this off.
% ip.addParamValue('minLifetime',5,@isscalar) %For cells with only NAs, we turn this off.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;
%Parse input, store in parameter structure
%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('AdhesionAnalysisProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(AdhesionAnalysisProcess(movieData));                                                                                                 
end
adhAnalProc = movieData.processes_{iProc};
paramsIn = parseProcessParams(adhAnalProc,paramsIn);

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(movieData,showAllTracks,plotEachTrack,varargin{:});
outputPath=paramsIn.outputPath;
saveAnalysis=paramsIn.saveAnalysis;
matchWithFA=paramsIn.matchWithFA;
minLifetime=paramsIn.minLifetime;
showAllTracks=paramsIn.showAllTracks;
plotEachTrack=paramsIn.plotEachTrack;
onlyEdge=paramsIn.onlyEdge;
iChans=paramsIn.iChan;
% Load the Paxillin channel

%% Data Set up
% Load the MovieData
% movieDataPath = [movieData '/movieData.mat'];
MD = movieData; %MovieData.load(movieDataPath);
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

iiformat = ['%.' '3' 'd'];
%     paxLevel = zeros(nFrames,1);
% SegmentationPackage = MD.getPackage(MD.getPackageIndex('SegmentationPackage'));
% minSize = round((500/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=.5um*.5um
minLifetime = min(nFrames,minLifetime);
markerSize = 2;
% tracks
% iPaxChannel = 1; % this should be intentionally done in the analysis level
% if ~trackNAProc.checkChannelOutput(iPaxChannel)
%     iPaxChannel = 2;
% end
% From here I have to go by channel
for iPaxChan=iChans
    disp(['-----------------------Working on channel: ' num2str(iPaxChan) ' out of channels ' num2str(iChans) '....-------------------------------------------'])
    % Set up the output file path per channel
    [~,iChanName]=fileparts(movieData.channels_(iPaxChan).channelPath_);
    outputFilePath = [paramsIn.OutputDirectory filesep outputPath filesep iChanName];
    dataPath = [outputFilePath filesep 'data'];
    paxPath = [outputFilePath filesep 'pax'];
    paxtifPath = [outputFilePath filesep 'paxtifs'];
    epsPath = [outputFilePath filesep 'eps'];
    figPath = [outputFilePath filesep 'figs'];
    if ~exist(paxtifPath,'dir') || ~exist(paxPath,'dir') || ~exist(figPath,'dir') || ~exist(epsPath,'dir') || ~exist(dataPath,'dir') 
        try
            mkdir(paxPath);
            mkdir(paxtifPath);
            mkdir(figPath);
            mkdir(epsPath);
            mkdir(dataPath);
        catch
            system(['mkdir -p ' paxPath]);
            system(['mkdir -p ' paxtifPath]);
            system(['mkdir -p ' figPath]);
            system(['mkdir -p ' epsPath]);
            system(['mkdir -p ' dataPath]);
        end
    end
    % filter out tracks that have lifetime less than 2 frames
    disp('loading NA tracks...')
    tic
    tracksNAorg = trackNAProc.loadChannelOutput(iPaxChan);
    toc
    % if ii>minLifetime
    SEL = getTrackSEL(tracksNAorg); %SEL: StartEndLifetime
    % Remove any less than 3-frame long track.
    isValid = SEL(:,3) >= minLifetime;
    tracksNAorg = tracksNAorg(isValid);
    % end
    detectedNAs = detectedNAProc.loadChannelOutput(iPaxChan);

    % re-express tracksNA so that each track has information for every frame
    disp('reformating NA tracks...')
    tic
    tracksNA = formatTracks(tracksNAorg,detectedNAs,nFrames); 
    toc

    % disp('loading segmented FAs...')
    iFASeg = 6;
    FASegProc = FASegPackage.processes_{iFASeg};
    bandArea = zeros(nFrames,1);
    NADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
    FADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
    numFAs = zeros(nFrames,1);
%     nChannels = numel(MD.channels_);
%     iChan = 0;
    % Finding which channel has a cell mask information
    maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
%     for k=1:nChannels
%         if maskProc.checkChannelOutput(k)
%             iChan = k;
%         end
%     end
    % Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered
    if onlyEdge
        disp(['Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered'])
        trackIdx = false(numel(tracksNA),1);
        bandwidthNA = 7; %um 
        bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
        for ii=1:nFrames
            % Cell Boundary Mask 
            mask = maskProc.loadChannelOutput(iPaxChan,ii);
            % mask for band from edge
            iMask = imcomplement(mask);
            distFromEdge = bwdist(iMask);
            bandMask = distFromEdge <= bandwidthNA_pix;

            maskOnlyBand = bandMask & mask;
            bandArea(ii) = sum(maskOnlyBand(:)); % in pixel

            % collect index of tracks which first appear at frame ii
            idxFirstAppear = arrayfun(@(x) x.startingFrame==ii,tracksNA);
            % now see if these tracks ever in the maskOnlyBand
            for k=find(idxFirstAppear)'
                if maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii)))
                    trackIdx(k) = true;
                end
            end
        end    
    else
        disp('Entire adhesion tracks are considered.')
        trackIdx = true(numel(tracksNA),1);
        mask = maskProc.loadChannelOutput(iPaxChan,1);
    %     bandwidthNA = 5; %um 
    %     bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    %     for ii=1:nFrames
    %         % Cell Boundary Mask 
    %         mask = maskProc.loadChannelOutput(iPaxChan,ii);
    %         % mask for band from edge
    %         iMask = imcomplement(mask);
    %         distFromEdge = bwdist(iMask);
    %         bandMask = distFromEdge <= bandwidthNA_pix;
    % 
    %         maskOnlyBand = bandMask & mask;
    %         bandArea(ii) = sum(maskOnlyBand(:)); % in pixel
    %         % filter tracks with naMasks
    %         % only deal with presence and status
    %         % Tracks in its emerging state ever overlap with bandMask are
    %         % considered.
    %         for k=1:numel(tracksNA)
    %             if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
    %                     ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
    %                     round(tracksNA(k).xCoord(ii)) < 1 || ...
    %                     round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
    %                     round(tracksNA(k).yCoord(ii)) < 1) || ...
    %                     ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
    %                 tracksNA(k).state{ii} = 'Out_of_Band';
    %                 tracksNA(k).presence(ii) = false;
    %                 if trackIdx(k)
    %                     trackIdx(k) = false;
    %                 end
    %             end
    %         end
    %     end
    end
    % get rid of tracks that have out of bands...
    tracksNA = tracksNA(trackIdx);
    %% add intensity of tracks including before and after NA status
    % get the movie stack
    [h,w]=size(mask);
    imgStack = zeros(h,w,nFrames);
    for ii=1:nFrames
        imgStack(:,:,ii)=MD.channels_(iPaxChan).loadImage(ii); 
    end
    % get the intensity
    disp('Reading intensities with additional tracking...')
    tic
    tracksNA = readIntensityFromTracks(tracksNA,imgStack,1,'extraLength',20); % 1 means intensity collection from pax image
    toc
    firstMask=maskProc.loadChannelOutput(iPaxChan,1);
    cropMaskStack = false(size(firstMask,1),size(firstMask,2),nFrames);
    numTracks=numel(tracksNA);
    progressText(0,'Matching with segmented adhesions:');
    %% Matching with segmented adhesions
    for ii=1:nFrames
        % Cell Boundary Mask 
        if length(iChans)>1
            %Combine the the multiple masks to one
            maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,ii),iChans,'UniformOutput',false);
            maskAll=reshape(cell2mat(maskEach),size(firstMask,1),size(firstMask,2),[]);
            mask = any(maskAll,3);
        elseif length(iChans)==1
            mask = maskProc.loadChannelOutput(iPaxChan,ii);
        end
        % Cell Boundary
        [B,~,nBD]  = bwboundaries(mask,'noholes');
        cropMaskStack(:,:,ii) = mask; %maskProc.loadChannelOutput(iPaxChan,ii);
        if matchWithFA
            % filter tracks with naMasks
    %         disp(['Processing ' num2str(ii) 'th frame out of ' num2str(nFrames) ' total frames...'])
            % focal contact (FC) analysis
            % Get the mask for FAs
            maskFAs = FASegProc.loadChannelOutput(iPaxChan,ii);
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
            for k=1:numel(tracksNA)
                if tracksNA(k).presence(ii)
                    if ~strcmp(tracksNA(k).state{ii} , 'NA') && ii>1
                        tracksNA(k).state{ii} = tracksNA(k).state{ii-1};
                    end
                    % decide if each track is associated with FC or FA
                    p = 0;
                    for jj=FCIdx'
                        p=p+1;
                        if any(round(tracksNA(k).xCoord(ii))==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))==Adhs(jj).PixelList(:,2))
                            tracksNA(k).state{ii} = 'FC';
                            tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
                            tracksNA(k).FApixelList{ii} = Adhs(jj).PixelList;
                            tracksNA(k).adhBoundary{ii} = adhBound{jj};
                            tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii)));
                        end
                    end
                    p = 0;
                    for jj=FAIdx'
                        p=p+1;
                        if any(round(tracksNA(k).xCoord(ii))==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))==Adhs(jj).PixelList(:,2))
                            tracksNA(k).state{ii} = 'FA';
                            tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
                            tracksNA(k).FApixelList{ii} = Adhs(jj).PixelList;
                            tracksNA(k).adhBoundary{ii} = adhBound{jj};
                            tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii)));
                        end
                    end
                elseif ii>tracksNA(k).endingFrame && (strcmp(tracksNA(k).state{tracksNA(k).endingFrame},'FA')...
                        || strcmp(tracksNA(k).state{tracksNA(k).endingFrame},'FC'))
                    % starting from indexed maskFAs, find out segmentation that is
                    % closest to the last track point.
                    subMaskFAs = maskFAs==tracksNA(k).faID(tracksNA(k).endingFrame);
                    if max(subMaskFAs(:))==0
                        tracksNA(k).state{ii} = 'ANA';
                        tracksNA(k).FApixelList{ii} = NaN;
                        tracksNA(k).adhBoundary{ii} = NaN;
                        continue
                    else
                        propSubMaskFAs = regionprops(subMaskFAs,'PixelList');
                        minDist = zeros(length(propSubMaskFAs),1);
                        for q=1:length(propSubMaskFAs)
                            minDist(q) = min(sqrt((propSubMaskFAs(q).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame))).^2 +...
                                (propSubMaskFAs(q).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame))).^2));
                        end
                        % find the closest segment
                        [~,subMaskFAsIdx] = min(minDist);
                        subAdhBound = bwboundaries(subMaskFAs,'noholes');    
                        [~,closestPixelID] = min(sqrt((propSubMaskFAs(subMaskFAsIdx).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame))).^2 +...
                            (propSubMaskFAs(subMaskFAsIdx).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame))).^2));

                        tracksNA(k).state{ii} = 'FC';
    %                     tracksNA(k).xCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,1);
    %                     tracksNA(k).yCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,2);
                        tracksNA(k).FApixelList{ii} = propSubMaskFAs(subMaskFAsIdx).PixelList;
                        tracksNA(k).adhBoundary{ii} = subAdhBound{subMaskFAsIdx};
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
        for k=1:numel(tracksNA)
            % distance to the cell edge
    %         if tracksNA(k).presence(ii)
            if ii>=tracksNA(k).startingFrameExtraExtra && ii<=tracksNA(k).endingFrameExtraExtra
                xCropped = tracksNA(k).xCoord(ii);
                yCropped = tracksNA(k).yCoord(ii);
                distToAdh = sqrt(sum((allBdPoints- ...
                    ones(size(allBdPoints,1),1)*[xCropped, yCropped]).^2,2));
                [minDistToBd,indMinBdPoint] = min(distToAdh);
                tracksNA(k).distToEdge(ii) = minDistToBd;
                tracksNA(k).closestBdPoint(ii,:) = allBdPoints(indMinBdPoint,:); % this is lab frame of reference. (not relative to adhesion position)
            end
        end

        paxImageCropped=MD.channels_(iPaxChan).loadImage(ii); 
        if showAllTracks
            h2=figure;
            %Scale bar 2 um
        %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
            % size of the region of interest
%             if ii==1
%                 [imSizeY,imSizeX] = size(paxImageCropped);
%             end

            paxImageCroppedInverted = imcomplement(paxImageCropped);
            minPax = min(paxImageCroppedInverted(:));
            maxPax = max(paxImageCroppedInverted(:));

    %         if ii==1
    %             minPax1 = 1*minPax;
    %             minPax2 = uint16(double(minPax)+double(0.25*(maxPax-minPax)));
    %             hPaxTemp = figure;
    %             subplot(1,2,1),imshow(paxImageCroppedInverted,[minPax1 maxPax]),title(['minPax1 = ' num2str(minPax1) ]);
    %             subplot(1,2,2),imshow(paxImageCroppedInverted,[minPax2 maxPax]),title(['minPax2 = ' num2str(minPax2) ]);
    %             minPax = input('type desired minPax for maximum of the image: ');
    %             close(hPaxTemp);
    %         end        
    %         imshow(paxImageCroppedInverted,[minPax maxPax]), hold on
            imshow(paxImageCroppedInverted,[minPax+0.7*(maxPax-minPax) maxPax]), hold on
            line([10 10+round(2000/MD.pixelSize_)],[15 15],'LineWidth',2,'Color',[0,0,0])

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
                plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
            end
            for k=1:numel(tracksNA)
                if tracksNA(k).presence(ii)
                    if strcmp(tracksNA(k).state{ii} , 'NA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii),tracksNA(k).yCoord(1:ii),'r', 'LineWidth', 0.5)
                        plot(tracksNA(k).xCoord(ii),tracksNA(k).yCoord(ii),'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{ii} , 'FC')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii),tracksNA(k).yCoord(1:ii),'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                        plot(tracksNA(k).xCoord(ii),tracksNA(k).yCoord(ii),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{ii} , 'FA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii),tracksNA(k).yCoord(1:ii),'k', 'LineWidth', 0.5)
                        plot(tracksNA(k).xCoord(ii),tracksNA(k).yCoord(ii),'ko','MarkerSize',markerSize, 'LineWidth', 0.5)
                    end
                end
            end

            print(h2, '-depsc2', strcat(epsPath,'/pax',num2str(ii,iiformat),'.eps'));
            print(h2, '-dtiff', strcat(paxtifPath,'/pax',num2str(ii,iiformat),'.tif'));
        %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
            hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
            close(h2)
            clear h2
        end
        if saveAnalysis
            imwrite(paxImageCropped,strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
        end
        progressText(ii/(nFrames-1),'Matching with segmented adhesions:');
    end
    %% disp('Intermediate saving before post analysis...')
    disp('Intermediate saving before post analysis...')
    save([dataPath filesep 'tracksNA.mat'], 'tracksNA','-v7.3')
%     save([dataPath filesep 'intermediateWorkspace.mat'], '-v7.3')
    %% protrusion/retraction information
    % time after protrusion onset (negative value if retraction, based
    % on the next protrusion onset) in frame, based on tracksNA.distToEdge
    % First I have to quantify when the protrusion and retraction onset take
    % place.
    for ii=1:numTracks
        idxZeros = tracksNA(ii).closestBdPoint(:)==0;
        tracksNA(ii).closestBdPoint(idxZeros)=NaN;
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
        presIdx = logical(tracksNA(k).presence);
        % get the instantaneous velocity
        % get the distance first
        curTrackVelLength=sum(presIdx)-1;
        distTrajec=zeros(curTrackVelLength,1);
        presIdxSeq = find(presIdx);
        for kk=1:curTrackVelLength
            real_kk = presIdxSeq(kk);
            distTrajec(kk) = sqrt(sum((tracksNA(k).closestBdPoint(real_kk+1,:)- ...
                tracksNA(k).closestBdPoint(real_kk,:)).^2,2));
            lastPointIntX = round(tracksNA(k).closestBdPoint(real_kk+1,1));
            lastPointIntY = round(tracksNA(k).closestBdPoint(real_kk+1,2));
            if cropMaskStack(lastPointIntY,lastPointIntX,real_kk) %if the last point is in the first mask, it is inward
                distTrajec(kk) = -distTrajec(kk);
            end
        end
        if any(distTrajec~=0)
            [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
            if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
                tracksNA(k).isProtrusion = true;
            else
                tracksNA(k).isProtrusion = false;
            end
            % average velocity (positive for protrusion)
            curProtVel = (Protrusion.Veloc); curProtVel(isnan(curProtVel))=0;
            curProtPersTime = (Protrusion.persTime); curProtPersTime(isnan(curProtPersTime))=0;
            curRetVel = (Retraction.Veloc); curRetVel(isnan(curRetVel))=0;
            curRetPersTime = (Retraction.persTime); curRetPersTime(isnan(curRetPersTime))=0;

            tracksNA(k).edgeVel = (mean(curProtVel.*curProtPersTime)-mean(curRetVel.*curRetPersTime))/mean([curProtPersTime;curRetPersTime]);
        else
            tracksNA(k).edgeVel = 0;
        end
        % lifetime information
        try
            sFextend=tracksNA(k).startingFrameExtraExtra;
            eFextend=tracksNA(k).endingFrameExtraExtra;
            sF=tracksNA(k).startingFrameExtra;
            eF=tracksNA(k).endingFrameExtra;
            if isempty(sF)
                sF=tracksNA(k).startingFrame;
            end
            if isempty(eF)
                eF=tracksNA(k).endingFrame;
            end
        catch
            sF=tracksNA(k).startingFrame;
            eF=tracksNA(k).endingFrame;
            tracksNA(k).lifeTime = eF-sF;    
        end
        tracksNA(k).lifeTime = eF-sF;    
        % Inital intensity slope for one min
        timeInterval = deltaT/60; % in min
        earlyPeriod = floor(1/timeInterval); % frames per minute
        lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
        lastFrameFromOne = lastFrame - sF+1;
    %     lastFrameFromOne = sF;
    %     lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
        [curR,curM] = regression(timeInterval*(1:lastFrameFromOne),tracksNA(k).amp(sF:lastFrame));
        tracksNA(k).ampSlope = curM; % in a.u./min
        tracksNA(k).ampSlopeR = curR; % Pearson's correlation coefficient

        curEndFrame = min(tracksNA(k).startingFrameExtra+periodFrames-1,tracksNA(k).endingFrame);
        curEarlyPeriod = curEndFrame - tracksNA(k).startingFrameExtra+1;
        [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:curEndFrame));
        tracksNA(k).earlyAmpSlope = curM; % in a.u./min

        curStartFrame = max(tracksNA(k).startingFrame,tracksNA(k).endingFrameExtra-periodFrames+1);
        curLatePeriod = tracksNA(k).endingFrameExtra - curStartFrame+1;
        [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksNA(k).ampTotal(curStartFrame:tracksNA(k).endingFrameExtra));
        tracksNA(k).lateAmpSlope = curMlate; % in a.u./min

        curEndFrame = min(sF+periodFrames-1,eF);
        curEarlyPeriod = curEndFrame - sF+1;
        [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).distToEdge(sF:curEndFrame));
        tracksNA(k).distToEdgeSlope = curMdist; % in a.u./min
        tracksNA(k).distToEdgeChange = (tracksNA(k).distToEdge(end)-tracksNA(k).distToEdge(tracksNA(k).startingFrame)); % in pixel
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
        fitobj = fit(tracksNA(k).xCoord(sF:eF)',tracksNA(k).yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
        x0=nanmedian(tracksNA(k).xCoord);
        y0=fitobj(x0);
        dx = 1;
        dy = fitobj.p1;
        trackLine = createLineGeom2d(x0,y0,dx,dy); % this is a geometric average linear line of the adhesion track

    %     trackLine = edgeToLine(edge);
        firstBdPoint = [tracksNA(k).closestBdPoint(sF,1) tracksNA(k).closestBdPoint(sF,2)];
        firstBdPointProjected = projPointOnLine(firstBdPoint, trackLine); % this is an edge boundary point at the first time point projected on the average line of track.
        % try to record advanceDist and edgeAdvanceDist for every single time
        % point ...
        for ii=sFextend:eFextend
            curBdPoint = [tracksNA(k).closestBdPoint(ii,1) tracksNA(k).closestBdPoint(ii,2)];
            curBdPointProjected = projPointOnLine(curBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.

            fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
            fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(ii)-firstBdPointProjected(1), tracksNA(k).yCoord(ii)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
            fromCurBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-curBdPointProjected(1), tracksNA(k).yCoord(sF)-curBdPointProjected(2)]; % a vector from the last edge point to the first track point
            fromCurBdPointToLastAdh = [tracksNA(k).xCoord(ii)-curBdPointProjected(1), tracksNA(k).yCoord(ii)-curBdPointProjected(2)]; % a vector from the last edge point to the last track point
            firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
            curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
            if firstBDproduct>0 && firstBDproduct>curBDproduct% both adhesion points are in the same side
                tracksNA(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
            else
                if curBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
                    tracksNA(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                    tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                    % this code is nice to check:
        %             figure, imshow(paxImgStack(:,:,tracksNA(k).endingFrame),[]), hold on, plot(tracksNA(k).xCoord,tracksNA(k).yCoord,'w'), plot(tracksNA(k).closestBdPoint(:,1),tracksNA(k).closestBdPoint(:,2),'r')
        %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
        %             plot(tracksNA(k).xCoord(tracksNA(k).startingFrame),tracksNA(k).yCoord(tracksNA(k).startingFrame),'yo'),plot(tracksNA(k).xCoord(tracksNA(k).endingFrame),tracksNA(k).yCoord(tracksNA(k).endingFrame),'mo')
                else % Neither products are positive. This means the track crossed both the first and last boundaries. These would show shear movement. Relative comparison is performed.
        %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
                    % Using actual BD points instead of projected ones because
                    % somehow the track might be tilted...
                    fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(sF,2)];
                    fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(ii)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(ii)-tracksNA(k).closestBdPoint(sF,2)];
                    fromCurBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(ii,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(ii,2)];
                    fromCurBdPointToLastAdh = [tracksNA(k).xCoord(ii)-tracksNA(k).closestBdPoint(ii,1), tracksNA(k).yCoord(ii)-tracksNA(k).closestBdPoint(ii,2)];
                    firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                    curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
                    if firstBDproduct>curBDproduct % First BD point is in more distant position from the two adhesion points than the current BD point is.
                        tracksNA(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                        tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                    else
                        tracksNA(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                        tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
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
            tracksNA(k).advanceDistChange2min(ii) = tracksNA(k).advanceDist(ii)-tracksNA(k).advanceDist(i2minBefore);
            tracksNA(k).edgeAdvanceDistChange2min(ii) = tracksNA(k).edgeAdvanceDist(ii)-tracksNA(k).edgeAdvanceDist(i2minBefore);
        end
        % Get the maximum of them. 
        tracksNA(k).maxAdvanceDistChange = max(tracksNA(k).advanceDistChange2min(sF:eF-1));
        tracksNA(k).maxEdgeAdvanceDistChange = max(tracksNA(k).edgeAdvanceDistChange2min(sF:eF-1));
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
        meanX = mean(tracksNA(k).xCoord(logical(tracksNA(k).presence)));
        meanY = mean(tracksNA(k).yCoord(logical(tracksNA(k).presence)));
        tracksNA(k).MSD=sum((tracksNA(k).xCoord(logical(tracksNA(k).presence))'-meanX).^2+...
            (tracksNA(k).yCoord(logical(tracksNA(k).presence))'-meanY).^2);
        tracksNA(k).MSDrate = tracksNA(k).MSD/tracksNA(k).lifeTime;
        progressText(k/(numTracks-1),'Post-analysis:');
    end
    %% saving
    save([dataPath filesep 'tracksNA.mat'], 'tracksNA','-v7.3')
    %% saving
    if saveAnalysis
        % saving
        %% NA FA Density analysis
        numNAs = zeros(nFrames,1);
        for ii=1:nFrames
            numNAs(ii) = sum(arrayfun(@(x) (x.presence(ii)==true), tracksNA));
            NADensity(ii) = numNAs(ii)/(bandArea(ii)*MD.pixelSize_^2*1e6);  % unit: number/um2 
            FADensity(ii) = numFAs(ii)/(bandArea(ii)*MD.pixelSize_^2*1e6);  % unit: number/um2 
        end
        save([dataPath filesep 'NAFADensity.mat'], 'NADensity','FADensity')

        %% Lifetime analysis
        p=0;
        idx = false(numel(tracksNA),1);
        for k=1:numel(tracksNA)
            % look for tracks that had a state of 'BA' and become 'NA'
            firstNAidx = find(strcmp(tracksNA(k).state,'NA'),1,'first');
            % see if the state is 'BA' before 'NA' state
            if (~isempty(firstNAidx) && firstNAidx>1 && strcmp(tracksNA(k).state(firstNAidx-1),'BA')) || (~isempty(firstNAidx) &&firstNAidx==1)
                p=p+1;
                idx(k) = true;
                tracksNA(k).emerging = true;
                tracksNA(k).emergingFrame = firstNAidx;
            else
                tracksNA(k).emerging = false;
            end        
        end
        %% Analysis of those whose force was under noise level: how long does it take
        % Analysis shows that force is already developed somewhat compared to
        % background. 
        % Filter out any tracks that has 'Out_of_ROI' in their status (especially
        % after NA ...)
        trNAonly = tracksNA(idx);
        indMature = false(numel(trNAonly));
        indFail = false(numel(trNAonly));
        p=0; q=0;

        for k=1:numel(trNAonly)
            if trNAonly(k).emerging 
                % maturing NAs
                if (any(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'FC')) || ...
                        any(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'FA'))) && ...
                        sum(trNAonly(k).presence)>8

                    trNAonly(k).maturing = true;
                    indMature(k) = true;
                    p=p+1;
                    % lifetime until FC
                    lifeTimeNAmaturing(p) = sum(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'NA'));
                    % it might be beneficial to store amplitude time series. But
                    % this can be done later from trackNAmature
                elseif sum(tracksNA(k).presence)<61 && sum(tracksNA(k).presence)>6
                % failing NAs
                    tracksNA(k).maturing = false;
                    indFail(k) = true;
                    q=q+1;
                    % lifetime until FC
                    lifeTimeNAfailing(q) = sum(strcmp(trNAonly(k).state(trNAonly(k).emergingFrame:end),'NA'));
                end
            end
        end
        maturingRatio = p/(p+q);
        tracksNAmaturing = trNAonly(indMature);
        tracksNAfailing = trNAonly(indFail);
        save([dataPath filesep 'allData.mat'], 'trNAonly', 'tracksNAfailing','tracksNAmaturing','maturingRatio','lifeTimeNAfailing','lifeTimeNAmaturing')
    else
        trNAonly = tracksNA;
    end

    %% Run this separately with loading allData.mat if something failed 
    if plotEachTrack
        paxImageCropped=MD.channels_(iPaxChan).loadImage(1); 
        % size of the region of interest
        [imSizeY,imSizeX] = size(paxImageCropped);
        r1 = 50;
        h2=figure;
        for k=1:numel(trNAonly)%[7 15 32 127 129]%afterFAK[1 2 18 30 36 11 12 14 39 41 44]% 
            % try to crop window around the track
            if isempty(trNAonly(k).startingFrame) 
                continue
            end
            if strcmp(trNAonly(k).state{trNAonly(k).startingFrame},'FC') || strcmp(trNAonly(k).state{trNAonly(k).startingFrame},'FA')
                continue
            end
            fstart = max(trNAonly(k).startingFrame-20,1);
            fend = min(trNAonly(k).endingFrame,nFrames);
            iSF = trNAonly(k).startingFrame;
            wRoi = min(trNAonly(k).xCoord(iSF),r1)...
                +min(imSizeX+1-(trNAonly(k).xCoord(iSF)),r1);
            hRoi = min(trNAonly(k).yCoord(iSF),r1)...
                +min(imSizeY+1-(trNAonly(k).yCoord(iSF)),r1);
            set(h2,'Units','inches')
            set(h2,'PaperPositionMode','auto')
            set(h2, 'Position', [1,1,wRoi/(hRoi), 1])

            eachPaxPath = [paxPath filesep '/track' num2str(k,iiformat)];
            eachEpsPath = [epsPath filesep '/track' num2str(k,iiformat)];
            if ~exist(eachPaxPath,'dir') || ~exist(eachEpsPath,'dir')
                mkdir(eachPaxPath);
                mkdir(eachEpsPath);
            end

            for j=fstart:fend
                ha1 = get(h2,'CurrentAxes');%subplot('position',[0  0.5  1  0.5]);
                if j==fstart
                    paxImageCropped = imread(strcat(paxPath,'/pax',num2str(j,iiformat),'.tif'));
                    xminROI = round(max(1,trNAonly(k).xCoord(iSF)-(r1-1))); 
                    xmaxROI = round(min(imSizeX+1,trNAonly(k).xCoord(iSF)+r1)); 
                    yminROI = round(max(1,trNAonly(k).yCoord(iSF)-(r1-1)));
                    ymaxROI = round(min(imSizeY+1,trNAonly(k).yCoord(iSF)+r1));

                    paxImageCropped2 = paxImageCropped(yminROI:ymaxROI,xminROI:xmaxROI);
                    invPaxImageCropped2 = imcomplement(paxImageCropped2);
                    lastPax = imread(strcat(paxPath,'/pax',num2str(nFrames,iiformat),'.tif'));
                    lastPaxCropped = lastPax(yminROI:ymaxROI,xminROI:xmaxROI);
                    invLastPaxCropped = imcomplement(lastPaxCropped);
                    pmax = max([invPaxImageCropped2(:); invLastPaxCropped(:)]);
                    pmin = min([invPaxImageCropped2(:); invLastPaxCropped(:)]);
                    if j~= iSF && strcmp(trNAonly(k).state{iSF} , 'NA') % remember the first NA's position
                        xFirst = trNAonly(k).xCoord(iSF);
                        yFirst = trNAonly(k).yCoord(iSF);
    %                     bkgAmpFirst = trNAonly(k).bkgAmp(iSF);
                        ampFirst = trNAonly(k).amp(iSF);
                        sigmaFirst = trNAonly(k).sigma(iSF);
                    end
                end
                if isempty(ha1)
                    imshow(imcomplement(paxImageCropped2),[pmin pmax]),hold on
                    ha1 = get(h2,'CurrentAxes');
                else
                    imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1); hold(ha1,'on')
                end
                if strcmp(trNAonly(k).state{j} , 'BA')
                    % drawing tracks
                    plot(ha1,xFirst-xminROI,yFirst-yminROI,'g', 'LineWidth', 0.5)
                    % remembering adhesion intensity and TF
                    ynmin = max(1,round(yFirst)-neighPix);
                    ynmax = min(size(tsMap,1),round(yFirst)+neighPix);
                    xnmin = max(1,round(xFirst)-neighPix);
                    xnmax = min(size(tsMap,2),round(xFirst)+neighPix);
                    %here, now I try to do gaussian fit that
                    %pointsourceDetection used
                    paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                    pstruct = fitGaussians2D(double(paxImageCropped), xFirst, yFirst, ampFirst*0.1, sigmaFirst*0.5, min(paxNeigh(:)),'xyAc');
                    trNAonly(k).amp(j) = pstruct.A(1);
                    if isnan(pstruct.A)
                        trNAonly(k).amp(j) = double(paxImageCropped(round(yFirst),round(xFirst)))-min(paxNeigh(:));
                    end
                elseif strcmp(trNAonly(k).state{j} , 'NA')
                    % drawing tracks
                    plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'r', 'LineWidth', 0.5)
                elseif strcmp(trNAonly(k).state{j} , 'FC')
                    % drawing tracks
                    plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                    adhBoundary = trNAonly(k).adhBoundary{j};
                    plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                elseif strcmp(trNAonly(k).state{j} , 'FA')
                    % drawing tracks
                    plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'k', 'LineWidth', 0.5)
                    adhBoundary = trNAonly(k).adhBoundary{j};
                    plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'k', 'LineWidth', 0.5) %adhesion boundary
                end
                curRenderer = get(h2,'Renderer');
                if ~strcmp(curRenderer,'painters')
                    set(h2,'Renderer','painters')
                end
                print('-depsc2', '-r150', strcat(eachEpsPath,'/trackFrame',num2str(j,iiformat),'.eps'));
                print('-dtiff', '-r150', strcat(eachPaxPath,'/trackFrame',num2str(j,iiformat),'.tif'));
                hold(ha1,'off')
            end
            seeMore = input('Do you want to look at video until the end of the movie? (y/(n))','s');
            if isempty(seeMore)
                seeMore = 'n';
            end
            if strcmp(seeMore,'y')
                if fend < nFrames
                    for j=fend+1:nFrames
                        %acquiring force
                        ynmin = max(1,round(trNAonly(k).yCoord(fend))-neighPix);
                        ynmax = min(size(tsMap,1),round(trNAonly(k).yCoord(fend))+neighPix);
                        xnmin = max(1,round(trNAonly(k).xCoord(fend))-neighPix);
                        xnmax = min(size(tsMap,2),round(trNAonly(k).xCoord(fend))+neighPix);

                        paxImageCropped = imread(strcat(paxPath,'/pax',num2str(j,iiformat),'.tif'));
                        xminROI = round(max(1,trNAonly(k).xCoord(iSF)-(r1-1))); 
                        xmaxROI = round(min(imSizeX+1,trNAonly(k).xCoord(iSF)+r1)); 
                        yminROI = round(max(1,trNAonly(k).yCoord(iSF)-(r1-1)));
                        ymaxROI = round(min(imSizeY+1,trNAonly(k).yCoord(iSF)+r1));

                        paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                        xLast = trNAonly(k).xCoord(j-1);
                        yLast = trNAonly(k).yCoord(j-1);
                        ampLast = trNAonly(k).amp(j-1);
                        sigmaLast = 2.1; %trNAonly(k).sigma(j-1);
                        pstruct = fitGaussians2D(double(paxImageCropped), xLast, yLast, ampLast, sigmaLast*0.5, min(paxNeigh(:)),'xyAc');
                        trNAonly(k).amp(j) = pstruct.A(1);
                        trNAonly(k).xCoord(j) = pstruct.x(1);
                        trNAonly(k).yCoord(j) = pstruct.y(1);
                        if isnan(pstruct.A)
                            trNAonly(k).amp(j) = double(paxImageCropped(round(yLast),round(xLast)))-min(paxNeigh(:));
                            trNAonly(k).xCoord(j) = trNAonly(k).xCoord(j-1);
                            trNAonly(k).yCoord(j) = trNAonly(k).yCoord(j-1);
                        end

                        paxImageCropped2 = paxImageCropped(yminROI:ymaxROI,xminROI:xmaxROI);
                        ha1 = get(h2,'CurrentAxes');%subplot('position',[0  0.5  1  0.5]);
                        if isempty(ha1)
                            imshow(imcomplement(paxImageCropped2),[pmin pmax]),hold on
                            ha1 = get(h2,'CurrentAxes');
                        else
                            imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1); hold(ha1,'on')
                        end
                        if strcmp(trNAonly(k).state{fend} , 'NA')
                            % drawing tracks
                            plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'r', 'LineWidth', 0.5)
                        elseif strcmp(trNAonly(k).state{fend} , 'FC')
                            % drawing tracks
                            plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                            adhBoundary = trNAonly(k).adhBoundary{j};
                            if ~isnan(adhBoundary)
                                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                            end
                        elseif strcmp(trNAonly(k).state{fend} , 'FA')
                            % drawing tracks
                            plot(ha1,trNAonly(k).xCoord(1:j)-xminROI,trNAonly(k).yCoord(1:j)-yminROI,'k', 'LineWidth', 0.5)
                            adhBoundary = trNAonly(k).adhBoundary{j};
                            if ~isnan(adhBoundary)
                                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                            end
                        end
                        curRenderer = get(h2,'Renderer');
                        if ~strcmp(curRenderer,'painters')
                            set(h2,'Renderer','painters')
                        end
                        print('-depsc2','-r150', strcat(eachEpsPath,'/trackFrame',num2str(j,iiformat),'.eps'));
                        print('-dtiff', '-r150', strcat(eachPaxPath,'/trackFrame',num2str(j,iiformat),'.tif'));
                        hold(ha1,'off')
                   end
                else
                    disp('It is already the end of movie!')
                end
            end
            % inquire if the estimated state is right
            disp(['The state of this adhesion is : track' num2str(k)])
            disp([num2cell(trNAonly(k).iFrame(fstart:fend)') trNAonly(k).state(fstart:fend)'])
            strFA = input('Is this state well describing what you see in the movie (no over-estimated FAs or noisy NAs)? ((y)/n)','s');
            while strcmp(strFA,'n')
                iFrames = input('Which frames do you want to change? n1:n2  ');
                state = input('What is the state in those range? (e.g. BA, NA, FC, FA)  ','s');
                for jj=iFrames
                    trNAonly(k).state{jj} = state;
                end
                disp(['Now, the state of this adhesion is :' num2str(k)])
                disp([num2cell(trNAonly(k).iFrame(fstart:fend)') trNAonly(k).state(fstart:fend)'])
                strFA = input('Is this state well describing what you see in the movie (no over-estimated FAs or noisy NAs)? ((y)/n)','s');
            end
            hold off
        end
    end
end
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
%     masGap=20;
%     for gap = masGap:-1:1;
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
%     end
%         if j<nFrames-9 && sum(newTracks(i).presence(j:j+9))==10 ...
%                 && sum(isnan(newTracks(i).xCoord(j:j+9)))==10 
%             gap = 10;
%             for kk=1:gap
%                 newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
%                 newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
%                 newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
%                 newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
%             end
%         elseif j<nFrames-8 && sum(newTracks(i).presence(j:j+8))==9 ...
%                 && sum(isnan(newTracks(i).xCoord(j:j+8)))==9 
%             gap = 9;
%             for kk=1:gap
%                 newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
%                 newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
%                 newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
%                 newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
%             end
%         elseif j<nFrames-7 && sum(newTracks(i).presence(j:j+7))==8 ...
%                 && sum(isnan(newTracks(i).xCoord(j:j+7)))==8 
%             gap = 8;
%             for kk=1:gap
%                 newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
%                 newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
%                 newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
%                 newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
%            end
%         elseif j<nFrames-6 && sum(newTracks(i).presence(j:j+6))==7 ...
%                 && sum(isnan(newTracks(i).xCoord(j:j+6)))==7 
%             gap = 7;
%             for kk=1:gap
%                 newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
%                 newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
%                 newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
%                 newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
%             end
%         elseif j<nFrames-5 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
%                && newTracks(i).presence(j+4) && newTracks(i).presence(j+5) && isnan(newTracks(i).xCoord(j)) ...
%                && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))...
%                && isnan(newTracks(i).xCoord(j+4)) && isnan(newTracks(i).xCoord(j+5))
%             newTracks(i).xCoord(j) = (6*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j) = (6*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j) = (6*newTracks(i).amp(j-1)+newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j) = (6*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+6))/7;
%             newTracks(i).xCoord(j+1) = (5*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j+1) = (5*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j+1) = (5*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j+1) = (5*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+6))/7;
%             newTracks(i).xCoord(j+2) = (4*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j+2) = (4*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j+2) = (4*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j+2) = (4*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+6))/7;
%             newTracks(i).xCoord(j+3) = (3*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j+3) = (3*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j+3) = (3*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j+3) = (3*newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+6))/7;
%             newTracks(i).xCoord(j+4) = (2*newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j+4) = (2*newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j+4) = (2*newTracks(i).amp(j-1)+5*newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j+4) = (2*newTracks(i).bkgAmp(j-1)+5*newTracks(i).bkgAmp(j+6))/7;
%             newTracks(i).xCoord(j+5) = (newTracks(i).xCoord(j-1)+6*newTracks(i).xCoord(j+6))/7;
%             newTracks(i).yCoord(j+5) = (newTracks(i).yCoord(j-1)+6*newTracks(i).yCoord(j+6))/7;
%             newTracks(i).amp(j+5) = (newTracks(i).amp(j-1)+6*newTracks(i).amp(j+6))/7;
%             newTracks(i).bkgAmp(j+5) = (newTracks(i).bkgAmp(j-1)+6*newTracks(i).bkgAmp(j+6))/7;
%         elseif j<nFrames-4 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
%                 && newTracks(i).presence(j+4) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) ...
%                 && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3)) && isnan(newTracks(i).xCoord(j+4))
%             newTracks(i).xCoord(j) = (5*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+5))/6;
%             newTracks(i).yCoord(j) = (5*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+5))/6;
%             newTracks(i).amp(j) = (5*newTracks(i).amp(j-1)+newTracks(i).amp(j+5))/6;
%             newTracks(i).bkgAmp(j) = (5*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+5))/6;
%             newTracks(i).xCoord(j+1) = (4*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+5))/6;
%             newTracks(i).yCoord(j+1) = (4*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+5))/6;
%             newTracks(i).amp(j+1) = (4*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+5))/6;
%             newTracks(i).bkgAmp(j+1) = (4*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+5))/6;
%             newTracks(i).xCoord(j+2) = (3*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+5))/6;
%             newTracks(i).yCoord(j+2) = (3*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+5))/6;
%             newTracks(i).amp(j+2) = (3*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+5))/6;
%             newTracks(i).bkgAmp(j+2) = (3*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+5))/6;
%             newTracks(i).xCoord(j+3) = (2*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+5))/6;
%             newTracks(i).yCoord(j+3) = (2*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+5))/6;
%             newTracks(i).amp(j+3) = (2*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+5))/6;
%             newTracks(i).bkgAmp(j+3) = (2*newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+5))/6;
%             newTracks(i).xCoord(j+4) = (newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+5))/6;
%             newTracks(i).yCoord(j+4) = (newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+5))/6;
%             newTracks(i).amp(j+4) = (newTracks(i).amp(j-1)+5*newTracks(i).amp(j+5))/6;
%             newTracks(i).bkgAmp(j+4) = (newTracks(i).bkgAmp(j-1)+5*newTracks(i).bkgAmp(j+5))/6;
%         elseif j<nFrames-3 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
%                 && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))
%             newTracks(i).xCoord(j) = (4*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+4))/5;
%             newTracks(i).yCoord(j) = (4*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+4))/5;
%             newTracks(i).amp(j) = (4*newTracks(i).amp(j-1)+newTracks(i).amp(j+4))/5;
%             newTracks(i).bkgAmp(j) = (4*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+4))/5;
%             newTracks(i).xCoord(j+1) = (3*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+4))/5;
%             newTracks(i).yCoord(j+1) = (3*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+4))/5;
%             newTracks(i).amp(j+1) = (3*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+4))/5;
%             newTracks(i).bkgAmp(j+1) = (3*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+4))/5;
%             newTracks(i).xCoord(j+2) = (2*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+4))/5;
%             newTracks(i).yCoord(j+2) = (2*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+4))/5;
%             newTracks(i).amp(j+2) = (2*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+4))/5;
%             newTracks(i).bkgAmp(j+2) = (2*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+4))/5;
%             newTracks(i).xCoord(j+3) = (newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+4))/5;
%             newTracks(i).yCoord(j+3) = (newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+4))/5;
%             newTracks(i).amp(j+3) = (newTracks(i).amp(j-1)+4*newTracks(i).amp(j+4))/5;
%             newTracks(i).bkgAmp(j+3) = (newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+4))/5;
%         elseif j<nFrames-2 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) ...
%                 && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2))
%             newTracks(i).xCoord(j) = (3*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+3))/4;
%             newTracks(i).yCoord(j) = (3*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+3))/4;
%             newTracks(i).amp(j) = (3*newTracks(i).amp(j-1)+newTracks(i).amp(j+3))/4;
%             newTracks(i).bkgAmp(j) = (3*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+3))/4;
%             newTracks(i).xCoord(j+1) = (2*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+3))/4;
%             newTracks(i).yCoord(j+1) = (2*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+3))/4;
%             newTracks(i).amp(j+1) = (2*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+3))/4;
%             newTracks(i).bkgAmp(j+1) = (2*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+3))/4;
%             newTracks(i).xCoord(j+2) = (newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+3))/4;
%             newTracks(i).yCoord(j+2) = (newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+3))/4;
%             newTracks(i).amp(j+2) = (newTracks(i).amp(j-1)+3*newTracks(i).amp(j+3))/4;
%             newTracks(i).bkgAmp(j+2) = (newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+3))/4;
%         elseif j<nFrames-1 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1))
%             newTracks(i).xCoord(j) = (2*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+2))/3;
%             newTracks(i).yCoord(j) = (2*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+2))/3;
%             newTracks(i).amp(j) = (2*newTracks(i).amp(j-1)+newTracks(i).amp(j+2))/3;
%             newTracks(i).bkgAmp(j) = (2*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+2))/3;
%             newTracks(i).xCoord(j+1) = (newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+2))/3;
%             newTracks(i).yCoord(j+1) = (newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+2))/3;
%             newTracks(i).amp(j+1) = (newTracks(i).amp(j-1)+2*newTracks(i).amp(j+2))/3;
%             newTracks(i).bkgAmp(j+1) = (newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+2))/3;
%         elseif newTracks(i).presence(j) && isnan(newTracks(i).xCoord(j))
%             newTracks(i).xCoord(j) = (newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+1))/2;
%             newTracks(i).yCoord(j) = (newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+1))/2;
%             newTracks(i).amp(j) = (newTracks(i).amp(j-1)+newTracks(i).amp(j+1))/2;
%             newTracks(i).bkgAmp(j) = (newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+1))/2;
%         end
%     end
    if isempty(newTracks(i).startingFrame)
        warning(['startingFrame is empty for track #' num2str(i)])
    end
end
end

% function pstruct_NAwithForce = findMagCurvature(tsMap,pstruct_NA,neighD)
%     nPoints = length(pstruct_NA.x);
%     pstruct_NAwithForce = pstruct_NA;
%     laplacian = [.5 1 .5; 1 -6 1; .5 1 .5];
%     for jj=1:nPoints
%         rowRange = round(pstruct_NA.y(jj))-neighD:round(pstruct_NA.y(jj))+neighD;
%         colRange = round(pstruct_NA.x(jj))-neighD:round(pstruct_NA.x(jj))+neighD;
%         pstruct_NAwithForce.fmag(jj) = max(max(tsMap(rowRange,colRange))); %force magnitude
%         pstruct_NAwithForce.fcurvature(jj) = sum(sum(tsMap(round(pstruct_NA.y(jj))-1:round(pstruct_NA.y(jj))+1,round(pstruct_NA.x(jj))-1:round(pstruct_NA.x(jj))+1)...
%                                                                             .* laplacian)); %force curvature
%     end
% end