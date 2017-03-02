function analyzeAdhesionMaturation(MD)
% function [trNAonly,indFail,indMature,lifeTimeNAfailing,lifeTimeNAmaturing,maturingRatio,NADensity,FADensity,focalAdhInfo] = analyzeAdhesionMaturation(pathForTheMovieDataFile, varargin)
% [tracksNA,lifeTimeNA] = analyzeAdhesionMaturation(pathForTheMovieDataFile,outputPath,showAllTracks,plotEachTrack)
% filter out NA tracks, obtain life time of each NA tracks
% input:    pathForTheMovieDataFile:       path to the movieData file (FA
%                   segmentation and NA    tracking package should be run beforehand)
%           outputPath                      outputPath
%           'onlyEdge' [false]              collect NA tracks that ever close to cell edge
%           'matchWithFA' [true]            For cells with only NAs, we turn this off.
%           'getEdgeRelatedFeatures'[true]  For cells with only NAs, we turn this off.
%           'reTrack' [true]     This is for 
%           'minLifetime' [5]               For cells with only NAs, we turn this off.
%           'iChan' [1]                        Channel with FA marker
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
% Andrew R. Jamieson Feb. 2017 - Updating to incorporate into MovieData Process GUI (Focal Adhesion Package)

%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('MD', @(x)(isa(x,'MovieData')));
ip.parse(MD);
%Get the indices of any previous processes
iProc = MD.getProcessIndex('AdhesionAnalysisProcess', 'nDesired', 1, 'askUser', false);
%Check if process exists
if isempty(iProc)
    error('No AdhesionAnalysisProcess in input movieData! please create the process and use the process.run method to run this function!')
end
%Parse input, store in parameter structure
thisProc = MD.processes_{iProc};
p = parseProcessParams(thisProc);

%% --------------- Parameters ---------- %%
matchWithFA = p.matchWithFA;
minLifetime = p.minLifetime;
onlyEdge = p.onlyEdge;
reTrack = p.reTrack;
getEdgeRelatedFeatures = p.getEdgeRelatedFeatures;
iChan = p.ChannelIndex;
bandwidthNA = p.bandwidthNA;

ApplyCellSegMask = p.ApplyCellSegMask;

% Load Respective Process objects
if ApplyCellSegMask
    maskProc = MD.getProcess(p.SegCellMaskProc);
end
detectedNAProc = MD.getProcess(p.detectedNAProc);
trackNAProc = MD.getProcess(p.trackFAProc);
FASegProc = MD.getProcess(p.FAsegProc);

%% ------------------ Config Output  ---------------- %%

% Set up the output file
%% Backup previous Analysis output
if exist(p.OutputDirectory,'dir')
    if p.backupOldResults
        disp('Backing up the original data')
        ii = 1;
        backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        if strcmp(computer('arch'), 'win64')
            mkdir(backupFolder);
        else
            try
                system(['mkdir -p ' backupFolder]);
            catch
                mkdir(backupFolder);
            end            
        end
        copyfile(p.OutputDirectory, backupFolder,'f')
    end
end
mkClrDir(p.OutputDirectory);

%% TODO -- What is the "canonical" way to define this for multi-input situations?
% Set up the input directories (input images)
inFilePaths = cell(4,numel(MD.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = MD.getChannelPaths(i);
    inFilePaths{2,i} = detectedNAProc.outFilePaths_{1,i};
    inFilePaths{3,i} = trackNAProc.outFilePaths_{1,i};
    inFilePaths{4,i} = FASegProc.outFilePaths_{1,i};
end
thisProc.setInFilePaths(inFilePaths);
    
% Set up the output files
outFilePaths = cell(4, numel(MD.channels_));
for i = p.ChannelIndex
    [abpath chanDirName ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_Chan' num2str(i) '_tracksNA'];
    outFilePaths{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
    dataPath_tracksNA = outFilePaths{1,i};

    outFilename = [chanDirName '_Chan' num2str(i) '_focalAdhInfo'];
    outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
    dataPath_focalAdhInfo = outFilePaths{2,i};

    outFilename = [chanDirName '_Chan' num2str(i) '_NAFADensity'];
    outFilePaths{3,i} = [p.OutputDirectory filesep outFilename '.mat'];
    dataPath_NAFADensity = outFilePaths{3,i};

    outFilename = [chanDirName '_Chan' num2str(i) '_allAnalysisFA'];
    outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '.mat'];
    dataPath_analysisAll = outFilePaths{4,i};
end

thisProc.setOutFilePaths(outFilePaths);

% Get whole frame number
nFrames = MD.nFrames_;
iPaxChannel = iChan;
%% TODO - Check with Sanity
iiformat = ['%.' '3' 'd'];
%% TODO - Check with Sangyoon on pixel size/minsize
% minSize = round((500/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=.5um*.5um
minLifetime = min(nFrames, minLifetime);
markerSize = 2;



%% For Debugging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadSaved_state = false;
if loadSaved_state
    [abpath outFilename ext] = fileparts(dataPath_tracksNA);
    load([backupFolder filesep outFilename ext]);
    [abpath outFilename ext] = fileparts(dataPath_focalAdhInfo);
    load([backupFolder filesep outFilename ext]);
    [abpath outFilename ext] = fileparts(dataPath_NAFADensity);
    load([backupFolder filesep outFilename ext]);
    [abpath outFilename ext] = fileparts(dataPath_analysisAll);
    load([backupFolder filesep outFilename ext]); 
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% tracks
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
tracksNA = formatTracks(tracksNAorg, detectedNAs, nFrames); 
toc

bandArea = zeros(nFrames,1);
NADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
FADensity = zeros(nFrames,1); % unit: number/um2 = numel(tracksNA)/(bandArea*MD.pixelSize^2*1e6);
numFCs = zeros(nFrames,1);
nChannels = numel(MD.channels_);
% minEcc = 0.7;


% Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered
if onlyEdge
    disp(['Filtering adhesion tracks based on cell mask. Only adhesions appearing at the edge are considered'])
    trackIdx = false(numel(tracksNA),1);
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    for ii=1:nFrames
        % Cell Boundary Mask 
        if ApplyCellSegMask
            mask = maskProc.loadChannelOutput(iChan,ii);
        else
            mask = true(MD.imSize_);
        end
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
    if ApplyCellSegMask
        mask = maskProc.loadChannelOutput(iChan,1);
    else
        mask = true(MD.imSize_);
    end
    for ii=1:nFrames
        % Cell Boundary Mask 
        if ApplyCellSegMask
            mask = maskProc.loadChannelOutput(iChan,ii);
        else
            mask = true(MD.imSize_);
        end
        % mask for band from edge
        maskOnlyBand = mask;
        bandArea(ii) = sum(maskOnlyBand(:)); % in pixel
        % filter tracks with naMasks
        % only deal with presence and status
        % Tracks in its emerging state ever overlap with bandMask are
        % considered.
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
                    ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
                    round(tracksNA(k).xCoord(ii)) < 1 || ...
                    round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
                    round(tracksNA(k).yCoord(ii)) < 1) || ...
                    ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
                tracksNA(k).state{ii} = 'Out_of_Band';
                tracksNA(k).presence(ii) = false;
                if trackIdx(k)
                    trackIdx(k) = false;
                end
            end
        end
    end
end
% get rid of tracks that have out of bands...
tracksNA = tracksNA(trackIdx);

    %% add intensity of tracks including before and after NA status
% get the movie stack
disp('Loading image stacks ...'); tic;
h = MD.imSize_(1); 
w = MD.imSize_(2);
imgStack = zeros(h,w,nFrames);
for ii = 1:nFrames
    imgStack(:,:,ii) = MD.channels_(iPaxChannel).loadImage(ii); 
end
toc;

% get the intensity
disp('Reading intensities with additional tracking...')
tic
% reTrack=false;
tracksNA = readIntensityFromTracks(tracksNA, imgStack, 1, 'extraLength',30,'movieData',MD,'retrack',reTrack); % 1 means intensity collection from pax image
toc
%% Filtering again after re-reading
disp('Filtering again after re-reading with cell mask ...')
tic
if onlyEdge
    trackIdx = true(numel(tracksNA),1);
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    for ii=1:nFrames
        % Cell Boundary Mask 
        if ApplyCellSegMask
            mask = maskProc.loadChannelOutput(iChan,ii);
        else
            mask = true(MD.imSize_);
        end
        % mask for band from edge
        iMask = imcomplement(mask);
        distFromEdge = bwdist(iMask);
        bandMask = distFromEdge <= bandwidthNA_pix;

        maskOnlyBand = bandMask & mask;
        bandArea(ii) = sum(maskOnlyBand(:)); % in pixel

        % collect index of tracks which first appear at frame ii
%         idxFirstAppear = arrayfun(@(x) x.startingFrameExtra==ii,tracksNA);
        % now see if these tracks ever in the maskOnlyBand
        for k=find(trackIdx)'
%         for k=find(idxFirstAppear)'
            if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
                    ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
                    round(tracksNA(k).xCoord(ii)) < 1 || ...
                    round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
                    round(tracksNA(k).yCoord(ii)) < 1) || ...
                    ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
                trackIdx(k) = false;
            end
        end
    end    
else
    disp('Entire adhesion tracks are considered.')
    trackIdx = true(numel(tracksNA),1);
    for ii=1:nFrames
        % Cell Boundary Mask 
        if ApplyCellSegMask
            mask = maskProc.loadChannelOutput(iChan,ii);
        else
            mask = true(MD.imSize_);
        end
        % mask for band from edge
        maskOnlyBand = mask;
        bandArea(ii) = sum(maskOnlyBand(:)); % in pixel
        % filter tracks with naMasks
        % only deal with presence and status
        % Tracks in its emerging state ever overlap with bandMask are
        % considered.
        for k=find(trackIdx)'
            if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
                    ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
                    round(tracksNA(k).xCoord(ii)) < 1 || ...
                    round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
                    round(tracksNA(k).yCoord(ii)) < 1) || ...
                    ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
                tracksNA(k).state{ii} = 'Out_of_Band';
                tracksNA(k).presence(ii) = false;
                if trackIdx(k)
                    trackIdx(k) = false;
                end
            end
        end
    end
end
toc
% get rid of tracks that have out of bands...
tracksNA = tracksNA(trackIdx);
%%%%%

%% Matching with adhesion setup
if ApplyCellSegMask
    firstMask=maskProc.loadChannelOutput(iChan,1);
else
    firstMask=true(MD.imSize_);
end

cropMaskStack = false(size(firstMask,1),size(firstMask,2),nFrames);
numTracks=numel(tracksNA);
progressText(0,'Matching with segmented adhesions:');
focalAdhInfo(nFrames,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'length',[],'meanFAarea',[],'medianFAarea',[]...
    ,'meanLength',[],'medianLength',[],'numberFA',[],'FAdensity',[],'cellArea',[],'ecc',[]);
FAInfo(nFrames,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'length',[]);

%% Matching with segmented adhesions
prevMask=[];
neighPix = 2;
for ii=1:nFrames
    % Cell Boundary Mask 
    if ApplyCellSegMask
        mask = maskProc.loadChannelOutput(iChan,ii);
        if ii>1 && max(mask(:))==0
            mask=prevMask;
            disp('Previous mask is used for cell edge because the current mask is empty.')
        else
            prevMask=mask;
        end
    else
        mask=true(MD.imSize_);
    end
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');
    cropMaskStack(:,:,ii) = mask;
    % Get the mask for FAs
    I=MD.channels_(iPaxChannel).loadImage(ii); 
    maskFAs = FASegProc.loadChannelOutput(iPaxChannel,ii);
    maskAdhesion = maskFAs>0 & mask;
    % FA Segmentation usually over-segments things. Need to chop them off
    % to smaller ones or filter insignificant segmentation out.
    xNA=arrayfun(@(x) x.xCoord(ii),tracksNA);
    yNA=arrayfun(@(x) x.yCoord(ii),tracksNA);
    maskAdhesion = refineAdhesionSegmentation(maskAdhesion,I,xNA,yNA,mask);
    % close once and dilate once
    % maskAdhesion = bwmorph(maskAdhesion,'close');
%         maskAdhesion = bwmorph(maskAdhesion,'thin',1);
%         Adhs = regionprops(maskAdhesion,'Area','Eccentricity','PixelIdxList','PixelList' );
    % Save focal adhesion information
    Adhs = regionprops(bwconncomp(maskAdhesion,4),'Centroid','Area','Eccentricity','PixelList','PixelIdxList','MajorAxisLength');
    numAdhs(ii) = numel(Adhs);
%         minFASize = round((1000/MD.pixelSize_)*(1000/MD.pixelSize_)); %adhesion limit=1um*1um
%         minFCSize = round((600/MD.pixelSize_)*(400/MD.pixelSize_)); %adhesion limit=0.6um*0.4um
    minFALength = round((2000/MD.pixelSize_)); %adhesion limit=2um
    minFCLength = round((600/MD.pixelSize_)); %adhesion limit=0.6um

%         fcIdx = arrayfun(@(x) x.Area<minFASize & x.Area>minFCSize, Adhs);
    fcIdx = arrayfun(@(x) x.MajorAxisLength<minFALength & x.MajorAxisLength>minFCLength, Adhs);
    FCIdx = find(fcIdx);
    adhBound = bwboundaries(maskAdhesion,4,'noholes');    

    % for larger adhesions
%         faIdx = arrayfun(@(x) x.Area>=minFASize, Adhs);
    faIdx = arrayfun(@(x) x.MajorAxisLength>=minFALength, Adhs);
    FAIdx =  find(faIdx);

    FCs = Adhs(fcIdx | faIdx);
    numFCs = length(FCs);
    for k=1:numFCs
        focalAdhInfo(ii).xCoord(k) = round(FCs(k).Centroid(1));
        focalAdhInfo(ii).yCoord(k) = round(FCs(k).Centroid(2));
        focalAdhInfo(ii).area(k) = FCs(k).Area;
        focalAdhInfo(ii).length(k) = FCs(k).MajorAxisLength;
        focalAdhInfo(ii).amp(k) = mean(I(FCs(k).PixelIdxList));
        focalAdhInfo(ii).ecc(k) = FCs(k).Eccentricity;
    end
    focalAdhInfo(ii).numberFA = numFCs;
    focalAdhInfo(ii).meanFAarea = mean(focalAdhInfo(ii).area);
    focalAdhInfo(ii).medianFAarea = median(focalAdhInfo(ii).area);
    focalAdhInfo(ii).meanLength = mean(focalAdhInfo(ii).length);
    focalAdhInfo(ii).medianLength = median(focalAdhInfo(ii).length);
    focalAdhInfo(ii).cellArea = sum(mask(:))*(MD.pixelSize_/1000)^2; % in um^2
    focalAdhInfo(ii).FAdensity = numFCs/focalAdhInfo(ii).cellArea; % number per um2
    focalAdhInfo(ii).numberFC = sum(fcIdx);
    focalAdhInfo(ii).FAtoFCratio = sum(faIdx)/sum(fcIdx);
    focalAdhInfo(ii).numberPureFA = sum(faIdx);

    FAs = Adhs(faIdx);
    numFAs = length(FAs);
    for k=1:numFAs
        FAInfo(ii).xCoord(k) = round(FAs(k).Centroid(1));
        FAInfo(ii).yCoord(k) = round(FAs(k).Centroid(2));
        FAInfo(ii).area(k) = FAs(k).Area;
        FAInfo(ii).length(k) = FAs(k).MajorAxisLength;
        FAInfo(ii).amp(k) = mean(I(FAs(k).PixelIdxList));
        FAInfo(ii).ecc(k) = FAs(k).Eccentricity;
    end

    if matchWithFA
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
                    [~, closestPixelID] = min(sqrt((propSubMaskFAs(subMaskFAsIdx).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame))).^2 +...
                        (propSubMaskFAs(subMaskFAsIdx).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame))).^2));

                    tracksNA(k).state{ii} = 'FC';
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
    if getEdgeRelatedFeatures
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
    end
    progressText(ii/(nFrames-1),'Matching with segmented adhesions:');
end

%% protrusion/retraction information
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.

if getEdgeRelatedFeatures
    for ii=1:numTracks
        idxZeros = tracksNA(ii).closestBdPoint(:)==0;
        tracksNA(ii).closestBdPoint(idxZeros)=NaN;
    end
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
%     presIdx = logical(tracksNA(k).presence);
    % get the instantaneous velocity
    % get the distance first
%     curTrackVelLength=sum(presIdx)-1;
%     presIdxSeq = find(presIdx);
    presIdxSeq = tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra;
    curTrackVelLength=length(presIdxSeq)-1;
    distTrajec=zeros(curTrackVelLength,1);
    if getEdgeRelatedFeatures
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

    % Assembly rate: Slope from emergence to maximum - added 10/27/2016 for
    % Michelle's analysis % This output might contain an error when the
    % value has noisy maximum. So it should be filtered with
    % earlyAmpSlope..
    splineParam=0.01; 
    tRange = tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra;
    curAmpTotal =  tracksNA(k).ampTotal(tRange);
    sd_spline= csaps(tRange,curAmpTotal,splineParam);
    sd=ppval(sd_spline,tRange);
    % Find the maximum
    [~,maxSdInd] = max(sd);
    maxAmpFrame = tRange(maxSdInd);
    [~,assemRate] = regression(tIntervalMin*(tRange(1:maxSdInd)),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:maxAmpFrame));
    tracksNA(k).assemRate = assemRate; % in a.u./min
    
    % Disassembly rate: Slope from maximum to end
    [~,disassemRate] = regression(tIntervalMin*(tRange(maxSdInd:end)),tracksNA(k).ampTotal(maxAmpFrame:tracksNA(k).endingFrameExtra));
    tracksNA(k).disassemRate = disassemRate; % in a.u./min

    curStartFrame = max(tracksNA(k).startingFrame,tracksNA(k).endingFrameExtra-periodFrames+1);
    curLatePeriod = tracksNA(k).endingFrameExtra - curStartFrame+1;
    [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksNA(k).ampTotal(curStartFrame:tracksNA(k).endingFrameExtra));
    tracksNA(k).lateAmpSlope = curMlate; % in a.u./min

    curEndFrame = min(sF+periodFrames-1,eF);
    curEarlyPeriod = curEndFrame - sF+1;
    if getEdgeRelatedFeatures
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
        try
            fitobj = fit(tracksNA(k).xCoord(sF:eF)',tracksNA(k).yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
        catch
            tracksNA(k)=readIntensityFromTracks(tracksNA(k),imgStack,1,'extraLength',30,'movieData',MD,'retrack',reTrack);
        end
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
    end
    %Mean Squared Displacement
    meanX = mean(tracksNA(k).xCoord(logical(tracksNA(k).presence)));
    meanY = mean(tracksNA(k).yCoord(logical(tracksNA(k).presence)));
    tracksNA(k).MSD=sum((tracksNA(k).xCoord(logical(tracksNA(k).presence))'-meanX).^2+...
        (tracksNA(k).yCoord(logical(tracksNA(k).presence))'-meanY).^2);
    tracksNA(k).MSDrate = tracksNA(k).MSD/tracksNA(k).lifeTime;
    progressText(k/(numTracks-1),'Post-analysis:');
end
end
%% saving
tableTracksNA = struct2table(tracksNA);
save(dataPath_tracksNA, 'tracksNA', 'tableTracksNA');
save(dataPath_focalAdhInfo, 'focalAdhInfo')

%% NA FA Density analysis
numNAs = zeros(nFrames,1);
numNAsInBand = zeros(nFrames,1);
trackIdx = true(numel(tracksNA),1);
bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
for ii=1:nFrames
    mask = maskProc.loadChannelOutput(iChan,ii);
    % mask for band from edge
    iMask = imcomplement(mask);
    distFromEdge = bwdist(iMask);
    bandMask = distFromEdge <= bandwidthNA_pix;

    maskOnlyBand = bandMask & mask;
    bandArea(ii) = sum(maskOnlyBand(:)); % in pixel

    % collect index of tracks which first appear at frame ii
%         idxFirstAppear = arrayfun(@(x) x.startingFrameExtra==ii,tracksNA);
    % now see if these tracks ever in the maskOnlyBand
    for k=find(trackIdx)'
%         for k=find(idxFirstAppear)'
        if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
                ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
                round(tracksNA(k).xCoord(ii)) < 1 || ...
                round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
                round(tracksNA(k).yCoord(ii)) < 1) || ...
                ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
            trackIdx(k) = false;
        end
    end
%         for k=find(trackIdxFC)' % I'll work on this later (1/10/17) SH
% %         for k=find(idxFirstAppear)'
%             if focalAdhInfo(ii).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
%                     ((round(tracksNA(k).xCoord(ii)) > size(maskOnlyBand,2) || ...
%                     round(tracksNA(k).xCoord(ii)) < 1 || ...
%                     round(tracksNA(k).yCoord(ii)) > size(maskOnlyBand,1) || ...
%                     round(tracksNA(k).yCoord(ii)) < 1) || ...
%                     ~maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii))))
%                 trackIdx(k) = false;
%             end
%         end
    numNAs(ii) = sum(arrayfun(@(x) (x.presence(ii)==true), tracksNA));
    numNAsInBand(ii) = sum(trackIdx);
    NADensity(ii) = numNAsInBand(ii)/(bandArea(ii)*MD.pixelSize_^2*1e-6);  % unit: number/um2 
    FADensity(ii) = numNAsInBand(ii)/(bandArea(ii)*MD.pixelSize_^2*1e-6);  % unit: number/um2 
end
save(dataPath_NAFADensity, 'NADensity','FADensity','bandwidthNA','numNAsInBand')

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
if matchWithFA
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
    save(dataPath_analysisAll);%, 'trNAonly', 'tracksNAfailing','tracksNAmaturing','maturingRatio','lifeTimeNAfailing','lifeTimeNAmaturing')
else
    trNAonly = tracksNA;
    indMature = [];
    indFail = [];
    lifeTimeNAfailing=[];
    lifeTimeNAmaturing =[];
    maturingRatio = [];
end
end

function newTracks = formatTracks(tracks,detectedNAs,nFrames)
% Format tracks structure into tracks with every frame

newTracks(numel(tracks),1) = struct('xCoord', [], 'yCoord', [],'state',[],'iFrame',[],'presence',[],'amp',[],'bkgAmp',[]);
% BA: before adhesion, NA: nascent adh, FC: focal complex, FA: focal adh,
% ANA: after NA (failed to be matured.
%OUTPUT tracksFinal   : Structure array where each element corresponds to a 
%                       compound track. Each element contains the following 
%                       fields:
%           .tracksFeatIndxCG: Connectivity matrix of features between
%                              frames, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = number of frames
%                              the compound track spans. Zeros indicate
%                              frames where track segments do not exist
%                              (either because those frames are before the
%                              segment starts or after it ends, or because
%                              of losing parts of a segment.
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist, like the zeros above.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a compound track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 = start of track segment, 2 = end of track segment;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN = start is a birth and end is a death,
%                                   number = start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
for i = 1:numel(tracks)
    % Get the x and y coordinate of all compound tracks
    startNA = true;
    endNA = true;
    for  jj = 1 : nFrames
        newTracks(i).iFrame(jj) = jj;
        if jj<tracks(i).seqOfEvents(1,1)     % before Adhesion
            newTracks(i).state{jj} = 'BA';
            newTracks(i).xCoord(jj) = NaN;
            newTracks(i).yCoord(jj) = NaN;
            newTracks(i).presence(jj) = false;
            newTracks(i).amp(jj) = NaN;
            newTracks(i).bkgAmp(jj) = NaN;
        elseif jj>tracks(i).seqOfEvents(end,1) % Adhesion gone!
            newTracks(i).state{jj} = 'ANA';
            newTracks(i).xCoord(jj) = NaN;
            newTracks(i).yCoord(jj) = NaN;
            newTracks(i).amp(jj) = NaN;
            newTracks(i).bkgAmp(jj) = NaN;
            newTracks(i).presence(jj) = false;
            if endNA
                newTracks(i).endingFrame = jj-1; % mark end of frame
                endNA = false;
            end
        elseif jj==tracks(i).seqOfEvents(2,1) % Adhesion occuring
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
                newTracks(i).startingFrame = jj; % Frame when adhesion starts
                startNA = false;
            end
            if endNA
                newTracks(i).endingFrame = jj; % Frame when adhesion ends
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
                newTracks(i).bkgAmp(jj) = detectedNAs(jj).bkg(tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(jj) = detectedNAs(jj).sigmaX(tracks(i).tracksFeatIndxCG(jj-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(jj) = true;
            if startNA
                newTracks(i).startingFrame = jj;
                startNA = false;
            end
        end
            
        if isfield(tracks, 'label')
            %% TODO -- what is this?
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

    if isempty(newTracks(i).startingFrame)
        warning(['startingFrame is empty for track #' num2str(i)])
    end
end
end