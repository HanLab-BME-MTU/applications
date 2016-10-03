function [tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,band,tmax,showAllTracks,plotEachTrack,tmaxEach,varargin)
% colocalizationTFMwithFA runs colocalization between peaks in TFM maps and peaks in paxillin using movieData.
% Basically this function make grayscale of TFM and pick up peaks of the
% map and see if how many of them are colocalized with significant paxillin
% signal or neighborhood of nascent adhesions.

% input:    pathForTheMovieDataFile:    path to the movieData file (TFM
% package should be run beforehand)
%           outputPath              outputPath
%           band:                       band width for cutting border
%           (default=4)
%           pointTF:                   true if you want to trace force at
%           a certain point (default = false)
% output:   images of heatmap stored in pathForTheMovieDataFile/heatmap
%           forceNAdetec,          force at detectable NAs
%           forceNAund,             force at undetectable NAs
%           forceFA,                     force magnitude at FA
%           peaknessRatio,      portion of forces that are themselves local maxima
%           DTMS                        Deviation of Traction Magnitude from Surrounding
% Note: detectable NAs are defined by coincidence of NA position with bead
% location and force node
% Sangyoon Han April 2013

% if nargin < 3
%     band = 16;
%     tmax=[];
%     plotEachTrack = false;
%     showAllTracks = true;
% elseif nargin <4
%     tmax=[];
%     plotEachTrack = false;
%     showAllTracks = true;
% elseif nargin<5
%     plotEachTrack = false;
%     showAllTracks = true;
% elseif nargin<6
%     plotEachTrack = false;
% end
% if isempty(pstruct_NAcell)
%     bFirst=true;
% else
%     bFirst=false;
% end
% bandwidth = 20; %um
ip =inputParser;
ip.addRequired('pathForTheMovieDataFile',@ischar)
ip.addOptional('band',4,@isscalar)
ip.addOptional('tmax',[],@(x)isscalar(x)||isempty(x))
ip.addOptional('showAllTracks',false,@(x)islogical(x)||isempty(x))
ip.addOptional('plotEachTrack',false,@(x)islogical(x)||isempty(x))
ip.addOptional('tmaxEach',[],@(x)isscalar(x)||isempty(x))
ip.addParamValue('onlyEdge',false,@islogical); % collect NA tracks that ever close to cell edge
ip.addParamValue('outputPath',[],@ischar)
ip.addParamValue('saveAnalysis',true,@islogical)
ip.addParamValue('matchWithFA',true,@islogical) %For cells with only NAs, we turn this off.
ip.addParamValue('minLifetime',3,@isscalar) %For cells with only NAs, we turn this off.
ip.addParamValue('showTrackID',false,@islogical) %For cells with only NAs, we turn this off.
ip.parse(pathForTheMovieDataFile,band,tmax,showAllTracks,plotEachTrack,tmaxEach,varargin{:});
outputPath=ip.Results.outputPath;
saveAnalysis=ip.Results.saveAnalysis;
matchWithFA=ip.Results.matchWithFA;
minLifetime=ip.Results.minLifetime;
% showAllTracks=ip.Results.showAllTracks;
% plotEachTrack=ip.Results.plotEachTrack;
onlyEdge=ip.Results.onlyEdge;
showTrackID = ip.Results.showTrackID;

%% Data Set up
% Load the MovieData
[pathStr,name,ext]= fileparts(pathForTheMovieDataFile);
if strcmp([name ext],'movieData.mat')
    movieDataPath = pathForTheMovieDataFile;
    pathForTheMovieDataFile = pathStr;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
end
MD = MovieData.load(movieDataPath);
% Get whole frame number
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Get FA segmentation package 
iFASegPackage = MD.getPackageIndex('FocalAdhesionSegmentationPackage');
% If not run, run the package
iFASeg = 6;
if isempty(iFASegPackage)
    MD.addPackage(FocalAdhesionSegmentationPackage(MD));
    iPack =  MD.getPackageIndex('FocalAdhesionSegmentationPackage');
    MD.getPackage(iPack).createDefaultProcess(iFASeg)
    params = MD.getPackage(iPack).getProcess(iFASeg).funParams_;
    params.ChannelIndex = 2; %paxillin
    params.SteerableFilterSigma = 72; % in nm
    params.OpeningRadiusXY = 0; % in nm
    params.MinVolTime = 1; %um2*s
    params.OpeningHeightT = 10; % sec
    MD.getPackage(iPack).getProcess(iFASeg).setPara(params);
    MD.getPackage(iPack).getProcess(iFASeg).run();
    MD.save;
    iFASegPackage = iPack;
end
FASegPackage=MD.getPackage(iFASegPackage);
% Load tracks
% iUTrack = MD.getPackageIndex('UTrackPackage');
iTrackProc = MD.getProcessIndex('TrackingProcess');%MD.getPackageIndex('UTrackPackage');
% If not run, run the package
trackNAProc = MD.getProcess(iTrackProc);

forceFC = zeros(nFrames,1);
forceFA = zeros(nFrames,1);
forceBGband(nFrames,1) = struct('mean',[],'err',[]);
forceBG(nFrames,1) = struct('mean',[],'err',[]);

% tracksNA = trackNAProc.loadChannelOutput;
% Load the displField
iDispFieldProc = 3;
displFieldProc=TFMPackage.processes_{iDispFieldProc};
maskArray = MD.getROIMask;
% Use mask of first frame to filter displacementfield
firstMask = maskArray(:,:,1);
displFieldOriginal=displFieldProc.loadChannelOutput;
displField = filterDisplacementField(displFieldOriginal,firstMask);

% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceField=forceFieldProc.loadChannelOutput;

%filter out force temporal noise 
% forceField=filterForceShortPeaks(forceField,60);

% Load Cell Segmentation
iMask = MD.getProcessIndex('MaskRefinementProcess');
if isempty(iMask)
    iMask = MD.getProcessIndex('ThresholdProcess');
end
maskProc = MD.getProcess(iMask);

% Set up the output file path
outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
paxPath = [outputFilePath filesep 'pax'];
forcemapPath = [outputFilePath filesep 'fMap'];
forcetifPath = [outputFilePath filesep 'forcetifs'];
paxtifPath = [outputFilePath filesep 'paxtifs'];
epsPath = [outputFilePath filesep 'eps'];
figPath = [outputFilePath filesep 'figs'];
if ~exist(forcetifPath,'dir') || ~exist(paxtifPath,'dir') || ~exist(paxPath,'dir') || ~exist(figPath,'dir') || ~exist(epsPath,'dir') || ~exist(forcemapPath,'dir')  || ~exist(dataPath,'dir') 
    mkdir(paxPath);
    mkdir(forcetifPath);
    mkdir(paxtifPath);
    mkdir(figPath);
    mkdir(epsPath);
    mkdir(forcemapPath);
    mkdir(dataPath);
end
%     set(h2, 'Position', [100+cropInfo(3)-cropInfo(1)*10/9 100 cropInfo(3)-cropInfo(1) cropInfo(4)-cropInfo(2)])
iiformat = ['%.' '3' 'd'];
minLifetime = min(nFrames,minLifetime/MD.timeInterval_);
markerSize = 4;
% tracks
iPaxChannel = 2; % this should be intentionally done in the analysis level
if ~trackNAProc.checkChannelOutput(iPaxChannel)
    iPaxChannel = 1;
end

% % Find out force at each tracksNA at this frame
% % find first the index of relevant tracks from seqOfEvents
% % tracksNA = trackNAProc.loadChannelOutput(iPaxChannel,ii);
% disp('Loading NA tracks...')
% tic
% tracksNAorg = trackNAProc.loadChannelOutput(iPaxChannel);
% toc
% % if ii>minLifetime
% SEL = getTrackSEL(tracksNAorg); %SEL: StartEndLifetime
% % Remove any less than 3-frame long track.
% isValid = SEL(:,3) >= minLifetime;
% tracksNAorg = tracksNAorg(isValid);
% end
% See if there is stage drift correction
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
else
    iBeadChan=1;
    iChan = 2;
end
outputFile=strcat(dataPath,filesep,'tracksNA.mat');
% See if you can use existing tracks
if exist(outputFile,'file')
    disp('tracksNA file is found. Using it ... If you do not want to reuse this, please backup the file and rename it to other name than tracksNA.mat.')
    tracksNAFile = load(outputFile,'tracksNA');
    tracksNA = tracksNAFile.tracksNA;
else
    % run analyzeAdhesionMaturation for obtaining tracks from paxillin channel
    tracksNA = analyzeAdhesionMaturation(pathForTheMovieDataFile,false,false,...
        'onlyEdge',onlyEdge,'saveAnalysis',false,'matchWithFA',matchWithFA,...
        'minLifetime',minLifetime,'outputPath',['Colocalization' filesep outputPath]);
end
% re-express tracksNA so that each track has information for every frame
if ~isempty(iSDCProc)
    if ~isfield(tracksNA,'SDC_applied')
        disp('Applying stage drift correction ...')
        tracksNA = applyDriftToTracks(tracksNA, T); % need some other function....formatNATracks(tracksNAorg,detectedNAs,nFrames,T); 
    else
        disp('Stage drift correction was already applied to tracksNA.')
    end
end

%% tracion map stack
% Build the interpolated TFM matrix first and then go through each track
% First build overall TFM images
% We don't need to build traction map again if this one is already built
% during force field calculation process.
tMapFromOutput = false;
try
    tMapIn=forceFieldProc.loadChannelOutput('output','tMap');
    tmaxAuto = 0.8*max(tMapIn{1}(:));
    tmin = min(tMapIn{1}(:));
    cropInfo = [ceil(min(forceField(1).pos(:,1))),ceil(min(forceField(1).pos(:,2))),floor(max(forceField(1).pos(:,1))),floor(max(forceField(1).pos(:,2)))];
    tMapFromOutput = true;
catch
    [tMapIn, tmaxAuto, tmin, cropInfo] = generateHeatmapShifted(forceField,displField,0);
end
if isempty(tmax)
    tmax = tmaxAuto;
end
% Insert traction map in forceField.pos 
if ~isempty(iSDCProc)
    disp('Stage drift application to traction maps and cropping...')
    tracImage=(SDCProc.loadChannelOutput(iBeadChan,1)); %movieData.channels_(2).loadImage(ii);
else
    tracImage=MD.getChannel(iBeadChan).loadImage(1); 
end
[h,w] = size(tracImage);
tMap = zeros(h,w,nFrames);
for ii=1:nFrames
    cur_tMap = zeros(size(tracImage));
    % starts with original size of beads
    cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{ii}(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3));
    tMap(:,:,ii) = cur_tMap;
end
save([forcemapPath filesep 'tMap.mat'],'tMap','-v7.3');
%% Filter out tracks that is out of traction field
idxTracks = true(numel(tracksNA),1);
disp('Filtering with TFM boundary...')
for ii=1:numel(tracksNA)
    if any(round(tracksNA(ii).xCoord)<=cropInfo(1) | round(tracksNA(ii).xCoord)>=cropInfo(3) ...
            | round(tracksNA(ii).yCoord)<=cropInfo(2) | round(tracksNA(ii).yCoord)>=cropInfo(4))
        idxTracks(ii) = false;
    end
end
tracksNA=tracksNA(idxTracks);
%% Read force from tMap
% get the intensity
disp('Reading traction...')
tic
tracksNA = readIntensityFromTracks(tracksNA,tMap,2); % 2 means traction magnitude collection from traction stacks
toc
%% Read force uncertainty
forceParams=forceFieldProc.funParams_;
if strcmp(forceParams.method, 'FastBEM')
    try
        tracksNA = readForceUncertaintyFromTracks(pathForTheMovieDataFile,'tracksNA',tracksNA);
    catch
        getForceConfidence(pathForTheMovieDataFile)
        tracksNA = readForceUncertaintyFromTracks(pathForTheMovieDataFile,'tracksNA',tracksNA);
    end
end
%% Save SDC-shifted paxillin image stack
disp('Loading adhesion channel image stacks...')
tic
paxImgStack = zeros(h,w,nFrames);
for ii=1:nFrames
    if ~isempty(iSDCProc)
        paxImage=(SDCProc.loadChannelOutput(iChan,ii)); %movieData.channels_(2).loadImage(ii);
    else
        paxImage=MD.getChannel(iChan).loadImage(ii); 
    end
    paxImgStack(:,:,ii) = paxImage;
end
save([paxPath filesep 'paxImgStack.mat'],'paxImgStack','-v7.3');
toc
% disp('loading segmented FAs...')
FASegProc = FASegPackage.processes_{iFASeg};
%% showing force map in only traction region
if showAllTracks
    h1 = figure('color','w');
    set(h1, 'Position', [100 50 (cropInfo(3)-cropInfo(1)-2*band+1)*1.25 cropInfo(4)-cropInfo(2)-2*band+1])
    h2=figure;
    set(h2, 'Position', [100 50+round(1.4*cropInfo(4)-cropInfo(2)) (cropInfo(3)-cropInfo(1)-2*band+1) cropInfo(4)-cropInfo(2)-2*band+1])
end
disp('Matching with segmented adhesions...')
tic
for ii=1:nFrames
    tsMap = tMap(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band,ii);
    %%-------------Adhesion detection----------------------
    % loading paxillin image
    paxImage = paxImgStack(:,:,ii);
    paxImageCropped = paxImage(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band);

    bwPI4 = maskProc.loadChannelOutput(iChan,ii);
    
    % Get the mask for FAs
    if FASegProc.checkChannelOutput(1) && FASegProc.checkChannelOutput(2)
        iPaxChannel_adh=2;
    else 
        iPaxChannel_adh=iPaxChannel;
    end
    maskFAs = FASegProc.loadChannelOutput(iPaxChannel_adh,ii);
    % Apply stage drift correction to the cell segmentation
    % Get limits of transformation array
    if ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        % Apply subpixel-wise registration to original masks
        I = padarray(maskFAs, [maxY, maxX]);
        maskFAs = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
        Ibw = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(Ibw, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    end
    cropMask = bwPI4(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band);
%     cropMaskStack(:,:,ii) = cropMask;
    [B,~,nBD]  = bwboundaries(cropMask,'noholes');

    % mask for band from edge
    iMask = imcomplement(cropMask);
%     distFromEdge = bwdist(iMask);
%     bandwidth_pix = round(bandwidth*1000/MD.pixelSize_);
%     bandMask = distFromEdge <= bandwidth_pix;

%     maskAdhesion = maskFAs>0;
%     % erode once and separate some dumbel shaped adhesion into two
%     maskAdhesion2 = bwmorph(maskAdhesion,'hbreak',1);
%     % Somehow  this doesn't work. I need find out better way of showing
%     % adhesions.
% 
%     maskAdhesion = maskAdhesion2(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band);
%     maskFAs = maskFAs(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band);

%     if matchWithFA
%         % focal contact (FC) analysis
%         conn=4;
%         CC = bwconncomp(maskAdhesion,conn);
%         Adhs = regionprops(CC,'Area','Eccentricity','PixelIdxList','PixelList','Centroid' );
% 
%     %     %filter out focal adhesions that are in 1 um band along the edge: they might
%     %     %be connected nascent adhesions - SH 20140708
%     %     indAdh = true(numel(Adhs),1);
%     %     for k=1:numel(Adhs)
%     % %         plot(Adhs(k).Centroid(1),Adhs(k).Centroid(2),'mo')
%     %         if bandMaskFA(round(Adhs(k).Centroid(2)),round(Adhs(k).Centroid(1)))
%     %             indAdh(k,1) = false;
%     %             maskAdhesion(Adhs(k).PixelIdxList) = false;
%     %         end
%     %     end
%     %     Adhs = Adhs(indAdh);
% 
%     %     propFAs = regionprops(maskFAs,'Area','Eccentricity','PixelIdxList','PixelList' );
%         minFASize = round((1000/MD.pixelSize_)*(1000/MD.pixelSize_)); %adhesion limit=1 um2
% %         minFCSize = round((600/MD.pixelSize_)*(400/MD.pixelSize_)); %adhesion limit=.24 um2
%         minFCSize = round((300/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=.09 um2
% 
%         fcIdx = arrayfun(@(x) x.Area<minFASize & x.Area>minFCSize, Adhs);
%         FCs = Adhs(fcIdx);
%         FCForce = arrayfun(@(x) tsMap(x.PixelIdxList),FCs,'UniformOutput',false);
%         forceFC(ii) = mean(cell2mat(FCForce));
%         FCIdx = find(fcIdx);
%         adhBound = bwboundaries(maskAdhesion,conn,'noholes');    
% 
%         % for larger adhesions
%         faIdx = arrayfun(@(x) x.Area>=minFASize, Adhs);
%         FAs = Adhs(faIdx);
%         FAForce = arrayfun(@(x) tsMap(x.PixelIdxList),FAs,'UniformOutput',false);
%         forceFA(ii) = mean(cell2mat(FAForce));
%         FAIdx =  find(faIdx);
%         neighPix = 2;
% 
%         % Reading traction force at each track location
%         for k=1:numel(tracksNA)
%             if tracksNA(k).presence(ii)
%                 % decide if each track is associated with FC or FA
%                 p = 0;
%                 for jj=FCIdx'
%                     p=p+1;
%                     if any(round(tracksNA(k).xCoord(ii))-cropInfo(1)==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))-cropInfo(2)==Adhs(jj).PixelList(:,2))
%                         tracksNA(k).state{ii} = 'FC';
%     %                     tracksNA(k).forceMag(ii) = forceFC(p);
%                         tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
%                         tracksNA(k).FApixelList{ii} = [Adhs(jj).PixelList(:,1)+cropInfo(1) Adhs(jj).PixelList(:,2)+cropInfo(2)];
%                         tracksNA(k).adhBoundary{ii} = [adhBound{jj}(:,2)+cropInfo(1) adhBound{jj}(:,1)+cropInfo(2)];
%                         tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii))-cropInfo(2),round(tracksNA(k).xCoord(ii))-cropInfo(1));
%                     end
%                 end
%                 p = 0;
%                 for jj=FAIdx'
%                     p=p+1;
%                     if any(round(tracksNA(k).xCoord(ii))-cropInfo(1)==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))-cropInfo(2)==Adhs(jj).PixelList(:,2))
%                         tracksNA(k).state{ii} = 'FA';
%     %                     tracksNA(k).forceMag(ii) = forceFA(p);
%                         tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
%                         tracksNA(k).FApixelList{ii} = [Adhs(jj).PixelList(:,1)+cropInfo(1) Adhs(jj).PixelList(:,2)+cropInfo(2)];
%                         tracksNA(k).adhBoundary{ii} = [adhBound{jj}(:,2)+cropInfo(1) adhBound{jj}(:,1)+cropInfo(2)];
% %                         tracksNA(k).FApixelList{ii} = Adhs(jj).PixelList;
% %                         tracksNA(k).adhBoundary{ii} = adhBound{jj};
%                         tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii))-cropInfo(2),round(tracksNA(k).xCoord(ii))-cropInfo(1));
%                     end
%                 end
%             elseif ii>tracksNA(k).endingFrame && (strcmp(tracksNA(k).state{tracksNA(k).endingFrame},'FA')...
%                     || strcmp(tracksNA(k).state{tracksNA(k).endingFrame},'FC'))
%                 % starting from indexed maskFAs, find out segmentation that is
%                 % closest to the last track point.
%                 subMaskFAs = maskFAs==tracksNA(k).faID(tracksNA(k).endingFrame);
%                 if max(subMaskFAs(:))==0
%                     tracksNA(k).state{ii} = 'ANA';
%                     tracksNA(k).FApixelList{ii} = NaN;
%                     tracksNA(k).adhBoundary{ii} = NaN;
%                     continue
%                 else
%                     propSubMaskFAs = regionprops(subMaskFAs,'PixelList');
%                     minDist = zeros(length(propSubMaskFAs),1);
%                     for q=1:length(propSubMaskFAs)
%                         minDist(q) = min(sqrt((propSubMaskFAs(q).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame) - cropInfo(1))).^2 +...
%                             (propSubMaskFAs(q).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame) - cropInfo(2))).^2));
%                     end
%                     % find the closest segment
%                     [~,subMaskFAsIdx] = min(minDist);
%                     subAdhBound = bwboundaries(subMaskFAs,'noholes');    
%                     [~,closestPixelID] = min(sqrt((propSubMaskFAs(subMaskFAsIdx).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame) - cropInfo(1))).^2 +...
%                         (propSubMaskFAs(subMaskFAsIdx).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame) - cropInfo(2))).^2));
% 
%                     tracksNA(k).state{ii} = 'FC';
% %                     tracksNA(k).xCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,1);
% %                     tracksNA(k).yCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,2);
%                     tracksNA(k).FApixelList{ii} = propSubMaskFAs(subMaskFAsIdx).PixelList;
%                     tracksNA(k).adhBoundary{ii} = subAdhBound{subMaskFAsIdx};
%                 end
%             else
%                 tracksNA(k).forceMag(ii) = NaN;
%             end
%         end
%     end
    % force at the background
    % mask for band from edge
    distOutEdge = bwdist(cropMask);
    bandOut_pix = round(1*1000/MD.pixelSize_);
    bandOutMask = distOutEdge <= bandOut_pix;
    bandOutMask = bandOutMask & iMask;

    bandOutForce = tsMap((bandOutMask(:)));
    forceBGband(ii).mean = mean(bandOutForce);
    forceBGband(ii).err = std(bandOutForce);
    
    BGForce = tsMap((iMask(:)));
    forceBG(ii).mean = mean(BGForce);
    forceBG(ii).err = std(BGForce);
    
    %% recording features
    % get the point on the boundary closest to the adhesion
%     allBdPoints = [];
%     for kk=1:nBD
%         boundary = B{kk};
%         allBdPoints = [allBdPoints; boundary(:,2), boundary(:,1)];
%     end
%     for k=1:numel(tracksNA)
%         % distance to the cell edge
%         if tracksNA(k).presence(ii)
%             xCropped = tracksNA(k).xCoord(ii)-cropInfo(1)-band+1;
%             yCropped = tracksNA(k).yCoord(ii)-cropInfo(2)-band+1;
%             distToAdh = sqrt(sum((allBdPoints- ...
%                 ones(size(allBdPoints,1),1)*[xCropped, yCropped]).^2,2));
%             [minDistToBd,indMinBdPoint] = min(distToAdh);
%             tracksNA(k).distToEdge(ii) = minDistToBd;
%             tracksNA(k).closestBdPoint(ii,:) = allBdPoints(indMinBdPoint,:)+[cropInfo(1)+band-1 cropInfo(2)+band-1]; % this is lab frame of reference. (not relative to adhesion position)
%             % figure, subplot(2,1,1),plot(tracksNA((ii)).closestBdPoint(:,1),tracksNA((ii)).closestBdPoint(:,2)); hold on
% %                 plot(tracksNA((ii)).xCoord,tracksNA((ii)).yCoord,'g')
% %                 subplot(2,1,2),  fRange=tracksNA(ii).startingFrame:(tracksNA(ii).endingFrame); figure, plot(fRange,tracksNA(ii).distToEdge(fRange));
% %             threeAdjPoints = allBdPoints(indMinBdPoint-1:indMinBdPoint+1,:); % this is clockwise by default. I assume indMinBdPoint>1 and indMinBdPoint<end
% %             adjLine = polyfit(threeAdjPoints(:,1),threeAdjPoints(:,2),1);
% %             tracksNA(k).bdAdjVec(ii,:) = [1,adjLine]; % boundary adjacent vector, this might be errorneous if the local boundary is vertical
%         end
%     end
    
    % Showing for debug (TFM)
    if showAllTracks
%         h1 = figure('color','w');
%         set(h1, 'Position', [100 50 (cropInfo(3)-cropInfo(1)-2*band+1)*1.25 cropInfo(4)-cropInfo(2)-2*band+1])
        figure(h1), hold off
        ax1 = subplot('Position',[0 0 0.8 1]);
        hold off
        imshow(tsMap,[tmin tmax]), colormap jet;
        hold on;
    %     plot(pstruct.x,pstruct.y,'mo') % all peaks
        for kk=1:nBD
            boundary = B{kk};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        end
        % unit vector plot
        [reg_grid_coarse,~,~,~]=createRegGridFromDisplField(displField(ii),1); %2=2 times fine interpolation
        [grid_mat_coarse,iu_mat_coarse,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid_coarse);
        pos_coarse = [reshape(grid_mat_coarse(:,:,1),[],1) reshape(grid_mat_coarse(:,:,2),[],1)]; %dense
        disp_vec_coarse = [reshape(iu_mat_coarse(:,:,1),[],1) reshape(iu_mat_coarse(:,:,2),[],1)]; 
        [~,if_mat_coarse,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid_coarse);
        force_vec_coarse = [reshape(if_mat_coarse(:,:,1),[],1) reshape(if_mat_coarse(:,:,2),[],1)]; 

        [~,tmat_coarse, ~, ~] = interp_vec2grid(pos_coarse+disp_vec_coarse, force_vec_coarse,[],grid_mat_coarse); %1:cluster size

        tmat_coarse(end-floor(band/4)-1:end,:,:)=[];
        tmat_coarse(:,end-floor(band/4)-1:end,:)=[];
        tmat_coarse(1:1+floor(band/4)+1,:,:)=[];
        tmat_coarse(:,1:1+floor(band/4)+1,:)=[];
        grid_mat_coarse(end-floor(band/4)-1:end,:,:)=[];
        grid_mat_coarse(:,end-floor(band/4)-1:end,:)=[];
        grid_mat_coarse(1:1+floor(band/4)+1,:,:)=[];
        grid_mat_coarse(:,1:1+floor(band/4)+1,:)=[];

        tmat_vecx = reshape(tmat_coarse(:,:,1),[],1);
        tmat_vecy = reshape(tmat_coarse(:,:,2),[],1);
        pos_vecx = reshape(grid_mat_coarse(:,:,1),[],1);
        pos_vecy = reshape(grid_mat_coarse(:,:,2),[],1);
        forceScale=0.1*max(sqrt(tmat_vecx.^2+tmat_vecy.^2));
        quiver(pos_vecx-cropInfo(1),pos_vecy-cropInfo(2), tmat_vecx./forceScale,tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);

%         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        if matchWithFA
            idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
            idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{ii}),tracksNA(idAdhLogic));
            idAdh = find(idAdhLogic);
            idAdhCur = idAdh(idAdhCur);
            arrayfun(@(x) plot(x.adhBoundary{ii}(:,1)-cropInfo(1)-band+1,x.adhBoundary{ii}(:,2)-cropInfo(2)-band+1, 'k', 'LineWidth', 0.5),tracksNA(idAdhCur))

%             for k = FCIdx'
%         %         eachFA = maskFAs==k;
%         %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
%         %         for kk=1:nEachFA
%         %             adhBoundary = adhBound{kk};
%         %             plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
%         %         end
%                 adhBoundary = adhBound{k};
%                 plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
%             end
% 
%             for k = FAIdx'
%         %         eachFA = maskFAs==k;
%         %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
%         %         for kk=1:nEachFA
%         %             adhBoundary = adhBound{kk};
%         %             plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
%         %         end
%                 adhBoundary = adhBound{k};
%                 plot(adhBoundary(:,2), adhBoundary(:,1),  'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
%             end
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if onlyEdge
                    plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'r', 'LineWidth', 0.5)
                    plot(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                else
                    if strcmp(tracksNA(k).state{ii} , 'NA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'r', 'LineWidth', 0.5)
                        plot(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                        if showTrackID
                            text(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1+4,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,['\leftarrow' num2str(k)],'Color','r','FontSize',6.5)
                        end
                    elseif strcmp(tracksNA(k).state{ii} , 'FC')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
    %                     plot(tracksNA(k).xCoord(ii)-cropInfo(1),tracksNA(k).yCoord(ii)-cropInfo(2),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{ii} , 'FA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) 
    %                     plot(tracksNA(k).xCoord(ii)-cropInfo(1),tracksNA(k).yCoord(ii)-cropInfo(2),'ko','MarkerSize',markerSize, 'LineWidth', 0.5)
                    end
                end
            end
        end
        % Scale bar 2000nm
        subplot('Position',[0.8 0.1 0.1 0.8])
        axis tight
        caxis([tmin tmax]), axis off
        hc = colorbar('West');
        set(hc,'Fontsize',12)
        hold on;

%         syFigureStyle(h1,ax1,cropInfo(4)-cropInfo(2)-2*band+1);
        print('-dtiff', '-r150', strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat),'.tif'));
    %     hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/forcePeakFig',num2str(ii,iiformat)),'-v7.3')
        print('-depsc2', '-r150', strcat(epsPath,'/forcePeak',num2str(ii,iiformat),'.eps'));
        hold off
    end    
    if showAllTracks % for pax image
%         h2=figure;
%         set(h2, 'Position', [100 50+round(1.4*cropInfo(4)-cropInfo(2)) (cropInfo(3)-cropInfo(1)-2*band+1) cropInfo(4)-cropInfo(2)-2*band+1])
        figure(h2), hold off
        %Scale bar 2 um
    %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
        darkeningFactor=0.5;
        paxImageCroppedInverted = imcomplement(paxImageCropped);
        minPax = min(paxImageCroppedInverted(:));
        maxPax = max(paxImageCroppedInverted(:));
        minPax = double(minPax)+double(darkeningFactor*(maxPax-minPax));
        imshow(paxImageCroppedInverted,[minPax maxPax]), hold on
        line([10 10+round(2000/MD.pixelSize_)],[15 15],'LineWidth',2,'Color',[0,0,0])
        
        for kk=1:nBD
            boundary = B{kk};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        end
    %     plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
        if matchWithFA
            idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
            idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{ii}),tracksNA(idAdhLogic));
            idAdh = find(idAdhLogic);
            idAdhCur = idAdh(idAdhCur);
            arrayfun(@(x) plot(x.adhBoundary{ii}(:,1)-cropInfo(1)-band+1,x.adhBoundary{ii}(:,2)-cropInfo(2)-band+1, 'k', 'LineWidth', 0.5),tracksNA(idAdhCur))
%             for k = FCIdx'
%                 adhBoundary = adhBound{k};
%                 plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
%             end
%             % for larger adhesions
%             for k = FAIdx'
%                 adhBoundary = adhBound{k};
%                 plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
%             end
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if onlyEdge
                    plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'r', 'LineWidth', 0.5)
                    plot(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                else
                    if strcmp(tracksNA(k).state{ii} , 'NA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'r', 'LineWidth', 0.5)
                        plot(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                        if showTrackID
                            text(tracksNA(k).xCoord(ii)-cropInfo(1)-band+1+4,tracksNA(k).yCoord(ii)-cropInfo(2)-band+1,['\leftarrow' num2str(k)],'Color','r','FontSize',6.5)
                        end
                    elseif strcmp(tracksNA(k).state{ii} , 'FC')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
    %                     plot(tracksNA(k).xCoord(ii)-cropInfo(1),tracksNA(k).yCoord(ii)-cropInfo(2),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{ii} , 'FA')
                        % drawing tracks
                        plot(tracksNA(k).xCoord(1:ii)-cropInfo(1)-band+1,tracksNA(k).yCoord(1:ii)-cropInfo(2)-band+1, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) 
    %                     plot(tracksNA(k).xCoord(ii)-cropInfo(1),tracksNA(k).yCoord(ii)-cropInfo(2),'ko','MarkerSize',markerSize, 'LineWidth', 0.5)
                    end
                end
            end
        end
%         syFigureStyle(h2,gca,cropInfo(4)-cropInfo(2)-2*band+1);

        print('-depsc2', '-r150', strcat(epsPath,'/pax',num2str(ii,iiformat),'.eps'));
        print('-dtiff', '-r150', strcat(paxtifPath,'/pax',num2str(ii,iiformat),'.tif'));
    %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
        hold off
%         close(h1)
%         clear h1
%         close(h2)
%         clear h2
    end
%     imwrite(uint16(round(tsMap*2^15/3500)),strcat(forcemapPath,'/force',num2str(ii,iiformat),'max',num2str(tmax),'.tif'));
    imwrite(uint16(paxImageCropped),strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
    if ii==1
        if tMapFromOutput
            save(strcat(dataPath,'/cropInfo.mat'),'cropInfo');
        else
            cropInfo = [cropInfo(1) cropInfo(2) cropInfo(3)-cropInfo(1) cropInfo(4)-cropInfo(2)];
            save(strcat(dataPath,'/cropInfo.mat'),'cropInfo');
        end
    end
end
toc
%% protrusion/retraction information - most of these are now done in analyzeAdhesionsMaturation
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.
% for ii=1:numel(tracksNA)
%     idxZeros = tracksNA(ii).closestBdPoint(:)==0;
%     tracksNA(ii).closestBdPoint(idxZeros)=NaN;
% end

disp('Post-analysis on adhesion movement and cross-correlation between fluorescence intensity and traction...')
deltaT = MD.timeInterval_; % sampling rate (in seconds, every deltaT seconds)
tIntervalMin = deltaT/60; % in min
periodMin = 1;
periodFrames = floor(periodMin/tIntervalMin); % early period in frames
timeInterval = deltaT/60; % in min
for k=1:numel(tracksNA)
    % cross-correlation scores
    presIdx = logical(tracksNA(k).presence);
%         [curCC,curLag] = xcorr(tracksNA(k).ampTotal(presIdx),tracksNA(k).forceMag(presIdx),'biased');
    maxLag = ceil(tracksNA(k).lifeTime/2);
    [curCC,curBounds,curLag,curP]  = nanCrossCorrelation(tracksNA(k).ampTotal(presIdx),tracksNA(k).forceMag(presIdx),'corrType','Pearson','maxLag',maxLag);
    tracksNA(k).CCscore = curCC;
    tracksNA(k).CCbounds = curBounds;
    tracksNA(k).CClag = curLag;
    tracksNA(k).CC_p = curP;
    [tracksNA(k).CCscoreMax,curMaxInd] = max(curCC);
    tracksNA(k).CCmaxLag = curLag(curMaxInd);
    % kurtosis around the peak to measure 'peakness' of the peak
    minIndKur = max(1, curMaxInd-10);
    maxIndKur = min(length(curCC), curMaxInd+10);
    tracksNA(k).CCkurtosis = kurtosis(curCC(minIndKur:maxIndKur));
    % the sharper the peak is, the higher CCkurtosis is.
    
%     [curAC,curBoundsAC,curLagAC,curAC_P]  = nanCrossCorrelation(tracksNA(k).ampTotal(presIdx),tracksNA(k).ampTotal(presIdx),'corrType','Pearson','maxLag',maxLag);
%     [curFAC,curBoundsFAC,curLagFAC,curFAC_P]  = nanCrossCorrelation(tracksNA(k).forceMag(presIdx),tracksNA(k).forceMag(presIdx),'corrType','Pearson','maxLag',maxLag);
    try
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
        tracksNA(k).lifeTime = eF-sF+1;    
    end
%     % get the instantaneous velocity
%     % get the distance first
%     curTrackVelLength=sum(presIdx)-1;
%     distTrajec=zeros(curTrackVelLength,1);
%     presIdxSeq = find(presIdx);
%     for kk=1:curTrackVelLength
%         real_kk = presIdxSeq(kk);
%         distTrajec(kk) = sqrt(sum((tracksNA(k).closestBdPoint(real_kk+1,:)- ...
%             tracksNA(k).closestBdPoint(real_kk,:)).^2,2));
%         lastPointIntX = round(tracksNA(k).closestBdPoint(real_kk+1,1));
%         lastPointIntY = round(tracksNA(k).closestBdPoint(real_kk+1,2));
%         if cropMaskStack(lastPointIntY-(floor(cropInfo(2))+band-1),lastPointIntX-(floor(cropInfo(1))+band-1),real_kk) %if the last point is in the first mask, it is inward
%             distTrajec(kk) = -distTrajec(kk);
%         end
%     end
%     if any(distTrajec~=0)
%         [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
%         if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
%             tracksNA(k).isProtrusion = true;
%         else
%             tracksNA(k).isProtrusion = false;
%         end
%         % average velocity (positive for protrusion)
%         curProtVel = (Protrusion.Veloc); curProtVel(isnan(curProtVel))=0;
%         curProtPersTime = (Protrusion.persTime); curProtPersTime(isnan(curProtPersTime))=0;
%         curRetVel = (Retraction.Veloc); curRetVel(isnan(curRetVel))=0;
%         curRetPersTime = (Retraction.persTime); curRetPersTime(isnan(curRetPersTime))=0;
% 
%         tracksNA(k).edgeVel = (mean(curProtVel.*curProtPersTime)-mean(curRetVel.*curRetPersTime))/mean([curProtPersTime;curRetPersTime]);
%     else
%         tracksNA(k).edgeVel = 0;
%     end
%     % lifetime information
%     tracksNA(k).lifeTime = tracksNA(k).endingFrame-tracksNA(k).startingFrame+1;    
%     % Inital intensity slope for one min
%     timeInterval = deltaT/60; % in min
    earlyPeriod = floor(1/timeInterval); % frames per minute
    lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
    lastFrameFromOne = lastFrame - sF+1;
%     [curR,curM] = regression(timeInterval*(1:lastFrameFromOne),tracksNA(k).amp(tracksNA(k).startingFrame:lastFrame));
%     tracksNA(k).ampSlope = curM; % in a.u./min
%     tracksNA(k).ampSlopeR = curR; % Pearson's correlation coefficient
    [curForceR,curForceM] = regression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame));
    tracksNA(k).forceSlope = curForceM; % in Pa/min
    tracksNA(k).forceSlopeR = curForceR; % Pearson's correlation coefficient
    
%     curEndFrame = min(tracksNA(k).startingFrameExtra+periodFrames-1,tracksNA(k).endingFrame);
%     curEarlyPeriod = curEndFrame - tracksNA(k).startingFrameExtra+1;
%     [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:curEndFrame));
%     tracksNA(k).earlyAmpSlope = curM; % in a.u./min
% 
%     curStartFrame = max(tracksNA(k).startingFrame,tracksNA(k).endingFrameExtra-periodFrames+1);
%     curLatePeriod = tracksNA(k).endingFrameExtra - curStartFrame+1;
%     [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksNA(k).ampTotal(curStartFrame:tracksNA(k).endingFrameExtra));
%     tracksNA(k).lateAmpSlope = curMlate; % in a.u./min
% 
%     curEndFrame = min(tracksNA(k).startingFrame+periodFrames-1,tracksNA(k).endingFrame);
%     curEarlyPeriod = curEndFrame - tracksNA(k).startingFrame+1;
%     [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).distToEdge(tracksNA(k).startingFrame:curEndFrame));
%     tracksNA(k).distToEdgeSlope = curMdist; % in a.u./min
%     tracksNA(k).distToEdgeChange = (tracksNA(k).distToEdge(end)-tracksNA(k).distToEdge(tracksNA(k).startingFrame)); % in pixel
%     sF=tracksNA(k).startingFrame;
%     eF=tracksNA(k).endingFrame;
%     % Determining protrusion/retraction based on closestBdPoint and [xCoord
%     % yCoord]. If the edge at one time point does not cross the adhesion
%     % tracks over entire life time, there is no problem. But that's not
%     % always the case: adhesion tracks crosses the cell boundary at first
%     % or last time point. And we don't know if the adhesion tracks are
%     % moving in the direction of protrusion or retraction. Thus, what I
%     % will do is to calculate the inner product of vectors from the
%     % adhesion tracks from the closestBdPoint at the first frame. If the
%     % product is positive, it means both adhesion points are in the same
%     % side. And if distance from the boundary to the last point is larger
%     % than the distance to the first point, it means the adhesion track is
%     % retracting. In the opposite case, the track is protruding (but it
%     % will happen less likely because usually the track would cross the
%     % first frame boundary if it is in the protrusion direction). 
%     
%     % We need to find out boundary points projected on the line of adhesion track 
% %     edge = [tracksNA(k).xCoord(sF) tracksNA(k).yCoord(sF) tracksNA(k).xCoord(eF) tracksNA(k).yCoord(eF)];
% %     [aLine,bLine,cLine,fitnessLine]=ortho2Dlinefit(tracksNA(k).xCoord,tracksNA(k).yCoord);
%     fitobj = fit(tracksNA(k).xCoord(sF:eF)',tracksNA(k).yCoord(sF:eF)','poly1');
%     x0=nanmedian(tracksNA(k).xCoord);
%     y0=fitobj(x0);
%     dx = 1;
%     dy = fitobj.p1;
%     trackLine = createLineGeom2d(x0,y0,dx,dy);
%     
% %     trackLine = edgeToLine(edge);
%     firstBdPoint = [tracksNA(k).closestBdPoint(sF,1) tracksNA(k).closestBdPoint(sF,2)];
%     firstBdPointProjected = projPointOnLine(firstBdPoint, trackLine);
%     lastBdPoint = [tracksNA(k).closestBdPoint(eF,1) tracksNA(k).closestBdPoint(eF,2)];
%     lastBdPointProjected = projPointOnLine(lastBdPoint, trackLine);
%     
%     firstAdhToFirstBdPoint = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)];
%     lastAdhToFirstBdPoint = [tracksNA(k).xCoord(eF)-firstBdPointProjected(1), tracksNA(k).yCoord(eF)-firstBdPointProjected(2)];
%     firstAdhToLastBdPoint = [tracksNA(k).xCoord(sF)-lastBdPointProjected(1), tracksNA(k).yCoord(sF)-lastBdPointProjected(2)];
%     lastAdhToLastBdPoint = [tracksNA(k).xCoord(eF)-lastBdPointProjected(1), tracksNA(k).yCoord(eF)-lastBdPointProjected(2)];
%     firstBDproduct=firstAdhToFirstBdPoint*lastAdhToFirstBdPoint';
%     lastBDproduct=firstAdhToLastBdPoint*lastAdhToLastBdPoint';
%     if firstBDproduct>0 && firstBDproduct>lastBDproduct% both adhesion points are in the same side
%         tracksNA(k).advanceDist = (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5 - ...
%                                                         (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%         tracksNA(k).edgeAdvanceDist = (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                         (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%     else
%         if lastBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
%             tracksNA(k).advanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                             (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5; % in pixel
%             tracksNA(k).edgeAdvanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                             (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%             % this code is nice to check:
% %             figure, imshow(paxImgStack(:,:,tracksNA(k).endingFrame),[]), hold on, plot(tracksNA(k).xCoord,tracksNA(k).yCoord,'w'), plot(tracksNA(k).closestBdPoint(:,1),tracksNA(k).closestBdPoint(:,2),'r')
% %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
% %             plot(tracksNA(k).xCoord(tracksNA(k).startingFrame),tracksNA(k).yCoord(tracksNA(k).startingFrame),'yo'),plot(tracksNA(k).xCoord(tracksNA(k).endingFrame),tracksNA(k).yCoord(tracksNA(k).endingFrame),'mo')
%         else
%             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
%             firstAdhToFirstBdPoint = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(sF,2)];
%             lastAdhToFirstBdPoint = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(sF,2)];
%             firstAdhToLastBdPoint = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(eF,2)];
%             lastAdhToLastBdPoint = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(eF,2)];
%             firstBDproduct=firstAdhToFirstBdPoint*lastAdhToFirstBdPoint';
%             lastBDproduct=firstAdhToLastBdPoint*lastAdhToLastBdPoint';
%             if firstBDproduct>lastBDproduct
%                 tracksNA(k).advanceDist = (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5 - ...
%                                                                 (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%                 tracksNA(k).edgeAdvanceDist = (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                                 (lastAdhToFirstBdPoint(1)^2 + lastAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%             else
%                 tracksNA(k).advanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                                 (lastAdhToLastBdPoint(1)^2 + lastAdhToLastBdPoint(2)^2)^0.5; % in pixel
%                 tracksNA(k).edgeAdvanceDist = (firstAdhToLastBdPoint(1)^2 + firstAdhToLastBdPoint(2)^2)^0.5 - ...
%                                                                 (firstAdhToFirstBdPoint(1)^2 + firstAdhToFirstBdPoint(2)^2)^0.5; % in pixel
%             end
%         end
%     end
    
end
disp('Saving...')
save(outputFile,'tracksNA','-v7.3');
disp('Done!')
if plotEachTrack
    r1 = 50;
    h2=figure;
    for k=1:round(numel(tracksNA)/3)%[7 15 32 127 129]%afterFAK[1 2 18 30 36 11 12 14 39 41 44]% 
        % try to crop window around the track
        if isempty(tracksNA(k).startingFrame) 
            continue
        end
        if strcmp(tracksNA(k).state{tracksNA(k).startingFrame},'FC') || strcmp(tracksNA(k).state{tracksNA(k).startingFrame},'FA')
            continue
        end
        fstart = max(tracksNA(k).startingFrame-20,1);
        fend = min(tracksNA(k).endingFrame,nFrames);
        iSF = tracksNA(k).startingFrame;
        wRoi = min(tracksNA(k).xCoord(iSF)-cropInfo(1),r1)...
            +min(cropInfo(3)-cropInfo(1)+1-(tracksNA(k).xCoord(iSF)-cropInfo(1)),r1);
        hRoi = min(tracksNA(k).yCoord(iSF)-cropInfo(2),r1)...
            +min(cropInfo(4)-cropInfo(2)+1-(tracksNA(k).yCoord(iSF)-cropInfo(2)),r1);
        set(h2,'Units','inches')
        set(h2,'PaperPositionMode','auto')
        set(h2, 'Position', [1,1,wRoi/(hRoi*2), 1])

        eachPaxPath = [paxPath filesep '/track' num2str(k,iiformat)];
        eachEpsPath = [epsPath filesep '/track' num2str(k,iiformat)];
        if ~exist(eachPaxPath,'dir') || ~exist(eachEpsPath,'dir')
            mkdir(eachPaxPath);
            mkdir(eachEpsPath);
        end
        for j=fstart:fend
            % show each track with paxillin images
            tsMap = imread(strcat(forcemapPath,'/force',num2str(j,iiformat),'max',num2str(tmax),'.tif'));
            tsMap = double(tsMap)*3500/(2^15); %converting to Pa
            paxImageCropped = imread(strcat(paxPath,'/pax',num2str(j,iiformat),'.tif'));
            xminROI = round(max(1,tracksNA(k).xCoord(iSF)-cropInfo(1)-(r1-1))); 
            xmaxROI = round(min(cropInfo(3)-cropInfo(1)+1,tracksNA(k).xCoord(iSF)-cropInfo(1)+r1)); 
            yminROI = round(max(1,tracksNA(k).yCoord(iSF)-cropInfo(2)-(r1-1)));
            ymaxROI = round(min(cropInfo(4)-cropInfo(2)+1,tracksNA(k).yCoord(iSF)-cropInfo(2)+r1));
            paxImageCropped2 = paxImageCropped(yminROI:ymaxROI,xminROI:xmaxROI);
            tsMapCropped = tsMap(yminROI:ymaxROI,xminROI:xmaxROI);
            
            ha1 = subplot('position',[0  0.5  1  0.5]);
            if j==fstart
                invPaxImageCropped2 = imcomplement(paxImageCropped2);
                lastPax = imread(strcat(paxPath,'/pax',num2str(nFrames,iiformat),'.tif'));
                lastPaxCropped = lastPax(yminROI:ymaxROI,xminROI:xmaxROI);
                invLastPaxCropped = imcomplement(lastPaxCropped);
                pmax = max([invPaxImageCropped2(:); invLastPaxCropped(:)]);
                pmin = min([invPaxImageCropped2(:); invLastPaxCropped(:)]);
                if j~= iSF && strcmp(tracksNA(k).state{iSF} , 'NA') % remember the first NA's position
                    xFirst = tracksNA(k).xCoord(iSF);
                    yFirst = tracksNA(k).yCoord(iSF);
%                     bkgAmpFirst = tracksNA(k).bkgAmp(iSF);
                    ampFirst = tracksNA(k).amp(iSF);
                    sigmaFirst = tracksNA(k).sigma(iSF);
                end
            end
            imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1),colormap(ha1,'gray');freezeColors; hold(ha1,'on')
            if strcmp(tracksNA(k).state{j} , 'BA')
                % drawing tracks
                plot(ha1,xFirst-cropInfo(1)-xminROI,yFirst-cropInfo(2)-yminROI,'g', 'LineWidth', 0.5)
                % remembering adhesion intensity and TF
                ynmin = max(1,round(yFirst)-cropInfo(2)-neighPix);
                ynmax = min(size(tsMap,1),round(yFirst)-cropInfo(2)+neighPix);
                xnmin = max(1,round(xFirst)-cropInfo(1)-neighPix);
                xnmax = min(size(tsMap,2),round(xFirst)-cropInfo(1)+neighPix);
                forceNeigh = tsMap(ynmin:ynmax,xnmin:xnmax);
                tracksNA(k).forceMag(j) = max(forceNeigh(:));    
                %here, now I try to do gaussian fit that
                %pointsourceDetection used
                paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                pstruct = fitGaussians2D(double(paxImageCropped), xFirst-cropInfo(1), yFirst-cropInfo(2), ampFirst*0.1, sigmaFirst*0.5, min(paxNeigh(:)),'xyAc');
                tracksNA(k).amp(j) = pstruct.A(1);
                if isnan(pstruct.A)
                    tracksNA(k).amp(j) = double(paxImageCropped(round(yFirst)-cropInfo(2),round(xFirst)-cropInfo(1)))-min(paxNeigh(:));
                end
            elseif strcmp(tracksNA(k).state{j} , 'NA')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'r', 'LineWidth', 0.5)
            elseif strcmp(tracksNA(k).state{j} , 'FC')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
            elseif strcmp(tracksNA(k).state{j} , 'FA')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'k', 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'k', 'LineWidth', 0.5) %adhesion boundary
            end
            ha2 = subplot('position',[0  0  1 0.5]);
            if j==fstart
                tmax2 = max(min(2500,max(tsMapCropped(:))*0.9),tmaxEach);
                tmin2 = min(min(tsMapCropped(:)),100);
            end
            imshow(tsMapCropped,[tmin2 tmax2],'Parent', ha2), colormap(ha2,'jet');freezeColors; hold(ha2,'on')
%             set(h2,'CurrentAxes',ha2)
            
            if strcmp(tracksNA(k).state{j}, 'NA')
                % drawing tracks
                plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'r', 'LineWidth', 0.5)
            elseif strcmp(tracksNA(k).state{j}, 'FC')
                if j==iSF
                    continue
                end
                % drawing tracks
                plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
            elseif strcmp(tracksNA(k).state{j}, 'FA')
                if j==iSF
                    continue
                end
                % drawing tracks
                plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'k', 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'k', 'LineWidth', 0.5) %adhesion boundary
            end
            if j==fend
                cb = colorbar('location','East');
                cbfreeze(cb)
                line([3 3+round(2000/MD.pixelSize_)],[hRoi-3 hRoi-3],'LineWidth',2,'Color',[0,0,0])
            end
            curRenderer = get(h2,'Renderer');
            if ~strcmp(curRenderer,'painters')
                set(h2,'Renderer','painters')
            end
            print('-depsc2', '-r150', strcat(eachEpsPath,'/trackFrame',num2str(j,iiformat),'.eps'));
            print('-dtiff', '-r150', strcat(eachPaxPath,'/trackFrame',num2str(j,iiformat),'.tif'));
            hold(ha1,'off')
            hold(ha2,'off')
        end
        seeMore = input('Do you want to look at video until the end of the movie? (y/(n))','s');
        if isempty(seeMore)
            seeMore = 'n';
        end
        if strcmp(seeMore,'y')
            if fend < nFrames
                for j=fend+1:nFrames
                    % show each track with paxillin images
                    tsMap = imread(strcat(forcemapPath,'/force',num2str(j,iiformat),'max',num2str(tmax),'.tif'));
                    tsMap = double(tsMap)*3500/(2^15); %converting to Pa
                    
                    %acquiring force
                    ynmin = max(1,round(tracksNA(k).yCoord(fend))-cropInfo(2)-neighPix);
                    ynmax = min(size(tsMap,1),round(tracksNA(k).yCoord(fend))-cropInfo(2)+neighPix);
                    xnmin = max(1,round(tracksNA(k).xCoord(fend))-cropInfo(1)-neighPix);
                    xnmax = min(size(tsMap,2),round(tracksNA(k).xCoord(fend))-cropInfo(1)+neighPix);
                    forceNeigh = tsMap(ynmin:ynmax,xnmin:xnmax);
                    tracksNA(k).forceMag(j) = max(forceNeigh(:));    
                   
                    paxImageCropped = imread(strcat(paxPath,'/pax',num2str(j,iiformat),'.tif'));
                    xminROI = round(max(1,tracksNA(k).xCoord(iSF)-cropInfo(1)-(r1-1))); 
                    xmaxROI = round(min(cropInfo(3)-cropInfo(1)+1,tracksNA(k).xCoord(iSF)-cropInfo(1)+r1)); 
                    yminROI = round(max(1,tracksNA(k).yCoord(iSF)-cropInfo(2)-(r1-1)));
                    ymaxROI = round(min(cropInfo(4)-cropInfo(2)+1,tracksNA(k).yCoord(iSF)-cropInfo(2)+r1));

                    paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                    xLast = tracksNA(k).xCoord(j-1);
                    yLast = tracksNA(k).yCoord(j-1);
                    ampLast = tracksNA(k).amp(j-1);
                    sigmaLast = 2.1; %tracksNA(k).sigma(j-1);
                    pstruct = fitGaussians2D(double(paxImageCropped), xLast-cropInfo(1), yLast-cropInfo(2), ampLast, sigmaLast*0.5, min(paxNeigh(:)),'xyAc');
                    tracksNA(k).amp(j) = pstruct.A(1);
                    tracksNA(k).xCoord(j) = pstruct.x(1)+cropInfo(1);
                    tracksNA(k).yCoord(j) = pstruct.y(1)+cropInfo(2);
                    if isnan(pstruct.A)
                        tracksNA(k).amp(j) = double(paxImageCropped(round(yLast)-cropInfo(2),round(xLast)-cropInfo(1)))-min(paxNeigh(:));
                        tracksNA(k).xCoord(j) = tracksNA(k).xCoord(j-1);
                        tracksNA(k).yCoord(j) = tracksNA(k).yCoord(j-1);
                    end
                    
                    paxImageCropped2 = paxImageCropped(yminROI:ymaxROI,xminROI:xmaxROI);
                    tsMapCropped = tsMap(yminROI:ymaxROI,xminROI:xmaxROI);
                    ha1 = subplot('position',[0  0.5  1  0.5]);
                    imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1),colormap(ha1,'gray');freezeColors; hold(ha1,'on')
                    if strcmp(tracksNA(k).state{fend} , 'NA')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'r', 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{fend} , 'FC')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'k', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    end
                    ha2 = subplot('position',[0  0  1 0.5]);
                    imshow(tsMapCropped,[tmin2 tmax2],'Parent', ha2), colormap(ha2,'jet');freezeColors; hold(ha2,'on')
                    if strcmp(tracksNA(k).state{fend} , 'NA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'r', 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{fend} , 'FC')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-cropInfo(1)-xminROI,tracksNA(k).yCoord(1:j)-cropInfo(2)-yminROI,'k', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    end
                    curRenderer = get(h2,'Renderer');
                    if ~strcmp(curRenderer,'painters')
                        set(h2,'Renderer','painters')
                    end
                    if j==nFrames
                        cb = colorbar('location','East');
                        cbfreeze(cb)
                        line([3 3+round(2000/MD.pixelSize_)],[hRoi-3 hRoi-3],'LineWidth',2,'Color',[0,0,0])
                    end
                    print('-depsc2','-r150', strcat(eachEpsPath,'/trackFrame',num2str(j,iiformat),'.eps'));
                    print('-dtiff', '-r150', strcat(eachPaxPath,'/trackFrame',num2str(j,iiformat),'.tif'));
                    hold(ha1,'off')
                    hold(ha2,'off')
               end
            else
                disp('It is already the end of movie!')
            end
        end
        % inquire if the estimated state is right
        disp(['The state of this adhesion is : track' num2str(k)])
        disp([num2cell(tracksNA(k).iFrame(fstart:fend)') tracksNA(k).state(fstart:fend)'])
        strFA = input('Is this state well describing what you see in the movie (no over-estimated FAs or noisy NAs)? ((y)/n)','s');
        while strcmp(strFA,'n')
            iFrames = input('Which frames do you want to change? n1:n2  ');
            state = input('What is the state in those range? (e.g. BA, NA, FC, FA)  ','s');
            for jj=iFrames
                tracksNA(k).state{jj} = state;
            end
            disp(['Now, the state of this adhesion is :' num2str(k)])
            disp([num2cell(tracksNA(k).iFrame(fstart:fend)') tracksNA(k).state(fstart:fend)'])
            strFA = input('Is this state well describing what you see in the movie (no over-estimated FAs or noisy NAs)? ((y)/n)','s');
        end
        hold off
    end
end
end
%

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