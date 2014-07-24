function [tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,showAllTracks,plotEachTrack,tmaxEach)
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

if nargin < 3
    band = 16;
    tmax=[];
    plotEachTrack = false;
    showAllTracks = true;
elseif nargin <4
    tmax=[];
    plotEachTrack = false;
    showAllTracks = true;
elseif nargin<5
    plotEachTrack = false;
    showAllTracks = true;
elseif nargin<6
    plotEachTrack = false;
end
% if isempty(pstruct_NAcell)
%     bFirst=true;
% else
%     bFirst=false;
% end
% bandwidth = 20; %um
%% Data Set up
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
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
iUTrack = MD.getPackageIndex('UTrackPackage');
% If not run, run the package
iFASeg = 6;
if isempty(iUTrack)
    MD.addPackage(UTrackPackage(MD));
    iUTrack =  MD.getPackageIndex('UTrackPackage');
    MD.getPackage(iUTrack).createDefaultProcess(1)
    params = MD.getPackage(iUTrack).getProcess(1).funParams_;
    params.ChannelIndex = 2; %paxillin
    params.RedundancyRadius(2) = 0.1; % in px
    params.filterSigma(2) = 1.1; % in nm
    MD.getPackage(iUTrack).getProcess(1).setPara(params);
    MD.getPackage(iUTrack).getProcess(1).run();
    MD.save;
    
    MD.getPackage(iUTrack).createDefaultProcess(2)
    params = MD.getPackage(iUTrack).getProcess(2).funParams_;
    params.ChannelIndex = 2; %paxillin
    params.probDim = 0.1; % in px
    params.gapCloseParam.timeWindow = 7;
    params.gapCloseParam.minTrackLen = 4;
    params.gapCloseParam.diagnostics = 0;
    params.gapCloseParam.mergeSplit = 0;
    params.costMatrices(1).parameters.minSearchRadius = 2;
    params.costMatrices(2).parameters.minSearchRadius = 2;
    
    MD.getPackage(iUTrack).getProcess(2).setPara(params);
    MD.getPackage(iUTrack).getProcess(2).run();
    MD.save;
    trackNAProc = MD.getProcess(MD.getProcessIndex('TrackingProcess'));
elseif isempty(MD.getPackage(iUTrack).getProcess(2))
    MD.getPackage(iUTrack).createDefaultProcess(2)
    params = MD.getPackage(iUTrack).getProcess(2).funParams_;
    params.ChannelIndex = 2; %paxillin
    params.probDim = 2; % in px
    params.gapCloseParam.timeWindow = 5;
    params.gapCloseParam.minTrackLen = 4;
    params.gapCloseParam.diagnostics = 0;
    params.gapCloseParam.mergeSplit = 0;
    params.costMatrices(1).parameters.minSearchRadius = 2;
    params.costMatrices(2).parameters.minSearchRadius = 2;
    
    MD.getPackage(iUTrack).getProcess(2).setPara(params);
    MD.getPackage(iUTrack).getProcess(2).run();
    MD.save;
    trackNAProc = MD.getProcess(MD.getProcessIndex('TrackingProcess'));
else
    trackNAProc = MD.getProcess(MD.getProcessIndex('TrackingProcess'));
end

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

%Find the maximum force.
tmin = 100000;
[reg_grid1,~,~,~]=createRegGridFromDisplField(displField(1),1); %2=2 times fine interpolation

% band width for cutting border
%     band=4;
if isempty(tmax)
    tmax = 0;
    for ii = 1:nFrames
       %Load the saved body force map.
        [~,fmat, ~, ~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1); %1:cluster size
        fnorm = (fmat(:,:,1).^2 + fmat(:,:,2).^2).^0.5;
        % Boundary cutting - I'll take care of this boundary effect later
        fnorm(end-round(band/2):end,:)=[];
        fnorm(:,end-round(band/2):end)=[];
        fnorm(1:1+round(band/2),:)=[];
        fnorm(:,1:1+round(band/2))=[];
        fnorm_vec = reshape(fnorm,[],1); 

        tmax = max(tmax,max(fnorm_vec));
        tmin = min(tmin,min(fnorm_vec));
    end
    display(['Force maximum = ' num2str(tmax)])
else
    for ii = 1:nFrames
       %Load the saved body force map.
        [~,fmat, ~, ~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1); %1:cluster size
        fnorm = (fmat(:,:,1).^2 + fmat(:,:,2).^2).^0.5;
        % Boundary cutting - I'll take care of this boundary effect later
        fnorm(end-round(band/2):end,:)=[];
        fnorm(:,end-round(band/2):end)=[];
        fnorm(1:1+round(band/2),:)=[];
        fnorm(:,1:1+round(band/2))=[];
        fnorm_vec = reshape(fnorm,[],1); 

        tmin = min(tmin,min(fnorm_vec));
    end
end

[reg_grid,~,~,~]=createRegGridFromDisplField(displField(ii),4); %4=4 times fine interpolation

%     h2 = figure;
hold off
%     set(h2, 'Position', [100+imSizeX*10/9 100 imSizeX imSizeY])
iiformat = ['%.' '3' 'd'];
%     paxLevel = zeros(nFrames,1);
minLifetime = min(nFrames,3);
markerSize = 3;
% tracks
iPaxChannel = 2; % this should be intentionally done in the analysis level
if ~trackNAProc.checkChannelOutput(iPaxChannel)
    iPaxChannel = 1;
end
% filter out tracks that have lifetime less than 2 frames

% Build the interpolated TFM matrix first and then go through each track
% ...

% Find out force at each tracksNA at this frame
% find first the index of relevant tracks from seqOfEvents
% tracksNA = trackNAProc.loadChannelOutput(iPaxChannel,ii);
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
detectedNAProc = MD.getProcess(MD.getProcessIndex('PointSourceDetectionProcess'));%MD.getPackage(iUTrack).getProcess(1);%
detectedNAs = detectedNAProc.loadChannelOutput(iPaxChannel);

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
end

% re-express tracksNA so that each track has information for every frame
tic
if ~isempty(iSDCProc)
    disp('reformating NA tracks and applying stage drift correction ...')
    tracksNA = formatNATracks(tracksNAorg,detectedNAs,nFrames,T); 
else
    disp('reformating NA tracks...')
    tracksNA = formatNATracks(tracksNAorg,detectedNAs,nFrames); 
end
toc

% Filter out tracks that is out of XI and YI
xmin = reg_grid(2+band,2+band,1);
ymin = reg_grid(2+band,2+band,2);
xmax = reg_grid(end-band-1,end-band-1,1);
ymax = reg_grid(end-band-1,end-band-1,2);

idxTracks = true(numel(tracksNA),1);
disp('filtering with TFM boundary...')
tic
for ii=1:numel(tracksNA)
%     for k=find(tracksNA(ii).presence)
%         if tracksNA(ii).xCoord(k)<xmin || tracksNA(ii).xCoord(k)>xmax ...
%                 || tracksNA(ii).yCoord(k)<ymin || tracksNA(ii).yCoord(k)>ymax
%             idxTracks(ii) = false;
%             continue;
%         end
%     end
    if any(round(tracksNA(ii).xCoord)<=xmin | round(tracksNA(ii).xCoord)>=xmax ...
            | round(tracksNA(ii).yCoord)<=ymin | round(tracksNA(ii).yCoord)>=ymax)
        idxTracks(ii) = false;
    end
end

toc
tracksNA=tracksNA(idxTracks);
trackIdx = true(numel(tracksNA),1);

% disp('loading segmented FAs...')
FASegProc = FASegPackage.processes_{iFASeg};
regSpacing = (reg_grid(2,1,1)-reg_grid(1,1,1));
for ii=1:nFrames
    [XI,YI]=meshgrid(reg_grid(1,1,1):reg_grid(end,end,1),reg_grid(1,1,2):reg_grid(end,end,2));
    reg_gridFine(:,:,1) = XI;
    reg_gridFine(:,:,2) = YI;
    [~,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_gridFine);
    [grid_mat,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_gridFine);

    pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 
    [grid_mat,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
        
    tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

    % Boundary cutting - I'll take care of this boundary effect later
    grid_mat(end-band*regSpacing:end,:,:)=[];
    grid_mat(:,end-band*regSpacing:end,:)=[];
    grid_mat(1:1+band*regSpacing,:,:)=[];
    grid_mat(:,1:1+band*regSpacing,:)=[];
    tnorm(end-band*regSpacing:end,:,:)=[];
    tnorm(:,end-band*regSpacing:end,:)=[];
    tnorm(1:1+band*regSpacing,:,:)=[];
    tnorm(:,1:1+band*regSpacing,:)=[];
    
    % size of the region of interest
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);

    % Making traction stress map, tsMap
%     [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
%     tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'linear');
    tsMap = tnorm;
    
    % Cell Boundary Mask 
%     maskProc = SegmentationPackage.processes_{2};
%     mask = maskProc.loadChannelOutput(2,ii);
%     cropMask = mask(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    % Cell Boundary

    %%-------------Adhesion detection----------------------
    % loading paxillin image
    if ~isempty(iSDCProc)
        paxImage=(SDCProc.loadChannelOutput(iChan,ii)); %movieData.channels_(2).loadImage(ii);
    else
        paxImage=(MD.channels_(2).loadImage(ii)); 
    end
    paxImageCropped = paxImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);

    % Get Cell mask
%     pId = double(paxImageCropped)/max(double(paxImageCropped(:)));
    % alpha = graythresh(pId);
    %estimate the intensity level to use for thresholding the image
%     level1 = graythresh(pId); %Otsu
%     if ii==1
%         alpha1 = 0.95*level1;
%         alpha2 = 0.75*level1;
%         figure,subplot(2,1,1),imshow(pId,[]),title('original image');
%         subplot(2,2,3),imshow(im2bw(pId,alpha1)),title(['0.95*level1, alpha = ' num2str(alpha1) ]);
%         subplot(2,2,4),imshow(im2bw(pId,alpha2)),title(['0.75*level1, alpha = ' num2str(alpha2) ]);
%         alpha = input('type desired alpha to threshold the image so that it encompass entire cell membrane: ');
%     end
%     pId = double(paxImage)/max(double(paxImageCropped(:)));
%     pId2 = filterGauss2D(pId, 3);
%     bwPI4 = im2bw(pId2,alpha);
    bwPI4 = maskProc.loadChannelOutput(iChan,ii);
    
    % Get the mask for FAs
%     maskAdhesion = blobSegmentThreshold(paxImageCropped,minSize,false,bandMask & cropMask);
    if FASegProc.checkChannelOutput(1) && FASegProc.checkChannelOutput(2)
        iPaxChannel_adh=2;
    else 
        iPaxChannel_adh=iPaxChannel;
    end

    maskFAs = FASegProc.loadChannelOutput(iPaxChannel_adh,ii);
    % Apply stage drift correction
    % Get limits of transformation array
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
    Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
    % Apply subpixel-wise registration to original masks
    I = padarray(maskFAs, [maxY, maxX]);
    maskFAs = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    I = padarray(bwPI4, [maxY, maxX]);
    bwPI4 = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    cropMask = bwPI4(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    [B,~,nBD]  = bwboundaries(cropMask,'noholes');

    % mask for band from edge
    iMask = imcomplement(cropMask);
    distFromEdge = bwdist(iMask);
%     bandwidth_pix = round(bandwidth*1000/MD.pixelSize_);
%     bandMask = distFromEdge <= bandwidth_pix;
    
    maskAdhesion = maskFAs>0;
    % erode once and separate some dumbel shaped adhesion into two
    maskAdhesion2 = bwmorph(maskAdhesion,'hbreak',1);
    % Somehow  this doesn't work. I need find out better way of showing
    % adhesions.
    
    maskAdhesion = maskAdhesion2(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    maskFAs = maskFAs(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
%     maskFAs = maskFAs(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    
%     % Nascent adhesion detection
%     sigmaPSF_NA = getGaussianPSFsigma(MD.numAperture_, 1, MD.pixelSize_*1e-9, 601*1e-9);

%     bandwidthNA = 3; %um after FAKi
    bandwidthNA = 5; %um before FAKi
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    bandMask = distFromEdge <= bandwidthNA_pix;

    bandwidthFA = 1.5; %um before FAKi
    bandwidthFA_pix = round(bandwidthFA*1000/MD.pixelSize_);
    bandMaskFA = distFromEdge <= bandwidthFA_pix;
    
%     naMask = bandMask & cropMask & ~maskAdhesion;
    maskOnlyBand = bandMask & cropMask;
    
    % filter tracks with naMasks
    disp(['Processing ' num2str(ii) 'th frame out of ' num2str(nFrames) ' total frames...'])
    % only deal with presence and status
    for k=1:numel(tracksNA)
        if tracksNA(k).presence(ii) && ~isnan(tracksNA(k).yCoord(ii)) && ...
                ~maskOnlyBand(round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2),round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1))
            tracksNA(k).state{ii} = 'Out_of_ROI';
            tracksNA(k).presence(ii) = false;
            if trackIdx(k)
                trackIdx(k) = false;
            end
        end
    end
   
    % focal contact (FC) analysis
    conn=4;
    CC = bwconncomp(maskAdhesion,conn);
    Adhs = regionprops(CC,'Area','Eccentricity','PixelIdxList','PixelList','Centroid' );
    
    %filter out adhesions that are in 1 um band along the edge: they might
    %be connected nascent adhesions - SH 20140708
    indAdh = true(numel(Adhs),1);
    for k=1:numel(Adhs)
%         plot(Adhs(k).Centroid(1),Adhs(k).Centroid(2),'mo')
        if bandMaskFA(round(Adhs(k).Centroid(2)),round(Adhs(k).Centroid(1)))
            indAdh(k,1) = false;
            maskAdhesion(Adhs(k).PixelIdxList) = false;
        end
    end
    Adhs = Adhs(indAdh);

%     propFAs = regionprops(maskFAs,'Area','Eccentricity','PixelIdxList','PixelList' );
    minFASize = round((1000/MD.pixelSize_)*(1000/MD.pixelSize_)); %adhesion limit=1 um2
    minFCSize = round((600/MD.pixelSize_)*(400/MD.pixelSize_)); %adhesion limit=.24 um2

    fcIdx = arrayfun(@(x) x.Area<minFASize & x.Area>minFCSize, Adhs);
    FCs = Adhs(fcIdx);
    FCForce = arrayfun(@(x) tsMap(x.PixelIdxList),FCs,'UniformOutput',false);
    forceFC(ii) = mean(cell2mat(FCForce));
    FCIdx = find(fcIdx);
    adhBound = bwboundaries(maskAdhesion,conn,'noholes');    
    
    % for larger adhesions
    faIdx = arrayfun(@(x) x.Area>=minFASize, Adhs);
    FAs = Adhs(faIdx);
    FAForce = arrayfun(@(x) tsMap(x.PixelIdxList),FAs,'UniformOutput',false);
    forceFA(ii) = mean(cell2mat(FAForce));
    FAIdx =  find(faIdx);
    neighPix = 2;

    % Reading traction force at each track location
    for k=1:numel(tracksNA)
        if tracksNA(k).presence(ii)
            ynmin = max(1,round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2)-neighPix);
            ynmax = min(size(tsMap,1),round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2)+neighPix);
            xnmin = max(1,round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1)-neighPix);
            xnmax = min(size(tsMap,2),round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1)+neighPix);
            forceNeigh = tsMap(ynmin:ynmax,xnmin:xnmax);
            tracksNA(k).forceMag(ii) = max(forceNeigh(:));    
            if ~strcmp(tracksNA(k).state{ii} , 'NA')
                tracksNA(k).state{ii} = tracksNA(k).state{ii-1};
            end
            % decide if each track is associated with FC or FA
            p = 0;
            for jj=FCIdx'
                p=p+1;
                if any(round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1)==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2)==Adhs(jj).PixelList(:,2))
                    tracksNA(k).state{ii} = 'FC';
%                     tracksNA(k).forceMag(ii) = forceFC(p);
                    tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
                    tracksNA(k).FApixelList{ii} = Adhs(jj).PixelList;
                    tracksNA(k).adhBoundary{ii} = adhBound{jj};
                    tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2),round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1));
                end
            end
            p = 0;
            for jj=FAIdx'
                p=p+1;
                if any(round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1)==Adhs(jj).PixelList(:,1) & round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2)==Adhs(jj).PixelList(:,2))
                    tracksNA(k).state{ii} = 'FA';
%                     tracksNA(k).forceMag(ii) = forceFA(p);
                    tracksNA(k).area(ii) = Adhs(jj).Area;% in pixel
                    tracksNA(k).FApixelList{ii} = Adhs(jj).PixelList;
                    tracksNA(k).adhBoundary{ii} = adhBound{jj};
                    tracksNA(k).faID(ii) = maskFAs(round(tracksNA(k).yCoord(ii))-grid_mat(1,1,2),round(tracksNA(k).xCoord(ii))-grid_mat(1,1,1));
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
                    minDist(q) = min(sqrt((propSubMaskFAs(q).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame) - grid_mat(1,1,1))).^2 +...
                        (propSubMaskFAs(q).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame) - grid_mat(1,1,2))).^2));
                end
                % find the closest segment
                [~,subMaskFAsIdx] = min(minDist);
                subAdhBound = bwboundaries(subMaskFAs,'noholes');    
                [~,closestPixelID] = min(sqrt((propSubMaskFAs(subMaskFAsIdx).PixelList(:,1)-(tracksNA(k).xCoord(tracksNA(k).endingFrame) - grid_mat(1,1,1))).^2 +...
                    (propSubMaskFAs(subMaskFAsIdx).PixelList(:,2)-(tracksNA(k).yCoord(tracksNA(k).endingFrame) - grid_mat(1,1,2))).^2));

                tracksNA(k).state{ii} = 'FC';
                tracksNA(k).xCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,1);
                tracksNA(k).yCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,2);
                tracksNA(k).FApixelList{ii} = propSubMaskFAs(subMaskFAsIdx).PixelList;
                tracksNA(k).adhBoundary{ii} = subAdhBound{subMaskFAsIdx};
            end
        else
            tracksNA(k).forceMag(ii) = NaN;
        end
    end
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
    
    % Showing for debug (TFM)
    if showAllTracks
        h1 = figure('color','w');
        set(h1, 'Position', [100 50 (imSizeX+1)*1.25 imSizeY+1])
        ax1 = subplot('Position',[0 0 0.8 1]);
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
        quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), tmat_vecx./forceScale,tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);

%         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        for k = FCIdx'
    %         eachFA = maskFAs==k;
    %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
    %         for kk=1:nEachFA
    %             adhBoundary = adhBound{kk};
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
        end

        for k = FAIdx'
    %         eachFA = maskFAs==k;
    %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
    %         for kk=1:nEachFA
    %             adhBoundary = adhBound{kk};
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1),  'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if strcmp(tracksNA(k).state{ii} , 'NA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'r', 'LineWidth', 0.5)
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksNA(k).state{ii} , 'FC')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
%                     plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksNA(k).state{ii} , 'FA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) 
%                     plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ko','MarkerSize',markerSize, 'LineWidth', 0.5)
                end
            end
        end
        % Scale bar 2000nm
        subplot('Position',[0.8 0.1 0.1 0.8])
        axis tight
        caxis([tmin tmax]), axis off
        hc = colorbar('West');
        set(hc,'Fontsize',16)
        hold on;

        syFigureStyle(h1,ax1,imSizeY+1);
        print('-dtiff', '-r150', strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat),'.tif'));
    %     hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/forcePeakFig',num2str(ii,iiformat)),'-v7.3')
        print('-depsc2', '-r150', strcat(epsPath,'/forcePeak',num2str(ii,iiformat),'.eps'));
    end    
    if showAllTracks % for pax image
        h2=figure;
        set(h2, 'Position', [100 50+round(1.4*imSizeY) (imSizeX+1) imSizeY+1])

        %Scale bar 2 um
    %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
        paxImageCroppedInverted = imcomplement(paxImageCropped);
        minPax = min(paxImageCroppedInverted(:));
        maxPax = max(paxImageCroppedInverted(:));

        if ii==1
%             minPax1 = 1*minPax;
%             minPax2 = uint16(double(minPax)+double(0.25*(maxPax-minPax)));
            minPax = uint16(double(minPax)+double(0.25*(maxPax-minPax)));
%             hPaxTemp = figure;
%             subplot(1,2,1),imshow(paxImageCroppedInverted,[minPax1 maxPax]),title(['minPax1 = ' num2str(minPax1) ]);
%             subplot(1,2,2),imshow(paxImageCroppedInverted,[minPax2 maxPax]),title(['minPax2 = ' num2str(minPax2) ]);
%             minPax = input('type desired minPax for maximum of the image: ');
%             close(hPaxTemp);
        end        
        imshow(paxImageCroppedInverted,[minPax maxPax]), hold on
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
            plot(adhBoundary(:,2), adhBoundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if strcmp(tracksNA(k).state{ii} , 'NA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'r', 'LineWidth', 0.5)
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ro','MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksNA(k).state{ii} , 'FC')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
%                     plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'o','Color',[255/255 153/255 51/255],'MarkerSize',markerSize, 'LineWidth', 0.5)
                elseif strcmp(tracksNA(k).state{ii} , 'FA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) 
%                     plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ko','MarkerSize',markerSize, 'LineWidth', 0.5)
                end
            end
        end
        syFigureStyle(h2,gca,imSizeY+1);

        print('-depsc2', '-r150', strcat(epsPath,'/pax',num2str(ii,iiformat),'.eps'));
        print('-dtiff', '-r150', strcat(paxtifPath,'/pax',num2str(ii,iiformat),'.tif'));
    %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
        close(h1)
        clear h1
        close(h2)
        clear h2
    end
    imwrite(uint16(round(tsMap*2^15/3500)),strcat(forcemapPath,'/force',num2str(ii,iiformat),'max',num2str(tmax),'.tif'));
    imwrite(paxImageCropped,strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
    if ii==1
        cropPosition = [grid_mat(1,1,1) grid_mat(1,1,2) imSizeX imSizeY];
        save(strcat(dataPath,'/cropInfo.mat'),'cropPosition');
    end
end
% get rid of tracks that have out of rois...
tracksNA = tracksNA(trackIdx);
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
        wRoi = min(tracksNA(k).xCoord(iSF)-grid_mat(1,1,1),r1)...
            +min(imSizeX+1-(tracksNA(k).xCoord(iSF)-grid_mat(1,1,1)),r1);
        hRoi = min(tracksNA(k).yCoord(iSF)-grid_mat(1,1,2),r1)...
            +min(imSizeY+1-(tracksNA(k).yCoord(iSF)-grid_mat(1,1,2)),r1);
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
            xminROI = round(max(1,tracksNA(k).xCoord(iSF)-grid_mat(1,1,1)-(r1-1))); 
            xmaxROI = round(min(imSizeX+1,tracksNA(k).xCoord(iSF)-grid_mat(1,1,1)+r1)); 
            yminROI = round(max(1,tracksNA(k).yCoord(iSF)-grid_mat(1,1,2)-(r1-1)));
            ymaxROI = round(min(imSizeY+1,tracksNA(k).yCoord(iSF)-grid_mat(1,1,2)+r1));
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
                plot(ha1,xFirst-grid_mat(1,1,1)-xminROI,yFirst-grid_mat(1,1,2)-yminROI,'g', 'LineWidth', 0.5)
                % remembering adhesion intensity and TF
                ynmin = max(1,round(yFirst)-grid_mat(1,1,2)-neighPix);
                ynmax = min(size(tsMap,1),round(yFirst)-grid_mat(1,1,2)+neighPix);
                xnmin = max(1,round(xFirst)-grid_mat(1,1,1)-neighPix);
                xnmax = min(size(tsMap,2),round(xFirst)-grid_mat(1,1,1)+neighPix);
                forceNeigh = tsMap(ynmin:ynmax,xnmin:xnmax);
                tracksNA(k).forceMag(j) = max(forceNeigh(:));    
                %here, now I try to do gaussian fit that
                %pointsourceDetection used
                paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                pstruct = fitGaussians2D(double(paxImageCropped), xFirst-grid_mat(1,1,1), yFirst-grid_mat(1,1,2), ampFirst*0.1, sigmaFirst*0.5, min(paxNeigh(:)),'xyAc');
                tracksNA(k).amp(j) = pstruct.A(1);
                if isnan(pstruct.A)
                    tracksNA(k).amp(j) = double(paxImageCropped(round(yFirst)-grid_mat(1,1,2),round(xFirst)-grid_mat(1,1,1)))-min(paxNeigh(:));
                end
            elseif strcmp(tracksNA(k).state{j} , 'NA')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'r', 'LineWidth', 0.5)
            elseif strcmp(tracksNA(k).state{j} , 'FC')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
            elseif strcmp(tracksNA(k).state{j} , 'FA')
                % drawing tracks
                plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
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
                plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'r', 'LineWidth', 0.5)
            elseif strcmp(tracksNA(k).state{j}, 'FC')
                if j==iSF
                    continue
                end
                % drawing tracks
                plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
            elseif strcmp(tracksNA(k).state{j}, 'FA')
                if j==iSF
                    continue
                end
                % drawing tracks
                plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
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
                    ynmin = max(1,round(tracksNA(k).yCoord(fend))-grid_mat(1,1,2)-neighPix);
                    ynmax = min(size(tsMap,1),round(tracksNA(k).yCoord(fend))-grid_mat(1,1,2)+neighPix);
                    xnmin = max(1,round(tracksNA(k).xCoord(fend))-grid_mat(1,1,1)-neighPix);
                    xnmax = min(size(tsMap,2),round(tracksNA(k).xCoord(fend))-grid_mat(1,1,1)+neighPix);
                    forceNeigh = tsMap(ynmin:ynmax,xnmin:xnmax);
                    tracksNA(k).forceMag(j) = max(forceNeigh(:));    
                   
                    paxImageCropped = imread(strcat(paxPath,'/pax',num2str(j,iiformat),'.tif'));
                    xminROI = round(max(1,tracksNA(k).xCoord(iSF)-grid_mat(1,1,1)-(r1-1))); 
                    xmaxROI = round(min(imSizeX+1,tracksNA(k).xCoord(iSF)-grid_mat(1,1,1)+r1)); 
                    yminROI = round(max(1,tracksNA(k).yCoord(iSF)-grid_mat(1,1,2)-(r1-1)));
                    ymaxROI = round(min(imSizeY+1,tracksNA(k).yCoord(iSF)-grid_mat(1,1,2)+r1));

                    paxNeigh = double(paxImageCropped(ynmin:ynmax,xnmin:xnmax));
                    xLast = tracksNA(k).xCoord(j-1);
                    yLast = tracksNA(k).yCoord(j-1);
                    ampLast = tracksNA(k).amp(j-1);
                    sigmaLast = 2.1; %tracksNA(k).sigma(j-1);
                    pstruct = fitGaussians2D(double(paxImageCropped), xLast-grid_mat(1,1,1), yLast-grid_mat(1,1,2), ampLast, sigmaLast*0.5, min(paxNeigh(:)),'xyAc');
                    tracksNA(k).amp(j) = pstruct.A(1);
                    tracksNA(k).xCoord(j) = pstruct.x(1)+grid_mat(1,1,1);
                    tracksNA(k).yCoord(j) = pstruct.y(1)+grid_mat(1,1,2);
                    if isnan(pstruct.A)
                        tracksNA(k).amp(j) = double(paxImageCropped(round(yLast)-grid_mat(1,1,2),round(xLast)-grid_mat(1,1,1)))-min(paxNeigh(:));
                        tracksNA(k).xCoord(j) = tracksNA(k).xCoord(j-1);
                        tracksNA(k).yCoord(j) = tracksNA(k).yCoord(j-1);
                    end
                    
                    paxImageCropped2 = paxImageCropped(yminROI:ymaxROI,xminROI:xmaxROI);
                    tsMapCropped = tsMap(yminROI:ymaxROI,xminROI:xmaxROI);
                    ha1 = subplot('position',[0  0.5  1  0.5]);
                    imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1),colormap(ha1,'gray');freezeColors; hold(ha1,'on')
                    if strcmp(tracksNA(k).state{fend} , 'NA')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'r', 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{fend} , 'FC')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    end
                    ha2 = subplot('position',[0  0  1 0.5]);
                    imshow(tsMapCropped,[tmin2 tmax2],'Parent', ha2), colormap(ha2,'jet');freezeColors; hold(ha2,'on')
                    if strcmp(tracksNA(k).state{fend} , 'NA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'r', 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{fend} , 'FC')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'Color',[255/255 153/255 51/255], 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
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
save(strcat(dataPath,'/cropInfo.mat'),'grid_mat','imSizeX','imSizeY');
save(strcat(dataPath,'/tracks.mat'),'tracksNA');
end

%% formatTracks functions
function newTracks = formatNATracks(tracks,detectedNAs,nFrames,T)
% Format tracks structure into tracks with every frame
if nargin<4
    T=zeros(nFrames,2); % T is a translation matrix
end
% Get limits of transformation array
maxX = ceil(max(abs(T(:, 2))));
maxY = ceil(max(abs(T(:, 1))));


newTracks(numel(tracks),1) = struct('xCoord', [], 'yCoord', [],'state',[],'iFrame',[],'presence',[],'amp',[],'bkgAmp',[]);
% BA: before adhesion, NA: nascent adh, FC: focal complex, FA: focal adh,
% ANA: after NA (failed to be matured.
for i = 1:numel(tracks)
    % Get the x and y coordinate of all compound tracks
    startNA = true;
    endNA = true;
    for  j = 1 : nFrames
        newTracks(i).iFrame(j) = j;
        if j<tracks(i).seqOfEvents(1,1)
            newTracks(i).state{j} = 'BA';
            newTracks(i).xCoord(j) = NaN;
            newTracks(i).yCoord(j) = NaN;
            newTracks(i).presence(j) = false;
            newTracks(i).amp(j) = NaN;
        elseif j>tracks(i).seqOfEvents(2,1)
            newTracks(i).state{j} = 'ANA';
            newTracks(i).xCoord(j) = NaN;
            newTracks(i).yCoord(j) = NaN;
            newTracks(i).amp(j) = NaN;
            newTracks(i).presence(j) = false;
            if endNA
                newTracks(i).endingFrame = j-1;
                endNA = false;
            end
        elseif j==tracks(i).seqOfEvents(2,1)
            newTracks(i).state{j} = 'NA';
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,2)+maxX;
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,1)+maxY;
            newTracks(i).amp(j) = tracks(i).tracksCoordAmpCG(1,4+8*(j-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(j) = NaN;
            else
                newTracks(i).bkgAmp(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).bkg(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).sigmaX(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(j) = true;
            if endNA
                newTracks(i).endingFrame = j;
                endNA = false;
            end
        else
            newTracks(i).state{j} = 'NA';
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,2)+maxX;
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,1)+maxY;
            newTracks(i).amp(j) = tracks(i).tracksCoordAmpCG(1,4+8*(j-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(j) = NaN;
            else
                newTracks(i).bkgAmp(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).bkg(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(j) = detectedNAs(j-tracks(i).seqOfEvents(1,1)+1).sigmaX(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(j) = true;
            if startNA
                newTracks(i).startingFrame = j;
                startNA = false;
            end
        end
            
        if isfield(tracks, 'label'),
            newTracks(iTrack).label = tracks(i).label;
        end
    end
    % go through frames again and fill NaNs with numbers at the gap
    % position
    for j=1:nFrames-1
        if j<nFrames-9 && sum(newTracks(i).presence(j:j+9))==10 ...
                && sum(isnan(newTracks(i).xCoord(j:j+9)))==10 
            gap = 10;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
            end
        elseif j<nFrames-8 && sum(newTracks(i).presence(j:j+8))==9 ...
                && sum(isnan(newTracks(i).xCoord(j:j+8)))==9 
            gap = 9;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
            end
        elseif j<nFrames-7 && sum(newTracks(i).presence(j:j+7))==8 ...
                && sum(isnan(newTracks(i).xCoord(j:j+7)))==8 
            gap = 8;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
            end
        elseif j<nFrames-6 && sum(newTracks(i).presence(j:j+6))==7 ...
                && sum(isnan(newTracks(i).xCoord(j:j+6)))==7 
            gap = 7;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
            end
        elseif j<nFrames-5 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
               && newTracks(i).presence(j+4) && newTracks(i).presence(j+5) && isnan(newTracks(i).xCoord(j)) ...
               && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))...
               && isnan(newTracks(i).xCoord(j+4)) && isnan(newTracks(i).xCoord(j+5))
            newTracks(i).xCoord(j) = (6*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j) = (6*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j) = (6*newTracks(i).amp(j-1)+newTracks(i).amp(j+6))/7;
            newTracks(i).xCoord(j+1) = (5*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+1) = (5*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+1) = (5*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+6))/7;
            newTracks(i).xCoord(j+2) = (4*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+2) = (4*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+2) = (4*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+6))/7;
            newTracks(i).xCoord(j+3) = (3*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+3) = (3*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+3) = (3*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+6))/7;
            newTracks(i).xCoord(j+4) = (2*newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+4) = (2*newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+4) = (2*newTracks(i).amp(j-1)+5*newTracks(i).amp(j+6))/7;
            newTracks(i).xCoord(j+5) = (newTracks(i).xCoord(j-1)+6*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+5) = (newTracks(i).yCoord(j-1)+6*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+5) = (newTracks(i).amp(j-1)+6*newTracks(i).amp(j+6))/7;
        elseif j<nFrames-4 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
                && newTracks(i).presence(j+4) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) ...
                && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3)) && isnan(newTracks(i).xCoord(j+4))
            newTracks(i).xCoord(j) = (5*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j) = (5*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j) = (5*newTracks(i).amp(j-1)+newTracks(i).amp(j+5))/6;
            newTracks(i).xCoord(j+1) = (4*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+1) = (4*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+1) = (4*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+5))/6;
            newTracks(i).xCoord(j+2) = (3*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+2) = (3*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+2) = (3*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+5))/6;
            newTracks(i).xCoord(j+3) = (2*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+3) = (2*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+3) = (2*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+5))/6;
            newTracks(i).xCoord(j+4) = (newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+4) = (newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+4) = (newTracks(i).amp(j-1)+5*newTracks(i).amp(j+5))/6;
        elseif j<nFrames-3 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
                && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))
            newTracks(i).xCoord(j) = (4*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j) = (4*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j) = (4*newTracks(i).amp(j-1)+newTracks(i).amp(j+4))/5;
            newTracks(i).xCoord(j+1) = (3*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+1) = (3*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+1) = (3*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+4))/5;
            newTracks(i).xCoord(j+2) = (2*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+2) = (2*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+2) = (2*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+4))/5;
            newTracks(i).xCoord(j+3) = (newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+3) = (newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+3) = (newTracks(i).amp(j-1)+4*newTracks(i).amp(j+4))/5;
        elseif j<nFrames-2 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) ...
                && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2))
            newTracks(i).xCoord(j) = (3*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j) = (3*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j) = (3*newTracks(i).amp(j-1)+newTracks(i).amp(j+3))/4;
            newTracks(i).xCoord(j+1) = (2*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j+1) = (2*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j+1) = (2*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+3))/4;
            newTracks(i).xCoord(j+2) = (newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j+2) = (newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j+2) = (newTracks(i).amp(j-1)+3*newTracks(i).amp(j+3))/4;
        elseif j<nFrames-1 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1))
            newTracks(i).xCoord(j) = (2*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+2))/3;
            newTracks(i).yCoord(j) = (2*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+2))/3;
            newTracks(i).amp(j) = (2*newTracks(i).amp(j-1)+newTracks(i).amp(j+2))/3;
            newTracks(i).xCoord(j+1) = (newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+2))/3;
            newTracks(i).yCoord(j+1) = (newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+2))/3;
            newTracks(i).amp(j+1) = (newTracks(i).amp(j-1)+2*newTracks(i).amp(j+2))/3;
        elseif newTracks(i).presence(j) && isnan(newTracks(i).xCoord(j))
            newTracks(i).xCoord(j) = (newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+1))/2;
            newTracks(i).yCoord(j) = (newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+1))/2;
            newTracks(i).amp(j) = (newTracks(i).amp(j-1)+newTracks(i).amp(j+1))/2;
        end
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