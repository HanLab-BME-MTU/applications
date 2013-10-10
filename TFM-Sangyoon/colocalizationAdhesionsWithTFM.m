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
% Get FA package (actually NA package)
NAPackage = MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Get FA segmentation package 
FASegPackage = MD.getPackage(MD.getPackageIndex('FocalAdhesionSegmentationPackage'));
% Load tracks
iTracking = 4;
trackNAProc = NAPackage.processes_{iTracking};

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

% Load the Paxillin channel

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
% SegmentationPackage = MD.getPackage(MD.getPackageIndex('SegmentationPackage'));
% minSize = round((500/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=.5um*.5um
minLifetime = 1;

% tracks
iPaxChannel = 2; % this should be intentionally done in the analysis level

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
detectedNAProc = NAPackage.processes_{3};
detectedNAs = detectedNAProc.loadChannelOutput(iPaxChannel);

% re-express tracksNA so that each track has information for every frame
disp('reformating NA tracks...')
tic
tracksNA = formatTracks(tracksNAorg,detectedNAs,nFrames); 
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
iFASeg = 6;
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
    iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(iSDCProc)
        SDCProc=MD.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(1)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        paxImage=(SDCProc.loadChannelOutput(2,ii)); %movieData.channels_(2).loadImage(ii);
    else
        paxImage=(MD.channels_(2).loadImage(ii)); 
    end
    paxImageCropped = paxImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);

    % Get Cell mask
    pId = double(paxImageCropped)/max(double(paxImageCropped(:)));
    % alpha = graythresh(pId);
    %estimate the intensity level to use for thresholding the image
    level1 = graythresh(pId); %Otsu
%     [~, level2] = cutFirstHistMode(pId,0); %Rosin
    alpha = 0.94*level1;
%     alpha = 0.82*level1; for after fak
%     alpha = 0.9*level2 + 0.87*level1;
%     tempH = figure; subplot(2,1,1), imshow(paxImageCropped,[]), title('original image');
%     subplot(2,2,3), imshow(im2bw(pId,level1)), title(['Otsu, alpha= ', num2str(level1)]);
%     subplot(2,2,4), imshow(im2bw(pId,level2)), title(['Rosin, alpha= ', num2str(level2)]);
%     
%     alpha = input(['Type desired alpha: (recommended: ' num2str(0.9*level2) ') :']);
%     close(tempH)
%     clear tempH

    pId = double(paxImage)/max(double(paxImageCropped(:)));
    pId2 = filterGauss2D(pId, 3);
    bwPI4 = im2bw(pId2,alpha);
%     bwPI = im2bw(pId,alpha);
%     bwPI2 = bwmorph(bwPI,'clean');
%     bwPI3 = bwmorph(bwPI2,'erode',5);
%     bwPI3 = bwmorph(bwPI3,'dilate',4);
%     bwPI4 = bwmorph(bwPI3,'close',10);
%     bwPI4 = bwmorph(bwPI4,'dilate',1);
    % bwPI5 = refineEdgeWithSteerableFilterGM(pId,bwPI4);
    % In case that there is still islands, pick only the largest chunk
%     [labelPI,nChunk] = bwlabel(bwPI4);
%     if nChunk>1
%         eachArea = zeros(nChunk,1);
%         for k=1:nChunk
%             currBWPI = labelPI==k;
%             eachArea(k) = sum(currBWPI(:));
%         end
%         [~,indCellArea] = max(eachArea);
%         bwPI4 = labelPI == indCellArea;
%     end
    cropMask = bwPI4(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    [B,~,nBD]  = bwboundaries(cropMask,'noholes');

    % mask for band from edge
    iMask = imcomplement(cropMask);
    distFromEdge = bwdist(iMask);
%     bandwidth_pix = round(bandwidth*1000/MD.pixelSize_);
%     bandMask = distFromEdge <= bandwidth_pix;

    % Get the mask for FAs
%     maskAdhesion = blobSegmentThreshold(paxImageCropped,minSize,false,bandMask & cropMask);
    maskFAs = FASegProc.loadChannelOutput(iPaxChannel,ii);
    maskAdhesion = maskFAs>0;
    maskAdhesion = maskAdhesion(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    maskFAs = maskFAs(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
%     maskFAs = maskFAs(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    
%     % Nascent adhesion detection
%     sigmaPSF_NA = getGaussianPSFsigma(MD.numAperture_, 1, MD.pixelSize_*1e-9, 601*1e-9);

%     bandwidthNA = 3; %um after FAKi
    bandwidthNA = 5; %um before FAKi
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    bandMask = distFromEdge <= bandwidthNA_pix;

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
   
%     % Get all nascent adhesions
%     if bFirst
%         pstruct_NA = pointSourceDetection(paxImageCropped, sigmaPSF_NA*1.3, 'alpha', 0.01,'Mask',naMask);
% 
%         % Showing for debug
% %         h0=figure; imshow(paxImageCropped,[])
% %         hold on;
% %         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.25) % cell boundary
% %         plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
% %         for k = 1:length(adhBound)
% %             adhBoundary = adhBound{k};
% %             plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1) %adhesion boundary
% %         end
% 
% %         % ginput test for not-detected adhesions
% %         disp('Click nascent adhesions and press Enter when you are done with clicking');
% %         NA_man_detected = ginput;
% %         if ~isempty(NA_man_detected)
% %             plot(NA_man_detected(:,1),NA_man_detected(:,2),'ro')
% %             nPoints = length(pstruct_NA.A);
% %             old_pstruct_NA = pstruct_NA;
% %             pstruct_NA.x(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,1);
% %             pstruct_NA.y(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,2);
% %             % for amplitude
% %             for k=1:size(NA_man_detected,1)
% %                 pstruct_NA.A(nPoints+k) = paxImageCropped(round(NA_man_detected(k,2)),round(NA_man_detected(k,1)));
% %             end
% %         end
%         pstruct_NAcell{ii} = pstruct_NA;
%     else
%         pstruct_NA = pstruct_NAcell{ii};
%     end
    % focal contact (FC) analysis
    Adhs = regionprops(maskAdhesion,'Area','Eccentricity','PixelIdxList','PixelList' );
    propFAs = regionprops(maskFAs,'Area','Eccentricity','PixelIdxList','PixelList' );
    minFASize = round((2000/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=1um*.5um
    minFCSize = round((800/MD.pixelSize_)*(300/MD.pixelSize_)); %adhesion limit=1um*.5um

    fcIdx = arrayfun(@(x) x.Area<minFASize & x.Area>minFCSize, Adhs);
    FCs = Adhs(fcIdx);
    FCForce = arrayfun(@(x) tsMap(x.PixelIdxList),FCs,'UniformOutput',false);
    forceFC(ii) = mean(cell2mat(FCForce));
    FCIdx = find(fcIdx);
    adhBound = bwboundaries(maskAdhesion,'noholes');    
    
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
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'b', 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'b', 'LineWidth', 0.5) %adhesion boundary
        end

        for k = FAIdx'
    %         eachFA = maskFAs==k;
    %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
    %         for kk=1:nEachFA
    %             adhBoundary = adhBound{kk};
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if strcmp(tracksNA(k).state{ii} , 'NA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'r')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ro')
                elseif strcmp(tracksNA(k).state{ii} , 'FC')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'b')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'bo')
                elseif strcmp(tracksNA(k).state{ii} , 'FA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'k')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ko')
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

        syFigureStyle(h1,ax1,2);
        print('-dtiff', '-r300', strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat),'.tif'));
    %     hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/forcePeakFig',num2str(ii,iiformat)),'-v7.3')
        print('-depsc2', '-r300', strcat(epsPath,'/forcePeak',num2str(ii,iiformat),'.eps'));
    end    
%     % Find detectable NAs from pstruct_NA
%     NAs = [grid_mat(1,1,1)+pstruct_NA.x' grid_mat(1,1,2)+pstruct_NA.y'];
%     deformedBeads = [displField(ii).pos(:,1)+displField(ii).vec(:,1) displField(ii).pos(:,2)+displField(ii).vec(:,2)];
%     idx = KDTreeBallQuery(deformedBeads,NAs, 3);
%     valid = ~cellfun(@isempty, idx);
% %     valid = (displField.vec(cell2mat(idx),1).^2+displField.vec(cell2mat(idx),2).^2).^0.5;
%     NAdetec = NAs(valid, :);
%     NAundetec = NAs(~valid, :);

%     plot(NAdetec(:,1)-grid_mat(1,1,1),NAdetec(:,2)-grid_mat(1,1,2),'ro')
%     plot(NAundetec(:,1)-grid_mat(1,1,1),NAundetec(:,2)-grid_mat(1,1,2),'ro')
%     forceNAdetec = zeros(length(NAdetec),1);
%     for i=1:size(NAdetec,1)
%         curF_NAdetec = tsMap(round(NAdetec(i,2))-grid_mat(1,1,2)-neighPix:round(NAdetec(i,2))-grid_mat(1,1,2)+neighPix,...
%             round(NAdetec(i,1))-grid_mat(1,1,1)-neighPix:round(NAdetec(i,1))-grid_mat(1,1,1)+neighPix);
%         forceNAdetec(i) = max(curF_NAdetec(:));
%     end
%     forceNAund = zeros(length(NAundetec),1);
%     for i=1:size(NAundetec,1)
%         curF_NAund = tsMap(round(NAundetec(i,2))-grid_mat(1,1,2)-neighPix:round(NAundetec(i,2))-grid_mat(1,1,2)+neighPix,...
%             round(NAundetec(i,1))-grid_mat(1,1,1)-neighPix:round(NAundetec(i,1))-grid_mat(1,1,1)+neighPix);
%         forceNAund(i) = max(curF_NAund(:));
%     end
%     forceNA = [forceNAdetec; forceNAund];
%     forceNAdetec = tsMap(round(NAdetec(:,2))-grid_mat(1,1,2),round(NAdetec(:,1))-grid_mat(1,1,1));
%     plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
%     plot(NA_man_detected(:,1),NA_man_detected(:,2),'m*')
        
    
%     idxCurTracksNA = arrayfun(@(x) x.seqOfEvents(1,1)<=ii & x.seqOfEvents(2,1)>=ii,tracksNA); % may not need this
    
    
%     % Get magnitude of force and the curvature of the force at each (stored
%     % in pstruct_NAwithForce.fmag and pstruct_NAwithForce.fcurvature)
%     pstruct_NAwithForce = findMagCurvature(tsMap,pstruct_NA,5);
%     
%     % Plot paxillin intensity with force magnitude
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fmag,'.')
%     xlabel('Pax Intensity (A.U.)')
%     ylabel('Stress magnitude (Pa)')
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fcurvature,'r.')
%     % Make a histogram with pstruct_NAwithForce
%     figure, hist(pstruct_NAwithForce.fmag,250:500:3000)

%     save(strcat(dataPath,'/pstructNAForce',num2str(ii,iiformat)), 'pstruct_NAwithForce');
%     %---- for FAs, do line integral of force gradient along focal adhesion
%     distVals=-3:3; %positive =inside
%     repDist = 2.5:-1:-2.5; % making positive=outside for plotting
%     [avgFprofile,avgFprofileStd,nPts] = intensityVsDistFromEdge(tsMap,maskAdhesion,distVals);
%     figure,plot(repDist,avgFprofile)
%     [~,avgFslope,~] = regression(repDist,avgFprofile);
% %     plot(repDist,mean(avgFprofile,1),'k','LineWidth',2)
%     FAstruct.avgFprofile = avgFprofile;
%     FAstruct.avgFprofileStd = avgFprofileStd;
%     FAstruct.avgFslope = avgFslope;
%     FAstruct.Finside = avgFprofile(end);
%     FAstruct.FinsideErr = avgFprofileStd(end)/sqrt(max(nPts)); %SEM
%     save(strcat(dataPath,'/FAstruct',num2str(ii,iiformat)),'FAstruct');
    if showAllTracks
        h2=figure;
        set(h2, 'Position', [100 50+round(1.4*imSizeY) (imSizeX+1) imSizeY+1])

        %Scale bar 2 um
    %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
        imshow(imcomplement(paxImageCropped),[]), hold on
        line([10 10+round(2000/MD.pixelSize_)],[15 15],'LineWidth',2,'Color',[0,0,0])
        
        for kk=1:nBD
            boundary = B{kk};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        end
    %     plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
        for k = FCIdx'
    %         eachFA = maskFAs==k;
    %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
    %         for kk=1:nEachFA
    %             adhBoundary = adhBound{kk};
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'b', 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'b', 'LineWidth', 0.5) %adhesion boundary
        end
        % for larger adhesions
        for k = FAIdx'
    %         eachFA = maskFAs==k;
    %         [adhBound,~,nEachFA] = bwboundaries(eachFA,'noholes');
    %         for kk=1:nEachFA
    %             adhBoundary = adhBound{kk};
    %             plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
    %         end
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 0.5) %adhesion boundary
        end
        for k=1:numel(tracksNA)
            if tracksNA(k).presence(ii)
                if strcmp(tracksNA(k).state{ii} , 'NA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'r')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ro')
                elseif strcmp(tracksNA(k).state{ii} , 'FC')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'b')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'bo')
                elseif strcmp(tracksNA(k).state{ii} , 'FA')
                    % drawing tracks
                    plot(tracksNA(k).xCoord(1:ii)-grid_mat(1,1,1),tracksNA(k).yCoord(1:ii)-grid_mat(1,1,2),'k')
                    plot(tracksNA(k).xCoord(ii)-grid_mat(1,1,1),tracksNA(k).yCoord(ii)-grid_mat(1,1,2),'ko')
                end
            end
        end
        syFigureStyle(h2,gca,2);

        print('-depsc2', '-r300', strcat(epsPath,'/pax',num2str(ii,iiformat),'.eps'));
        print('-dtiff', '-r300', strcat(paxtifPath,'/pax',num2str(ii,iiformat),'.tif'));
    %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
        close(h1)
        clear h1
        close(h2)
        clear h2
    end
    imwrite(uint16(round(tsMap*2^15/3500)),strcat(forcemapPath,'/force',num2str(ii,iiformat),'max',num2str(tmax),'.tif'));
    imwrite(paxImageCropped,strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
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
                plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'b', 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
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
                plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'b', 'LineWidth', 0.5)
                adhBoundary = tracksNA(k).adhBoundary{j};
                plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
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
                        plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'b', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha1,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha1,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
                        end
                    end
                    ha2 = subplot('position',[0  0  1 0.5]);
                    imshow(tsMapCropped,[tmin2 tmax2],'Parent', ha2), colormap(ha2,'jet');freezeColors; hold(ha2,'on')
                    if strcmp(tracksNA(k).state{fend} , 'NA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'r', 'LineWidth', 0.5)
                    elseif strcmp(tracksNA(k).state{fend} , 'FC')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'b', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
                        end
                    elseif strcmp(tracksNA(k).state{fend} , 'FA')
                        % drawing tracks
                        plot(ha2,tracksNA(k).xCoord(1:j)-grid_mat(1,1,1)-xminROI,tracksNA(k).yCoord(1:j)-grid_mat(1,1,2)-yminROI,'k', 'LineWidth', 0.5)
                        adhBoundary = tracksNA(k).adhBoundary{j};
                        if ~isnan(adhBoundary)
                            plot(ha2,adhBoundary(:,2)-xminROI, adhBoundary(:,1)-yminROI, 'b', 'LineWidth', 0.5) %adhesion boundary
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
function newTracks = formatTracks(tracks,detectedNAs,nFrames)
% Format tracks structure into tracks with every frame

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
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)));
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)));
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
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)));
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)));
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