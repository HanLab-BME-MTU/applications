function [tracksNAfailing,tracksNAmaturing,lifeTimeNAfailing,lifeTimeNAmaturing,maturingRatio] = analyzeAdhesionMaturation(pathForTheMovieDataFile,outputPath,showAllTracks,plotEachTrack)
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

% Sangyoon Han April 2014

%% Data Set up
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
MD = MovieData.load(movieDataPath);
% Get whole frame number
nFrames = MD.nFrames_;
% Get U-Track package (actually NA package)
NAPackage = MD.getPackage(MD.getPackageIndex('UTrackPackage'));
% Get FA segmentation package 
FASegPackage = MD.getPackage(MD.getPackageIndex('FocalAdhesionSegmentationPackage'));
% Load tracks
iTracking = 2;
trackNAProc = NAPackage.processes_{iTracking};

% Load the Paxillin channel

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
minLifetime = min(nFrames,3);
markerSize = 2;
% tracks
iPaxChannel = 1; % this should be intentionally done in the analysis level

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
detectedNAProc = NAPackage.processes_{1};
detectedNAs = detectedNAProc.loadChannelOutput(iPaxChannel);

% re-express tracksNA so that each track has information for every frame
disp('reformating NA tracks...')
tic
tracksNA = formatTracks(tracksNAorg,detectedNAs,nFrames); 
toc

trackIdx = true(numel(tracksNA),1);

% disp('loading segmented FAs...')
iFASeg = 6;
FASegProc = FASegPackage.processes_{iFASeg};
for ii=1:nFrames
    % Cell Boundary Mask 
    maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
    mask = maskProc.loadChannelOutput(1,ii);
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');

    %%-------------Adhesion detection----------------------

    % mask for band from edge
    iMask = imcomplement(mask);
    distFromEdge = bwdist(iMask);

    % Get the mask for FAs
    maskFAs = FASegProc.loadChannelOutput(iPaxChannel,ii);
    maskAdhesion = maskFAs>0;
    
    bandwidthNA = 5; %um 
    bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
    bandMask = distFromEdge <= bandwidthNA_pix;

    maskOnlyBand = bandMask & mask;
    paxImageCropped=MD.channels_(1).loadImage(ii); 
    % size of the region of interest
    if ii==1
        [imSizeY,imSizeX] = size(paxImageCropped);
    end

    % filter tracks with naMasks
    disp(['Processing ' num2str(ii) 'th frame out of ' num2str(nFrames) ' total frames...'])
    % only deal with presence and status
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
   
    % focal contact (FC) analysis
    Adhs = regionprops(maskAdhesion,'Area','Eccentricity','PixelIdxList','PixelList' );
%     propFAs = regionprops(maskFAs,'Area','Eccentricity','PixelIdxList','PixelList' );
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
                tracksNA(k).xCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,1);
                tracksNA(k).yCoord(ii) = propSubMaskFAs(subMaskFAsIdx).PixelList(closestPixelID,2);
                tracksNA(k).FApixelList{ii} = propSubMaskFAs(subMaskFAsIdx).PixelList;
                tracksNA(k).adhBoundary{ii} = subAdhBound{subMaskFAsIdx};
            end
        end
    end
    
    if showAllTracks
        h2=figure;
        %Scale bar 2 um
    %     paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
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
    imwrite(paxImageCropped,strcat(paxPath,'/pax',num2str(ii,iiformat),'.tif'));
end
% get rid of tracks that have out of bands...
tracksNA = tracksNA(trackIdx);
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
tracksNAfailing = trNAonly(indMature);
tracksNAmaturing = trNAonly(indFail);

if plotEachTrack
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
            imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1),colormap(ha1,'gray');freezeColors; hold(ha1,'on')
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
                    imshow(imcomplement(paxImageCropped2),[pmin pmax],'Parent', ha1),colormap(ha1,'gray');freezeColors; hold(ha1,'on')
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