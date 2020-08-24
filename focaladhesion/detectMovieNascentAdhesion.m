function [nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion(pathForTheMovieDataFile,bandwidth,iAdhChan,plotGraph,indMask,plotGraphTFM)
% detectMovieNascentAdhesion detect objects by fitting isotropic Gaussians
%
% SYNOPSIS adapted from detectMovieSubResFeatures(movieData,paramsIn)
%
% INPUT   
%   pathForTheMovieDataFile - path to a MovieData object describing the movie to be processed
%   bandwidth               - bandwidth in micron from the edge to qunatify
%                               nascent adhesions (default: 5 um)
%   iAdhChan                - the channel number for adhesion channel used
%   plotGraph               - true if you want to plot graphs and store
%                             them
%   indMask                 - the channel index of a mask to be added to
%                             mask in iAdhChan channel. If you add it, you are analyzing more
%                             features (e.g. CCPs)
%   plotGraphTFM            - true if you want to plot graphs about TFM
% 
% OUTPUT   
%   adhesionInfo  - stored in Adhesion Quantification folder containing:
%   nascentAdhInfo: locations of nascent adhesions, area of lamellipodial band
%         xCoord,yCoord:       (x,y) of NAs
%         amp:                      intensities
%         numberNA:             the number of detected NAs
%         bandArea:              area of band in um^2
%         NAperArea:            the number of NAs per band area.

%   focalAdhInfo: focal adhesion number, area and length
%         xCoord,yCoord       (x,y) of all FAs
%         area:                     semented area of all FAs
%         length:                   major axis length of all FAs
%         amp:                      mean intensities of all FAs
%         numberFA              the number of FAs
%         meanFAarea           mean area of FAs
%         medianFAarea         median area of FAs
%         meanLength            mean length of FAs;
%         medianLength          median length of FAs;

% Example:  [nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion('/some path/Nadia/WT',8)
% Sangyoon Han
% Last updated April 2014
%% ----------- Input ----------- %%
if nargin==1
    bandwidth = 5;
    plotGraph = true;
    plotGraphTFM = false;
end
if nargin<3
    iAdhChan = 1; %assumed
    plotGraph = true;
    indMask = [];
    plotGraphTFM = false;
end
if nargin<4
    plotGraph = true;
    indMask = [];
    plotGraphTFM = false;
end
if nargin<5
    indMask = []; %indMask is the channel index of a mask to be added to mask in iPax channel.
    plotGraphTFM = false;
end
% Load the MovieData
if isa(pathForTheMovieDataFile, 'MovieData')
    movieData = pathForTheMovieDataFile;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    movieData = MovieData.load(movieDataPath,false);
end

%% --------------- Initialization ---------------%%
% Set up the output directories
outputFilePath = [movieData.outputDirectory_ filesep 'Adhesion Quantification'];
imgPath = [outputFilePath filesep 'imgs'];
dataPath = [outputFilePath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];

if ~exist(figPath,'dir')
    mkClrDir(imgPath);
    mkClrDir(dataPath);
    mkClrDir(tifPath);
    mkClrDir(figPath);
end

nChans=numel(movieData.channels_);
% e.g. see if there is TFM package
iTFM =  movieData.getPackageIndex('TFMPackage');
% e.g. see if the additional channel is used for TFM or not.
%      if there is a channel, get the channel.
if ~isempty(iTFM)
    tfmPack = movieData.getPackage(iTFM);
    forceProc = tfmPack.getProcess(4);
    dispProc = tfmPack.getProcess(2);
    
    dispFunParam = dispProc.funParams_; % forceProc.checkChannelOutput;
    iChanTFM = dispFunParam.ChannelIndex;
    % iAdhChan can be mistakenly input. It is reasonable to change it to
    % the other channel if it overlaps with iChanTFM
    if iAdhChan==iChanTFM && nChans==2
        iAdhChan = setdiff(1:nChans,[iChanTFM]);
        disp(['Adhesion channel index is automatically changed to ' num2str(iAdhChan) ' since ' ...
            num2str(iChanTFM) ' overlapped with it.']);
    end
    iAdditionalChan = setdiff(1:nChans,[iChanTFM, iAdhChan]);
    % Read TFM map
    tMap=forceProc.loadChannelOutput('output','tMapUnshifted'); % in Pa per pixel (1pix x 1pix)
    forceFA=cell(1,movieData.nFrames_);
    forceFC=cell(1,movieData.nFrames_);
    forceNA=cell(1,movieData.nFrames_);
    forceBGinCell=cell(1,movieData.nFrames_);
    forceBGoutCell=cell(1,movieData.nFrames_);
    tMax = 2000;
else
    iAdditionalChan = setdiff(1:nChans,[iAdhChan]);
end

% psfSigma = 1.2; %hard coded for nascent adhesion
psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(iAdhChan).emissionWavelength_*1e-9);
if isempty(psfSigma)
    movieData.getChannel(iAdhChan).emissionWavelength_ = 568; 
    disp('Adhesion channel emmission wave length is set to be 568 nm. If this value is incorrect, please set up your fluorophore of the adhesion channel')
    psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(iAdhChan).emissionWavelength_*1e-9);    
end
%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting isotropic Gaussians...')

roiMask = movieData.getROIMask;
nascentAdhInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'NAdensity',[],'bandMask',[],'numberNA',[],'bandArea',[]);
focalAdhInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'length',[],'meanFAarea',[],'medianFAarea',[]...
    ,'meanLength',[],'medianLength',[],'numberFA',[],'FAdensity',[],'cellArea',[]...
    ,'maskFA',[],'boundFA',[],'ecc',[],'width',[],'labelFA',[],'FAdensityPeri',[],'FAdensityInside',[]);
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
pixSize = movieData.pixelSize_;
minSize = round((500/pixSize)*(500/pixSize)); %adhesion limit=0.25 um2
minLengthFC = 500/pixSize;
minLengthFA = 2000/pixSize;
minEcc = 0.7;
psAlpha = 0.0001;%This will ensure 

disp(['Paxillin channel was assumed to be in channel ' num2str(iAdhChan) '.'])
disp('Results will be saved under:')
disp(outputFilePath);

for j=1:movieData.nFrames_
    I=double(movieData.channels_(iAdhChan).loadImage(j));
    noMask=false;
    try
        maskProc = movieData.getProcess(movieData.getProcessIndex('MaskRefinementProcess'));
%         mask = maskProc.loadChannelOutput(iPax,j);
        % if there are masks for more than one channels, combine them.
        
        if length(maskProc.checkChannelOutput)>1
            %Combine the the multiple masks to one
            maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,j),find(maskProc.checkChannelOutput),'UniformOutput',false);
            maskAll=reshape(cell2mat(maskEach),size(I,1),size(I,2),[]);
            mask = any(maskAll,3);
        elseif nChans==1 && nChans==iAdhChan
            mask = maskProc.loadChannelOutput(iAdhChan,j); % 1 is CCP channel
        end
    catch
        mask = movieData.roiMask;
        noMask=true;
    end
    if ~isempty(indMask)
        maskComp = maskProc.loadChannelOutput(indMask,j); % 1 is CCP channel
        mask = mask | maskComp;
    end
%         maskProc = movieData.processes_{2};
%         mask = maskProc.loadChannelOutput(1,j);
    maskAdhesion = blobSegmentThreshold(I,minSize,0,mask & roiMask(:,:,j));
    maskAdhesionFine = blobSegmentThreshold(I,0); % mask for all adhesions without minSize

%     maskAdhesionC = imcomplement(maskAdhesion);
    % mask for band from edge
    if noMask
        ultimateMask = roiMask(:,:,j) & mask & maskAdhesionFine;
%         ultimateMask = roiMask(:,:,j) & maskAdhesionC & mask & maskAdhesionFine;
    elseif ~isempty(bandwidth)
        iMask = imcomplement(mask);
        distFromEdge = bwdist(iMask);
        bandwidth_pix = round(bandwidth*1000/pixSize);
        bandMask = distFromEdge <= bandwidth_pix;

%         ultimateMask = bandMask & roiMask(:,:,j) & maskAdhesionC & mask & maskAdhesionFine;
        ultimateMask = bandMask & roiMask(:,:,j) & mask & maskAdhesionFine;
    else
%         ultimateMask = roiMask(:,:,j) & maskAdhesionC & mask; % & maskAdhesionFine;
        ultimateMask = roiMask(:,:,j) & mask; % & maskAdhesionFine;
    end
    [pstruct, maskNAs] = pointSourceDetection(I, psfSigma,  ...
        'Alpha',psAlpha,'Mask',ultimateMask);
        % filter out points where overall intensity is in the noise level
%     pixelIntenMargin = I(~mask);
%     maxIntBg=quantile(pixelIntenMargin,0.9999);
%     psInt = pstruct.A+pstruct.c;
%     idxSigCCP = psInt>maxIntBg;
    if ~isempty(pstruct)
        xNA=pstruct.x;
        yNA=pstruct.y;
        maskAdhesion2 = refineAdhesionSegmentation(maskAdhesion,I,xNA,yNA); %,mask);
    else
        xNA=[];
        yNA=[];
        maskAdhesion2 = maskAdhesion;
    end
%     labelAdhesion = bwlabel(maskAdhesion);
    Adhs = regionprops(maskAdhesion2,'Centroid','Area','Eccentricity','PixelIdxList','MajorAxisLength','MinorAxisLength');
%         minFASize = round((2000/movieData.pixelSize_)*(500/movieData.pixelSize_)); %adhesion limit=1um*.5um

    adhEccIdx = arrayfun(@(x) x.Eccentricity>minEcc, Adhs);
    FAlengthAll = arrayfun(@(x) x.MajorAxisLength, Adhs);
%     maxLength=mean(FAlengthAll)+5*std(FAlengthAll);
    adhLengthIdxFC = FAlengthAll>minLengthFC;
    AdhsFCFA = Adhs(adhEccIdx & adhLengthIdxFC);
    labelAdhesion = zeros(size(maskAdhesion2));
    for kk=1:numel(AdhsFCFA)
        labelAdhesion(AdhsFCFA(kk).PixelIdxList)=kk;
    end
    maskAdhesion2 = logical(labelAdhesion);
    adhLengthIdxFA = FAlengthAll>minLengthFA;
    AdhsFA = Adhs(adhEccIdx & adhLengthIdxFA);
    labelAdhesionFA = zeros(size(maskAdhesion2));
    for kk=1:numel(AdhsFA)
        labelAdhesionFA(AdhsFA(kk).PixelIdxList)=kk;
    end
    
    AdhsFC = Adhs(~(adhEccIdx & adhLengthIdxFA));
    labelAdhesionFC = zeros(size(maskAdhesion2));
    for kk=1:numel(AdhsFC)
        labelAdhesionFC(AdhsFC(kk).PixelIdxList)=kk;
    end
    
    maskAdhesionFA = logical(labelAdhesionFA);
    maskAdhesionFC = maskAdhesion2 .* ~maskAdhesionFA;

    indInside=maskVectors(xNA,yNA,maskAdhesion2);
    indTrueNAs=~indInside;
    if ~isempty(pstruct)
        idxSigCCP = pstruct.A>0 & indTrueNAs';

        nascentAdhInfo(j).xCoord = [round(pstruct.x(idxSigCCP)'), round(pstruct.x_pstd(idxSigCCP)')];
        nascentAdhInfo(j).yCoord = [round(pstruct.y(idxSigCCP)'), round(pstruct.y_pstd(idxSigCCP)')];
        nascentAdhInfo(j).amp = [pstruct.A(idxSigCCP)', pstruct.A_pstd(idxSigCCP)'];
        if ~isempty(bandwidth)
            nascentAdhInfo(j).bandMask = bandMask;
        end
        nascentAdhInfo(j).numberNA = length(pstruct.x(idxSigCCP));
        if ~isempty(bandwidth)
            nascentAdhInfo(j).bandArea = sum(sum(bandMask))*(pixSize/1000)^2; % in um^2
        else
            nascentAdhInfo(j).bandArea = sum(ultimateMask(:))*(pixSize/1000)^2; % in um^2
        end
        nascentAdhInfo(j).NAdensity = nascentAdhInfo(j).numberNA/nascentAdhInfo(j).bandArea; % number per um2
    else
        nascentAdhInfo(j).xCoord = [];
        nascentAdhInfo(j).yCoord = [];
        nascentAdhInfo(j).amp = [];
        nascentAdhInfo(j).numberNA =0;
        nascentAdhInfo(j).NAdensity = NaN; % number per um2
    end
    % FA info
    % focal contact (FC) analysis
%         FCIdx = find(adhIdx);
    numAdhs = length(Adhs);
    boundFAs=bwboundaries(maskAdhesion2);
    for k=1:numAdhs
        focalAdhInfo(j).xCoord(k,1) = round(Adhs(k).Centroid(1));
        focalAdhInfo(j).yCoord(k,1) = round(Adhs(k).Centroid(2));
        focalAdhInfo(j).area(k,1) = Adhs(k).Area;
        focalAdhInfo(j).length(k,1) = Adhs(k).MajorAxisLength;
        focalAdhInfo(j).width(k,1) = Adhs(k).MinorAxisLength;
        focalAdhInfo(j).amp(k,1) = mean(I(Adhs(k).PixelIdxList));
        focalAdhInfo(j).ecc(k,1) = Adhs(k).Eccentricity;
    end
    focalAdhInfo(j).boundFA= boundFAs;
    focalAdhInfo(j).numberFA = numAdhs;
    focalAdhInfo(j).meanFAarea = mean(focalAdhInfo(j).area);
    focalAdhInfo(j).medianFAarea = median(focalAdhInfo(j).area);
    focalAdhInfo(j).meanLength = mean(focalAdhInfo(j).length);
    focalAdhInfo(j).medianLength = median(focalAdhInfo(j).length);
    focalAdhInfo(j).cellArea = sum(mask(:))*(pixSize/1000)^2; % in um^2
    focalAdhInfo(j).FAdensity = numAdhs/focalAdhInfo(j).cellArea; % number per um2
    focalAdhInfo(j).maskFA = maskAdhesion2;
    focalAdhInfo(j).labelFA = labelAdhesion;
    % FADensity at the cell periphery
    bandwidthNA = 5; %um
    bandwidthNA_pix = round(bandwidthNA*1000/movieData.pixelSize_);
    % Cell Boundary Mask 
    % mask for band from edge
    iMask = imcomplement(mask);
    distFromEdge = bwdist(iMask);
    bandMask = distFromEdge <= bandwidthNA_pix;

    maskOnlyBand = bandMask & mask;
    bandArea = sum(maskOnlyBand(:))*(pixSize/1000)^2; % in um^2

    % now see if these tracks ever in the maskOnlyBand
    indFAsAtEdge = false(numAdhs,1);
    for k=1:numAdhs
        if maskOnlyBand(round(focalAdhInfo(j).yCoord(k)),round(focalAdhInfo(j).xCoord(k)))
            indFAsAtEdge(k) = true;
        end
    end
    
    focalAdhInfo(j).FAdensityPeri = sum(indFAsAtEdge)/bandArea; % number per um2
    focalAdhInfo(j).FAdensityInside = (numAdhs-sum(indFAsAtEdge))/(focalAdhInfo(j).cellArea-bandArea); % number per um2
    
    % New feature: the other channel reading
    % 1. Check if there is the other channel
    % e.g. see if there is more than one channels
    if nChans >1 && (~isempty(iTFM) || ~isempty(iAdditionalChan))
        % 2. Go over each FA and FC and NA
        numFAs = numel(AdhsFA);
        numFCs = numAdhs - numFAs;
        numNAs = nascentAdhInfo(j).numberNA;
        if ~isempty(iTFM)
            curTmap = tMap(:,:,j);
            % FA first
            for ii=1:numFAs
                % 3. Get the specific segmentation-based mask
                curAdhMask = labelAdhesionFA==ii;
                % 4. Dilate the mask
                curAdhMaskDilated = bwmorph(curAdhMask,'dilate',1);
                % 5. Read from tMap
                forceFA{j}(ii) = mean(curTmap(curAdhMaskDilated));
            end
            % FC second
            for ii=1:numFCs
                % 3. Get the specific segmentation-based mask
                curAdhMask = labelAdhesionFC==ii;
                % 4. Dilate the mask
                curAdhMaskDilated = bwmorph(curAdhMask,'dilate',1);
                % 5. Read from tMap
                forceFC{j}(ii) = mean(curTmap(curAdhMaskDilated));
            end
            % NA next
            seNA = strel('disk',2);
            for ii=1:numNAs
                % 3. Get the specific segmentation-based mask
                curAdhMask = false(size(labelAdhesionFC));
                cutA=0; if nascentAdhInfo(j).yCoord(ii,1)-2 <=0, cutA= -nascentAdhInfo(j).yCoord(ii,1)+3; end
                cutC=0; if nascentAdhInfo(j).xCoord(ii,1)-2 <=0, cutC= -nascentAdhInfo(j).xCoord(ii,1)+3; end
                cutB=0; if nascentAdhInfo(j).yCoord(ii,1)+2 >size(curAdhMask,1), cutB= nascentAdhInfo(j).yCoord(ii,1)+2-size(curAdhMask,1); end
                cutD=0; if nascentAdhInfo(j).xCoord(ii,1)+2 >size(curAdhMask,2), cutD= nascentAdhInfo(j).xCoord(ii,1)+2-size(curAdhMask,2); end
                curAdhMask(nascentAdhInfo(j).yCoord(ii,1)-2+cutA:nascentAdhInfo(j).yCoord(ii,1)+2-cutB,...
                    nascentAdhInfo(j).xCoord(ii,1)-2+cutC:nascentAdhInfo(j).xCoord(ii,1)+2-cutD) = ...
                    seNA.Neighborhood(1+cutA:end-cutB,1+cutC:end-cutD);
                % 4. Dilate the mask
                curAdhMaskDilated = bwmorph(curAdhMask,'dilate',1);
                % 5. Read from tMap
                forceNA{j}(ii) = mean(curTmap(curAdhMaskDilated));
            end
            % Finally BG
            maskAdhesion2dilated = bwmorph(maskAdhesion2,'dilate',5);
            cellMask = roiMask(:,:,j) & mask;
            bgMask = roiMask(:,:,j) & ~mask;
            areaThres = round(mean(focalAdhInfo(j).area));
            curForceBGinCell = curTmap(~maskAdhesion2dilated & cellMask)';
            for ii=1:floor(length(curForceBGinCell)/areaThres)
                if ii<floor(length(curForceBGinCell)/areaThres)
                    forceBGinCell{j}(ii) = mean(curForceBGinCell((ii-1)*areaThres+1:(ii)*areaThres));
                else
                    forceBGinCell{j}(ii) = mean(curForceBGinCell((ii-1)*areaThres+1:end));
                end
            end
            
            curForceBGoutCell = curTmap(~maskAdhesion2dilated & bgMask)';
            for ii=1:floor(length(curForceBGoutCell)/areaThres)
                if ii<floor(length(curForceBGoutCell)/areaThres)
                    forceBGoutCell{j}(ii) = mean(curForceBGoutCell((ii-1)*areaThres+1:(ii)*areaThres));
                else
                    forceBGoutCell{j}(ii) = mean(curForceBGoutCell((ii-1)*areaThres+1:end));
                end
            end
            
        end
        % For the other channel
        if ~isempty(iAdditionalChan)
            % 6. Get the Ibg (mean intensity of the band in the channel image)
            chan2 = movieData.getChannel(iAdditionalChan);
            curImg2 = chan2.loadImage(j);
            maskAdhesionDilated = bwmorph(maskAdhesion,'dilate',5);
            cellMask = roiMask(:,:,j) & mask;
            M_bg = curTmap.*~maskAdhesionDilated .* cellMask;
            for ii=1:numAdhs
                % 3. Get the specific segmentation-based mask
                curAdhMask = labelAdhesion==ii;

                % 4. Dilate the mask
                curAdhMaskDilated = bwmorph(curAdhMask,'dilate',1);

                Ibg2d = M_bg .* double(curImg2);
                Ibg = mean(Ibg2d(:));

                % 7. Get the Iabs
                Iabs = mean(curImg2(curAdhMaskDilated));

                % 8. Get the ampTheOther
                ampTheOther = Iabs - Ibg;

                % 9. Save it to focalAdhInfo
                focalAdhInfo(j).ampTheOther(ii) = ampTheOther;
            end
        end
        
        % 10. colocalization analysis
        if ~isempty(iTFM)
            curTmap = tMap(:,:,j);
%             hScatter = figure; plot(I(:),curTmap(:),'.')
            cellMask2 = bwmorph(cellMask,'dilate',5) & roiMask(:,:,j);
            [~,xBins] = histcounts(I(cellMask2));  [~,yBins] = histcounts(curTmap(cellMask2));
            % Correlation - Pearson
            rP = corr(I(cellMask2), curTmap(cellMask2));
            % Mander's overlap coefficient, MOC
            products=(double(I(cellMask2)).*double(curTmap(cellMask2)));
            redsq=double(I(cellMask2)).^2;
            greensq=double(curTmap(cellMask2)).^2;

            MOC=sum(products(:))/sqrt(sum(redsq(:))*sum(greensq(:)));
            
            if plotGraphTFM
                hHist2D = figure; hold on
                densityplot(I(cellMask2), curTmap(cellMask2), xBins, yBins,'DisplayFunction', @log);
                colormap jet
                hCBar = colorbar;
                hCBar.Label.String = 'Counts (10^x )';
                ylabel('Traction (Pa)')
                xlabel('Adhesion Intensity (A.U.)')
                hold on
                text(xBins(round(0.09*length(xBins))),yBins(round(0.95*length(yBins))),....
                    ['Pearson''s corr: ' num2str(rP)], 'Color', 'w')
                text(xBins(round(0.09*length(xBins))),yBins(round(0.85*length(yBins))),...
                    ['Mander''s coeff: ' num2str(MOC)], 'Color', 'w')
                hHist2D.Units='inch';
                hHist2D.Position(3)=3; hHist2D.Position(4)=2.5;

                hgexport(hHist2D,[figPath filesep 'scatter2Dhist'],hgexport('factorystyle'),'Format','eps')
                hgsave(hHist2D,[figPath filesep 'scatter2Dhist'],'-v7.3')
                print(hHist2D,[figPath filesep 'scatter2Dhist'],'-dtiff')
                close(hHist2D)
            end

            % Save them
            tableCorr=table([rP; MOC],'RowNames',{'rP', 'MOC'});
            writetable(tableCorr,[dataPath filesep 'corrValues.csv'],'WriteRowNames',true)
            focalAdhInfo.rP(j) = rP;
            focalAdhInfo.MOC(j) = MOC;

            % 11. Now we are segmenting TFM map (force blob) and calculate
            % fraction of it inside the cell or outside.
            cellMask3 = bwmorph(cellMask,'dilate',100) & bwmorph(roiMask(:,:,j),'erode',15);
            maskForceBlob = blobSegmentThresholdTFM(curTmap,minSize,false,cellMask3);
            % fraction of the blob mask inside
            blobPixelsAll = sum(maskForceBlob(:));
            blobInside = sum(maskForceBlob(cellMask));
            blobOutside = sum(maskForceBlob(~cellMask));
            fractionBlobInside = blobInside/blobPixelsAll;
            fractionBlobOutside = blobOutside/blobPixelsAll;
            % Plot them
            focalAdhInfo.fractionBlobInside(j) = fractionBlobInside;
            focalAdhInfo.fractionBlobOutside(j) = fractionBlobOutside;
            if plotGraphTFM
                hBarFrac = figure; bar(categorical({'Inside','Outside'}), [fractionBlobInside fractionBlobOutside])
                hBarFrac.Units='inch';
                hBarFrac.Position(3)=3; hBarFrac.Position(4)=2.5;
                ylim([0 1]); title({'fraction of force blobs'; 'inside the cell or outside'})
                % Overlay them
                % cell mask
                combI(:,:,1) = 0.4*maskForceBlob.*cellMask + 0.6*cellMask;
                combI(:,:,2) = 0.7*maskForceBlob.*~cellMask + 0.6*cellMask - 0.5*maskForceBlob.*cellMask;
                combI(:,:,3) = -0.5*maskForceBlob.*cellMask + 0.6*cellMask;
                hBlob=figure; imshow(combI,[]); hold on
                hgexport(hBlob,[figPath filesep 'forceBlobOverlay'],hgexport('factorystyle'),'Format','eps')
                hgsave(hBlob,[figPath filesep 'forceBlobOverlay'],'-v7.3')
                print(hBlob,[figPath filesep 'forceBlobOverlay'],'-dtiff')
                close(hBarFrac)

                % Show force blobs on TFM image
                hTFMBlob=figure;
                imshow(curTmap, [0 tMax]), colormap jet
                hC=colorbar('east');
                hC.Color='w'; hC.Position(3:4)=[0.03 0.8];hC.Position(2)=0.1;
                hold on
                [B,~,nBD]  = bwboundaries(cellMask,'noholes');
                for kk=1:nBD
                    boundary = B{kk};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2) % cell boundary
                end
                [forceBlobIndiv,~,nFCs] = bwboundaries(maskForceBlob,4,'noholes');    
                for kk=1:nFCs
                    boundary = forceBlobIndiv{kk};
                    plot(boundary(:,2), boundary(:,1), 'Color','m', 'LineWidth', 2) % cell boundary
                end
                hgexport(hTFMBlob,[figPath filesep 'tractionMapWithForceBlobs'],hgexport('factorystyle'),'Format','eps')
                hgsave(hTFMBlob,[figPath filesep 'tractionMapWithForceBlobs'],'-v7.3')
                print(hTFMBlob,[figPath filesep 'tractionMapWithForceBlobs'],'-dtiff')
                close(hTFMBlob)
            end
        end        
    end

    % plotting detected adhesions
    if plotGraph
        h1=figure;
        maxI=max(I(:)); minI=min(I(:));
        dI = (double(I)-minI)/(maxI-minI);
        combI(:,:,1) = ~maskAdhesionFA.*dI+maskAdhesionFA.*(0.4*dI+double(maskAdhesionFA)*.5);
        combI(:,:,2) = ~maskAdhesionFA.*~maskAdhesionFC.*dI...
            +maskAdhesionFC.*(0.4*dI+double(maskAdhesionFC)*.6)...
            +maskAdhesionFA.*(0.4*dI);
        combI(:,:,3) = ~maskAdhesionFA.*dI+maskAdhesionFA.*(0.4*dI+double(maskAdhesionFA)*.0);%+double(ultimateMask)*.5;
        imshow(combI,[]); hold on
        if ~isempty(pstruct)
            plot(pstruct.x(idxSigCCP),pstruct.y(idxSigCCP),'yo')
        end
        hgexport(h1,strcat(tifPath,'/imgNAFA',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/imgNAFA',num2str(j,jformat)),'-v7.3')
        hold off
        
        if ~isempty(iTFM) && plotGraphTFM
            imshow(curTmap, [0 2000]), colormap jet
            hC=colorbar('east');
            hC.Color='w'; hC.Position(3:4)=[0.03 0.8];hC.Position(2)=0.1;
            hold on
            [B,~,nBD]  = bwboundaries(cellMask,'noholes');
            for kk=1:nBD
                boundary = B{kk};
                plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2) % cell boundary
            end
            [adhBoundFC,~,nFCs] = bwboundaries(maskAdhesionFC,4,'noholes');    
            for kk=1:nFCs
                boundary = adhBoundFC{kk};
                plot(boundary(:,2), boundary(:,1), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5) % cell boundary
            end
            [adhBoundFA,~,nFAs] = bwboundaries(maskAdhesionFA,4,'noholes');    
            for k = 1:nFAs
                adhBoundary = adhBoundFA{k};
                plot(adhBoundary(:,2), adhBoundary(:,1), 'm', 'LineWidth', 0.5) %adhesion boundary
            end
            if ~isempty(pstruct)
                plot(pstruct.x(idxSigCCP),pstruct.y(idxSigCCP),'yo')
            end
            hgexport(h1,[figPath filesep 'tractionMapWithAdhesions'],hgexport('factorystyle'),'Format','eps')
            hgsave(h1,[figPath filesep 'tractionMapWithAdhesions'],'-v7.3')
            print(h1,[figPath filesep 'tractionMapWithAdhesions'],'-dtiff')
        end
        close(h1)
    end
    if ~isempty(pstruct)
        disp(['detected ' num2str(numel(pstruct.x)) ' nascent adhesions, '  num2str(numFCs) ' focal complexes and ' num2str(numFAs) ' focal adhesions  for ' num2str(j) 'th frame.'])
    else
        disp(['detected zero nascent adhesion and ' num2str(numAdhs) ' focal adhesions for ' num2str(j) 'th frame.'])
    end
end

% Plot the force result
if ~isempty(iTFM)
    forceGroup = {forceFA, forceFC, forceNA, forceBGinCell, forceBGoutCell};
    nameList = {'FA', 'FC', 'NA', 'BG_inside', 'BG_outside'};
    forceGroupCell = cellfun(@(x) cell2mat(x),forceGroup,'unif',false);
    if plotGraphTFM
        h1=figure; 

        boxPlotCellArray(forceGroupCell,nameList,1,false,false);
        ylabel('Force magnitude (Pa)')
        title('Force per adhesion type')
        hgexport(h1,[figPath filesep 'forcePerAdhesionType'],hgexport('factorystyle'),'Format','eps')
        hgsave(h1,[figPath filesep 'forcePerAdhesionType'],'-v7.3')
        print(h1,[figPath filesep 'forcePerAdhesionType'],'-dtiff')
        close(h1)
    end
    % Save them
    tableForcePerAdhesionType=table(forceGroupCell','RowNames',nameList);
    writetable(tableForcePerAdhesionType,[dataPath filesep 'forcePerAdhesionType.csv'],'WriteRowNames',true)

    % Save them
    tableBlobFrac=table([focalAdhInfo.fractionBlobInside; focalAdhInfo.fractionBlobOutside],'RowNames',{'Inside','Outside'});
    writetable(tableBlobFrac,[dataPath filesep 'tableBlobFrac.csv'],'WriteRowNames',true)
end


save([dataPath filesep 'AdhInfo.mat'],'nascentAdhInfo','focalAdhInfo','-v7.3');

disp('Finished detecting objects...')

% detectMovieNascentAdhesion('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130410 paxilin crop/Cell12')