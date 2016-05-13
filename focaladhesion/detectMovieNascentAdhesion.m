function [nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion(pathForTheMovieDataFile,bandwidth,iPax,plotGraph,indMask)
% detectMovieNascentAdhesion detect objects by fitting isotropic Gaussians
%
% SYNOPSIS adapted from detectMovieSubResFeatures(movieData,paramsIn)
%
% INPUT   
%   pathForTheMovieDataFile - path to a MovieData object describing the movie to be processed
%   bandwidth               - bandwidth in micron from the edge to qunatify
%                               nascent adhesions (default: 5 um)
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
end
if nargin<3
    iPax = 1; %assumed
    plotGraph = true;
    indMask = [];
end
if nargin<4
    plotGraph = true;
    indMask = [];
end
if nargin<5
    indMask = []; %indMask is the channel index of a mask to be added to mask in iPax channel.
end
% Load the MovieData
if isa(pathForTheMovieDataFile, 'MovieData')
    movieData = pathForTheMovieDataFile;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    movieData = MovieData.load(movieDataPath,false);
end

%% --------------- Initialization ---------------%%
nChan=numel(movieData.channels_);
% psfSigma = 1.2; %hard coded for nascent adhesion
psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(iPax).emissionWavelength_*1e-9);

% Set up the output directories
outputFilePath = [movieData.outputDirectory_ filesep 'Adhesion Quantification'];
imgPath = [outputFilePath filesep 'imgs'];
dataPath = [outputFilePath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];

if ~exist(figPath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(tifPath);
    mkdir(figPath);
end

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting isotropic Gaussians...')

roiMask = movieData.getROIMask;
nascentAdhInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'NAdensity',[],'bandMask',[],'numberNA',[],'bandArea',[]);
focalAdhInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'length',[],'meanFAarea',[],'medianFAarea',[]...
    ,'meanLength',[],'medianLength',[],'numberFA',[],'FAdensity',[],'cellArea',[],'maskFA',[],'boundFA',[],'ecc',[]);
if plotGraph
    h1=figure;
end
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
pixSize = movieData.pixelSize_;
minSize = round((500/pixSize)*(500/pixSize)); %adhesion limit=0.25 um2
minEcc = 0.7;

disp(['Paxillin channel was assumed to be in channel ' num2str(iPax) '.'])
disp('Results will be saved under:')
disp(outputFilePath);

for j=1:movieData.nFrames_
    I=double(movieData.channels_(iPax).loadImage(j));
    noMask=false;
    try
        maskProc = movieData.getProcess(movieData.getProcessIndex('MaskRefinementProcess'));
%         mask = maskProc.loadChannelOutput(iPax,j);
        mask = maskProc.loadChannelOutput(iPax,j); % 1 is CCP channel
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
    maskAdhesion = blobSegmentThreshold(I,minSize,0,mask);
    maskAdhesionFine = blobSegmentThreshold(I,0); % mask for all adhesions without minSize

    labelAdhesion = bwlabel(maskAdhesion);
    Adhs = regionprops(maskAdhesion,'Centroid','Area','Eccentricity','PixelIdxList','MajorAxisLength');
%         minFASize = round((2000/MD.pixelSize_)*(500/MD.pixelSize_)); %adhesion limit=1um*.5um

    adhIdx = arrayfun(@(x) x.Eccentricity>minEcc, Adhs);
    Adhs = Adhs(adhIdx);
    labelAdhesion = false(size(maskAdhesion));
    for kk=1:numel(Adhs)
        labelAdhesion(Adhs(kk).PixelIdxList)=true;
    end
    maskAdhesion = logical(labelAdhesion);
    maskAdhesionC = imcomplement(maskAdhesion);
    % mask for band from edge
    if noMask
        ultimateMask = roiMask(:,:,j) & maskAdhesionC & mask & maskAdhesionFine;
    elseif ~isempty(bandwidth)
        iMask = imcomplement(mask);
        distFromEdge = bwdist(iMask);
        bandwidth_pix = round(bandwidth*1000/pixSize);
        bandMask = distFromEdge <= bandwidth_pix;

        ultimateMask = bandMask & roiMask(:,:,j) & maskAdhesionC & mask & maskAdhesionFine;
    else
        ultimateMask = roiMask(:,:,j) & maskAdhesionC & mask; % & maskAdhesionFine;
    end
    pstruct = pointSourceDetection(I, psfSigma,  ...
        'Alpha',0.05,'Mask',ultimateMask);
        % filter out points where overall intensity is in the noise level
    pixelIntenMargin = I(~mask);
    maxIntBg=quantile(pixelIntenMargin,0.99);
    psInt = pstruct.A+pstruct.c;
    idxSigCCP = psInt>maxIntBg;

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

    % FA info
    % focal contact (FC) analysis
%         FCIdx = find(adhIdx);
    numAdhs = length(Adhs);
    boundFAs=bwboundaries(maskAdhesion);
    for k=1:numAdhs
        focalAdhInfo(j).xCoord(k) = round(Adhs(k).Centroid(1));
        focalAdhInfo(j).yCoord(k) = round(Adhs(k).Centroid(2));
        focalAdhInfo(j).area(k) = Adhs(k).Area;
        focalAdhInfo(j).length(k) = Adhs(k).MajorAxisLength;
        focalAdhInfo(j).amp(k) = mean(I(Adhs(k).PixelIdxList));
        focalAdhInfo(j).ecc(k) = Adhs(k).Eccentricity;
    end
    focalAdhInfo(j).boundFA= boundFAs;
    focalAdhInfo(j).numberFA = numAdhs;
    focalAdhInfo(j).meanFAarea = mean(focalAdhInfo(j).area);
    focalAdhInfo(j).medianFAarea = median(focalAdhInfo(j).area);
    focalAdhInfo(j).meanLength = mean(focalAdhInfo(j).length);
    focalAdhInfo(j).medianLength = median(focalAdhInfo(j).length);
    focalAdhInfo(j).cellArea = sum(mask(:))*(pixSize/1000)^2; % in um^2
    focalAdhInfo(j).FAdensity = numAdhs/focalAdhInfo(j).cellArea; % number per um2
    focalAdhInfo(j).maskFA = maskAdhesion;

    % plotting detected adhesions
    if plotGraph
        dI = double(I)/max(max(I));
        combI(:,:,1) = dI;
        combI(:,:,2) = dI+double(maskAdhesion)*.5;
        combI(:,:,3) = dI;%+double(ultimateMask)*.5;
        imshow(combI,[]); hold on
        plot(pstruct.x,pstruct.y,'ro')
        hgexport(h1,strcat(tifPath,'/imgNAFA',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/imgNAFA',num2str(j,jformat)),'-v7.3')
        hold off
        close(h1)
    end
    disp(['detected ' num2str(length(pstruct.x)) ' nascent adhesions and ' num2str(numAdhs) ' focal adhesions for ' num2str(j) 'th frame.'])
end
save([dataPath filesep 'AdhInfo.mat'],'nascentAdhInfo','focalAdhInfo','-v7.3');

disp('Finished detecting objects...')

% detectMovieNascentAdhesion('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130410 paxilin crop/Cell12')