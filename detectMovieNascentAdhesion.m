function [nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion(pathForTheMovieDataFile,bandwidth)
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
end

% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
movieData = MovieData.load(movieDataPath,false);

%% --------------- Initialization ---------------%%
nChan=numel(movieData.channels_);
psfSigma = 2; %hard coded for nascent adhesion
% psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(1).emissionWavelength_*1e-9);

% Set up the output directories
outputFilePath = [pathForTheMovieDataFile filesep 'Adhesion Quantification'];
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
    'amp',[],'area',[],'NAperArea',[],'bandMask',[],'numberNA',[],'bandArea',[]);
focalAdhInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[],'area',[],'length',[],'meanFAarea',[],'medianFAarea',[]...
    ,'meanLength',[],'medianLength',[],'numberFA',[]);
h1=figure;
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
minSize = round((600/movieData.pixelSize_)*(400/movieData.pixelSize_)); %adhesion limit=1um*.5um

iPax = 1; %assumed
    disp(['Paxillin channel was assumed to be in channel ' num2str(iPax) '.'])
    disp('Results will be saved under:')
    disp(outputFilePath);
    
    for j=1:movieData.nFrames_
        I=double(movieData.channels_(iPax).loadImage(j));
        maskProc = movieData.getProcess(movieData.getProcessIndex('MaskRefinementProcess'));
        mask = maskProc.loadChannelOutput(1,j);
%         maskProc = movieData.processes_{2};
%         mask = maskProc.loadChannelOutput(1,j);
        maskAdhesion = blobSegmentThreshold(I,minSize,0,mask);
        maskAdhesionFine = blobSegmentThreshold(I,0); % mask for all adhesions without minSize
        maskAdhesionC = imcomplement(maskAdhesion);
        % mask for band from edge
        iMask = imcomplement(mask);
        distFromEdge = bwdist(iMask);
        bandwidth_pix = round(bandwidth*1000/movieData.pixelSize_);
        bandMask = distFromEdge <= bandwidth_pix;
        
        ultimateMask = bandMask & roiMask(:,:,j) & maskAdhesionC & mask & maskAdhesionFine;
        pstruct = pointSourceDetection(I, psfSigma,  ...
            'Alpha',0.05,'Mask',ultimateMask);
        nascentAdhInfo(j).xCoord = [round(pstruct.x'), round(pstruct.x_pstd')];
        nascentAdhInfo(j).yCoord = [round(pstruct.y'), round(pstruct.y_pstd')];
        nascentAdhInfo(j).amp = [pstruct.A', pstruct.A_pstd'];
        nascentAdhInfo(j).bandMask = bandMask;
        nascentAdhInfo(j).numberNA = length(pstruct.x);
        nascentAdhInfo(j).bandArea = sum(sum(bandMask))*(movieData.pixelSize_/1000)^2; % in um^2
        nascentAdhInfo(j).NAperArea = nascentAdhInfo(j).numberNA/nascentAdhInfo(j).bandArea;
        
        % FA info
        % focal contact (FC) analysis
        Adhs = regionprops(maskAdhesion,'Centroid','Area','Eccentricity','PixelIdxList','MajorAxisLength');
%         minFASize = round((2000/MD.pixelSize_)*(500/MD.pixelSize_)); %adhesion limit=1um*.5um

%         adhIdx = arrayfun(@(x) x.Area<minFASize & x.Eccentricity<0.95, Adhs);
%         FCs = Adhs(adhIdx);
%         FCIdx = find(adhIdx);
        numAdhs = length(Adhs);
        for k=1:numAdhs
            focalAdhInfo(j).xCoord(k) = round(Adhs(k).Centroid(1));
            focalAdhInfo(j).yCoord(k) = round(Adhs(k).Centroid(2));
            focalAdhInfo(j).area(k) = Adhs(k).Area;
            focalAdhInfo(j).length(k) = Adhs(k).MajorAxisLength;
            focalAdhInfo(j).amp(k) = mean(I(Adhs(k).PixelIdxList));
        end
        focalAdhInfo(j).numberFA = numAdhs;
        focalAdhInfo(j).meanFAarea = mean(focalAdhInfo(j).area);
        focalAdhInfo(j).medianFAarea = median(focalAdhInfo(j).area);
        focalAdhInfo(j).meanLength = mean(focalAdhInfo(j).length);
        focalAdhInfo(j).medianLength = median(focalAdhInfo(j).length);

        % plotting detected adhesions
        dI = double(I)/max(max(I));
        combI(:,:,1) = dI;
        combI(:,:,2) = dI+double(maskAdhesion)*.5;
        combI(:,:,3) = dI+double(ultimateMask)*.5;
        imshow(combI,[]); hold on
        plot(pstruct.x,pstruct.y,'ro')
        hgexport(h1,strcat(tifPath,'/imgNAFA',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/imgNAFA',num2str(j,jformat)),'-v7.3')
        hold off
        disp(['detected ' num2str(length(pstruct.x)) ' nascent adhesions and ' num2str(numAdhs) 'focal adhesions for ' num2str(j) 'th frame.'])
    end
    save([dataPath filesep 'AdhInfo.mat'],'nascentAdhInfo','focalAdhInfo','-v7.3');

disp('Finished detecting objects...')

% detectMovieNascentAdhesion('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130410 paxilin crop/Cell12')