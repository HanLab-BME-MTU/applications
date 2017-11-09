function [] = readTheOtherChannelFromTracks(movieData,varargin)
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

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;
%Parse input, store in parameter structure
%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('TheOtherChannelReadingProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(TheOtherChannelReadingProcess(movieData));                                                                                                 
end
theOtherChanReadProc = movieData.processes_{iProc};
paramsIn = parseProcessParams(theOtherChanReadProc,paramsIn);

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(movieData,showAllTracks,plotEachTrack,varargin{:});
outputPath=paramsIn.outputPath;
iChanMaster=paramsIn.iChanMaster;
iChanSlave=paramsIn.iChanMaster;

%% Data Set up
% Set up the output file path for master channel
[~,iChanName]=fileparts(movieData.channels_(iChanMaster).channelPath_);
outputFilePath = [paramsIn.OutputDirectory filesep outputPath filesep iChanName];
mkClrDir(outputFilePath)
outputFile=strcat(outputFilePath,filesep,'tracksNA.mat');

%% Reading tracks from master channel
% Use the previous analysis folder structure
iAAProc = movieData.getProcessIndex('AdhesionAnalysisProcess');
adhAnalProc = movieData.getProcess(iAAProc);
paramsInAA=adhAnalProc.funParams_;
outputFilePathMaster = [paramsInAA.OutputDirectory filesep paramsInAA.outputPath filesep iChanName];
dataPathMaster = [outputFilePathMaster filesep 'data'];
outputFilePrev=strcat(dataPathMaster,filesep,'tracksNA.mat');
% See if you can use existing tracks
if exist(outputFilePrev,'file')
    disp('tracksNA file is found. Using it ... If you do not want to reuse this, please backup the file and rename it to other name than tracksNA.mat.')
    tracksNAFile = load(outputFilePrev,'tracksNA');
    tracksNA = tracksNAFile.tracksNA;
else
    disp('AdhesionAnalysisProcess must be run before running this.')
end
%% the other channel map stack - iChanSlave
% Build the interpolated TFM matrix first and then go through each track
% First build overall TFM images
% We don't need to build traction map again if this one is already built
% during force field calculation process.
[h,w]=size(movieData.channels_(iChanSlave).loadImage(1));
nFrames = movieData.nFrames_;
imgStack = zeros(h,w,nFrames);
for ii=1:nFrames
    imgStack(:,:,ii)=movieData.channels_(iChanSlave).loadImage(ii); 
end
%% Read force from imgStack
% get the intensity
disp('Reading the other channel...')
tic
tracksNA = readIntensityFromTracks(tracksNA,imgStack,2); % 2 means traction magnitude collection from traction stacks
toc

%% protrusion/retraction information - most of these are now done in analyzeAdhesionsMaturation
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.

disp('Post-analysis on adhesion movement and cross-correlation between fluorescence intensity and the other channel...')
deltaT = movieData.timeInterval_; % sampling rate (in seconds, every deltaT seconds)
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
    earlyPeriod = floor(1/timeInterval); % frames per minute
    lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
    lastFrameFromOne = lastFrame - sF+1;
    [curForceR,curForceM] = regression(tIntervalMin*(1:lastFrameFromOne),tracksNA(k).forceMag(sF:lastFrame));
    tracksNA(k).forceSlope = curForceM; % in Pa/min
    tracksNA(k).forceSlopeR = curForceR; % Pearson's correlation coefficient
end
disp('Saving...')
save(outputFile,'tracksNA','-v7.3');
disp('Done!')
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