function [] = readTheOtherChannelFromTracks(MD,varargin)
% colocalizationTFMwithFA runs colocalization between peaks in TFM maps and peaks in paxillin using MD.
% Basically this function make grayscale of TFM and pick up peaks of the
% map and see if how many of them are colocalized with significant paxillin
% signal or neighborhood of nascent adhesions.

% input:    pathForTheMovieDataFile:    path to the MD file (TFM
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
ip.addRequired('MD', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(MD,varargin{:});
paramsIn=ip.Results.paramsIn;
%Parse input, store in parameter structure
%Get the indices of any previous threshold processes from this function                                                                              
iProc = MD.getProcessIndex('TheOtherChannelReadingProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(MD.processes_)+1;
    MD.addProcess(TheOtherChannelReadingProcess(MD));                                                                                                 
end
theOtherChanReadProc = MD.processes_{iProc};
p = parseProcessParams(theOtherChanReadProc,paramsIn);

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(MD,showAllTracks,plotEachTrack,varargin{:});
% iChanMaster=p.ChannelIndex;
iChanSlave=p.iChanSlave;

%% Reading tracks from master channel
% Use the previous analysis folder structure
% try
%     iAdhProc = MD.getProcessIndex('AdhesionClassificationProcess');
%     adhAnalProc = MD.getProcess(iAdhProc);
%     % numChans = numel(p.ChannelIndex);
%     tracksNA = load(adhAnalProc.outFilePaths_{5,p.ChannelIndex},'tracksNA');
%     tracksNA = tracksNA.tracksNA;
% catch
% We read tracksNA from AdhesionAnalysisProcess because one in
% Classification step has the identical one and now the step doesn't store
% tracksNA.
iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
adhAnalProc = MD.getProcess(iAdhProc);
% numChans = numel(p.ChannelIndex);
s = load(adhAnalProc.outFilePaths_{1,p.ChannelIndex},'metaTrackData');
metaTrackData = s.metaTrackData;
fString = ['%0' num2str(floor(log10(metaTrackData.numTracks))+1) '.f'];
numStr = @(trackNum) num2str(trackNum,fString);
trackIndPath = @(trackNum) [metaTrackData.trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
for ii=metaTrackData.numTracks:-1:1
    curTrackObj = load(trackIndPath(ii),'curTrack');
    tracksNA(ii) = curTrackObj.curTrack;
end
% end    
%% the other channel map stack - iChanSlave
% Build the interpolated TFM matrix first and then go through each track
% First build overall TFM images
% We don't need to build traction map again if this one is already built
% during force field calculation process.

% Go through each channel and if one is related to TFM, we skip it. It will
% be dealt in the next process, 'ColocalizationWithTFMProcess'.
% Find out which channel was used for TFMPackage
% (Since it is 1 that is used for beads, I'll just check if there is
% TFMPackage run).
iTFMPackage = MD.getPackageIndex('TFMPackage');
if ~isempty(iTFMPackage)
%     curTfmPackage = MD.getPackage(iTFMPackage);
    iForceProc = MD.getProcessIndex('ForceFieldCalculationProcess');
    if ~isempty(iForceProc)
        forceProc = MD.getProcess(iForceProc);
        if forceProc.success_
            iChanSlave=setdiff(iChanSlave,1);
        end
    end
end
%% Data Set up
% Set up the output file path for master channel
outputFile = cell(1, numel(MD.channels_));
for i = p.ChannelIndex
    [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_Chan' num2str(i) '_tracksNA'];
    outputFile{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
end
theOtherChanReadProc.setOutFilePaths(outputFile);
mkClrDir(p.OutputDirectory);
%% Reading
iReadingCode=4;
for iCurChan = iChanSlave
    iReadingCode = iReadingCode+1;
    [h,w]=size(MD.channels_(iCurChan).loadImage(1));
    nFrames = MD.nFrames_;
    imgStack = zeros(h,w,nFrames);
    for ii=1:nFrames
        imgStack(:,:,ii)=MD.channels_(iCurChan).loadImage(ii); 
    end
    %% Read force from imgStack
    % get the intensity
    disp('Reading the other channel...')
    tic
    if iReadingCode==5
        fieldsAll = fieldnames(tracksNA);
        unnecFields = setdiff(fieldsAll,{'xCoord','yCoord','startingFrameExtraExtra',...
            'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
            'endingFrame','amp'});
        essentialTracks = rmfield(tracksNA,unnecFields);
        addedTracksNA = readIntensityFromTracks(essentialTracks,imgStack,iReadingCode); % 5 means ampTotal2 from the other channel
    else
        addedTracksNA = readIntensityFromTracks(addedTracksNA,imgStack,iReadingCode); % 5 means ampTotal2 from the other channel
    end
    tracksAmpTotal = rmfield(addedTracksNA,{'xCoord','yCoord','startingFrameExtraExtra',...
    'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
    'endingFrame','amp'});
    toc
end
%% protrusion/retraction information - most of these are now done in analyzeAdhesionsMaturation
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.

disp('Saving...')
save(outputFile{1,p.ChannelIndex},'tracksAmpTotal','-v7.3'); % the later channel has the most information.
disp('Done!')
end
%
