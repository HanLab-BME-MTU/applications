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
p.OutputDirectory=[MD.getPath filesep 'FocalAdhesionPackage' filesep 'TheOtherChannelReading'];
theOtherChanReadProc.setPara(p)

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(MD,showAllTracks,plotEachTrack,varargin{:});
% iChanMaster=p.ChannelIndex;
iChanSlave=p.iChanSlave;

%% Data Set up
% Set up the output file path for master channel
% outputFile = cell(2, numel(MD.channels_));
% for i = [p.ChannelIndex iChanSlave]
%     [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
%     outFilename = [chanDirName '_Chan' num2str(i) '_tracksNA'];
%     outputFile{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_ampTotal2PerNAFCFA'];
%     outputFile{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
% end
% theOtherChanReadProc.setOutFilePaths(outputFile);
nChan = cellfun(@numel,p.ChannelIndex);
nInput = sum(nChan);
imSize = MD.imSize_;

mkClrDirWithBackup(p.OutputDirectory);
nFrames = MD.nFrames_;

%Set up and store the output directories for the window samples.
outFilePaths =cell(numel(p.ProcessIndex),numel(MD.channels_),2);
for i=1:numel(p.ProcessIndex)
    iProc = p.ProcessIndex{i};
    if isempty(iProc)
        pString='Raw images - channel ';
    else
        parentOutput = MD.processes_{iProc}.getDrawableOutput;
        iOutput = strcmp(p.OutputName{i},{parentOutput.var});
        pString=[parentOutput(iOutput).name ' - channel '];
    end
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        outFilePaths{i,iChan,1} =  [p.OutputDirectory filesep pString num2str(iChan) '.mat'];
        outFilePaths{i,iChan,2} =  [p.OutputDirectory filesep pString num2str(iChan) '-intensities.mat'];
    end
end
theOtherChanReadProc.setOutFilePaths(outFilePaths);

%% Loading tracks from master channel
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
pAdh = adhAnalProc.funParams_;
tracksNA=adhAnalProc.loadChannelOutput(pAdh.ChannelIndex,'output','tracksNA');
%% Code matching
% Depending on the input we need to separate the codes for the image stack
% to read from. 
iPureChan = 0; %counting pure other channel(s) 
codeSet = []; %for storing codes
chanSet = []; %for sotring corresponding channel id
outputNameSet = {}; %for output names
procSet = [];

for i=1:numel(p.ProcessIndex)
    iProc = p.ProcessIndex{i};
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        iInput = sum(nChan(1:i-1))+j;

        if ~isempty(iProc)
            procName = MD.processes_{iProc}.getName;
            if strcmpi(procName,'Ratioing')
                iReadingCode = 3; %fret
                codeSet = [codeSet iReadingCode];
                chanSet = [chanSet iChan];
                outputNameSet = [outputNameSet p.OutputName(i)];
                procSet = [procSet i];
            else
                error('Code for TFM or qFSM reading is not established yet. Unselect them from the setting.')
            end
            % stack2sample(:,:,iInput) = MD.processes_{iProc}.loadChannelOutput(iChan,iFrame,'output',p.OutputName{i});
        else
            iPureChan=iPureChan+1;
            if iPureChan==1
                iReadingCode=5;
                codeSet = [codeSet iReadingCode];
                chanSet = [chanSet iChan];
                outputNameSet = [outputNameSet p.OutputName(i)];
                procSet = [procSet i];
            elseif iPureChan==2
                iReadingCode=6;
                codeSet = [codeSet iReadingCode];
                chanSet = [chanSet iChan];
                outputNameSet = [outputNameSet p.OutputName(i)];
                procSet = [procSet i];
            else
                error('the current code doesn not support more than 2 channels to read. Please choose less channels in the setting')
            end           
            % stack2sample(:,:,iInput) = MD.channels_(iChan).loadImage(iFrame);
        end
    end
end 
%% Slave to Master registration using PIV
% p.doMultimodalRegistration = true;
% p.doFAregistration = false;
if p.doFAregistration && iPureChan==1
    % Get the images
    iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
    FAPackage=MD.packages_{iFAPack}; iSDCProc=1;
    SDCProc=FAPackage.processes_{iSDCProc};
    if ~isempty(SDCProc)
        mainI = SDCProc.loadChannelOutput(p.ChannelIndex,1);
        subI = SDCProc.loadChannelOutput(p.iChanSlave,1);
    else
        mainI = MD.getChannel(p.ChannelIndex).loadImage(1);
        subI = MD.getChannel(p.iChanSlave).loadImage(1);
    end
    figure,imshowpair(mainI,subI,'montage')
    title('Unregistered')
    imshowpair(mainI,subI)
    title('Unregistered')
    
    % Improving registration 
    [optimizer,metric] = imregconfig('multimodal');
    subRegisteredDefault = imregister(subI,mainI,'similarity',optimizer,metric);

    figure,imshowpair(subRegisteredDefault,mainI)
    title('A: Default Registration')
%     disp(optimizer)
%     disp(metric)
    optimizer.InitialRadius = optimizer.InitialRadius/3.5;
    subRegisteredAdjustedInitialRadius = imregister(subI,mainI,'similarity',optimizer,metric);
%     figure,imshowpair(subRegisteredAdjustedInitialRadius,mainI)
%     title('B: Adjusted InitialRadius')

    optimizer.MaximumIterations = 300;
    subRegisteredAdjustedInitialRadius300 = imregister(subI,mainI,'similarity',optimizer,metric);
%     figure,imshowpair(subRegisteredAdjustedInitialRadius300,mainI)
%     title('C: Adjusted InitialRadius, MaximumIterations = 300')    
    
    % Step 4: Use Initial Conditions to Improve Registration
    tformSimilarity = imregtform(subI,mainI,'similarity',optimizer,metric);
    Rfixed = imref2d(size(mainI));
    movingRegisteredRigid = imwarp(subI,tformSimilarity,'OutputView',Rfixed);
%     figure, imshowpair(movingRegisteredRigid, mainI)
%     title('D: Registration Based on Similarity Transformation Model')
    
    subRegisteredAffineWithIC = imregister(subI,mainI,'similarity',optimizer,metric,...
    'InitialTransformation',tformSimilarity);
%     figure, imshowpair(subRegisteredAffineWithIC,mainI)
%     title('E: Registration from Affine Model Based on Similarity Initial Condition')
    % Apply to the slave channel
    ref_obj = imref2d(size(subI));
    imageFileNames = MD.getImageFileNames;
%     newSubI = imwarp(subI, invCurXform, 'OutputView', ref_obj);
    
    % Overwrite the side channel
    if ~isempty(SDCProc)
        % back up
        backupFolder = [SDCProc.outFilePaths_{1,p.iChanSlave} ' Original']; 
        if ~exist(backupFolder,'dir')
            mkdir(backupFolder);
        end
        % Anonymous functions for reading input/output
        outFile=@(chan,frame) [SDCProc.outFilePaths_{1,chan} filesep imageFileNames{chan}{frame}];
        ii=1;
        mainI = SDCProc.loadChannelOutput(p.ChannelIndex,ii);
        subI = SDCProc.loadChannelOutput(p.iChanSlave,ii);
        newSubI = imregister(subI,mainI,'similarity',optimizer,metric,...
            'InitialTransformation',tformSimilarity);
        copyfile(outFile(p.iChanSlave, ii),backupFolder,'f')
        imwrite(uint16(newSubI), outFile(p.iChanSlave, ii));
            
        hFig = figure;
        hAx  = subplot(1,2,1);
        C = imfuse(mainI,subI,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        imshow(C,[],'Parent', hAx);

        title('Original result: red, ch 1, green, ch 2')
        hAx2  = subplot(1,2,2);
        C2 = imfuse(mainI,newSubI,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        imshow(C2,[],'Parent', hAx2);
        title('Aligned result: red, ch 1, green, adjusted ch 2')

        iMainChan=pAdh.ChannelIndex; iSlaveChan = iPureChan;
        parfor ii=2:nFrames
            mainI = SDCProc.loadChannelOutput(iMainChan,ii);
            subI = SDCProc.loadChannelOutput(iSlaveChan,ii);
            newSubI = imregister(subI,mainI,'similarity',optimizer,metric,...
                'InitialTransformation',tformSimilarity);
            copyfile(outFile(p.iChanSlave, ii),backupFolder,'f')
            imwrite(uint16(newSubI), outFile(p.iChanSlave, ii));
        end
        disp(['Channel has been overwritten in ' SDCProc.outFilePaths_{1,p.iChanSlave}])
    else
        backupFolder = [MD.getChannel(p.iChanSlave).getPath ' Original']; 
        if ~exist(backupFolder,'dir')
            mkdir(backupFolder);
        end
        % Anonymous functions for reading input/output
        outFile=@(chan,frame) [MD.getChannel(chan).getPath filesep imageFileNames{chan}{frame}];
        for ii=1:nFrames
            subI = MD.getChannel(p.iChanSlave).loadImage(ii);
            mainI = MD.getChannel(p.ChannelIndex).loadImage(ii);
            newSubI = imregister(subI,mainI,'similarity',optimizer,metric,...
                'InitialTransformation',tformSimilarity);
            copyfile(outFile(p.iChanSlave, ii),backupFolder)
            imwrite(uint16(newSubI), outFile(p.iChanSlave, ii));
        end
        disp(['Channel has been overwritten in ' MD.getChannel(chan).getPath])
    end
else
    iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
    FAPackage=MD.packages_{iFAPack}; iSDCProc=1;
    SDCProc=FAPackage.processes_{iSDCProc};
end

%% Image Stack building
%iReadingCode is used in readIntensity function
% 1: ampTotal, 
% 2: forceMag, 
% 3: fret, 
% 4: flowSpeed, 
% 5: ampTotal2 and 
% 6: ampTotal3
readingCodeNames ={'ampTotal','forceMag','fret','flowSpeed','ampTotal2','ampTotal3'};
kk=0;
for iReadingCode = codeSet
    kk=kk+1;
    imgStack = zeros(imSize(1),imSize(2),nFrames);
    if iReadingCode > 4 % if it is about raw image
        if ~isempty(SDCProc)
            for ii=1:nFrames
                imgStack(:,:,ii)=SDCProc.loadChannelOutput(chanSet(kk),ii);
            end
        else
            for ii=1:nFrames
                imgStack(:,:,ii)=MD.channels_(chanSet(kk)).loadImage(ii); 
            end
        end
    else
        for ii=1:nFrames
            imgStack(:,:,ii) = MD.processes_{iProc}.loadChannelOutput(chanSet(kk),ii,'output',outputNameSet{i});
        end        
    end
    %% Read signal from imgStack
    % get the intensity
    disp('Reading the other channel...')
    tic
    if iReadingCode==5
        fieldsAll = fieldnames(tracksNA);
        unnecFields = setdiff(fieldsAll,{'xCoord','yCoord','startingFrameExtraExtra',...
            'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
            'endingFrame','amp','sigma','state','FAPixelList'});
        essentialTracks = rmfield(tracksNA,unnecFields);
        % Trying background subtraction here
        disp('Performing smart background subtraction ...')
        tic
        imgClass = class(imgStack);
        imgStackBS = zeros(size(imgStack));
        parfor ii=1:nFrames
            curImg = imgStack(:,:,ii);
            % Getting adhesion mask
            blobMask = blobSegmentThresholdGeneral(curImg,'rosin',1,1,3); %1,'I');
            % Make adhesion-deleted image
            blobMask2 = bwmorph(blobMask,'dilate',3);
            %imageBlobSubtracted = blobMask .* curImg;
            % Get X Y for blobMask
            [rows,columns] = ind2sub(size(curImg),find(~blobMask2));
            fnInterp = scatteredInterpolant(columns,rows,curImg(~blobMask2),'nearest','nearest');
            [X, Y] = meshgrid(1:size(curImg,2),1:size(curImg,1));%ind2sub(size(curImg),find(curImg>0));
            % Interpolate images
            imgFilled = fnInterp(X,Y);
            % blur this img
            imageBackground = filterGauss2D(imgFilled,50);
            %calculate noise-filtered and background-subtracted image
            imgStackBS(:,:,ii) = curImg - cast(imageBackground,imgClass);
        end
        toc
        
        addedTracksNA = readIntensityFromTracks(essentialTracks,imgStack,iReadingCode,'imgStackBS',imgStackBS); % 5 means ampTotal2 from the other channel
        intensitiesInNAs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==2)),tracksNA,addedTracksNA);
        intensitiesInFCs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==3)),tracksNA,addedTracksNA);
        intensitiesInFAs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==4)),tracksNA,addedTracksNA);
        intensityGroup={intensitiesInNAs, intensitiesInFCs, intensitiesInFAs};
    
        amplitudeInNAs = arrayfun(@(x,y) nanmean(y.amp2(x.state==2)),tracksNA,addedTracksNA);
        amplitudeInFCs = arrayfun(@(x,y) nanmean(y.amp2(x.state==3)),tracksNA,addedTracksNA);
        amplitudeInFAs = arrayfun(@(x,y) nanmean(y.amp2(x.state==4)),tracksNA,addedTracksNA);
        amplitudeGroup={amplitudeInNAs, amplitudeInFCs, amplitudeInFAs};
        save(outFilePaths{procSet(kk),chanSet(kk), 2},'intensityGroup','amplitudeGroup');
    else
        addedTracksNA = readIntensityFromTracks(essentialTracks,imgStack,iReadingCode); % 6 means ampTotal3 from the other channel amplitudeInNAs = arrayfun(@(x,y) nanmean(y.amp2(x.state==2)),tracksNA,addedTracksNA);
        intensitiesInNAs = arrayfun(@(x) nanmean(x.(readingCodeNames{iReadingCode})(x.state==2)),addedTracksNA);
        intensitiesInFCs = arrayfun(@(x) nanmean(x.(readingCodeNames{iReadingCode})(x.state==3)),addedTracksNA);
        intensitiesInFAs = arrayfun(@(x) nanmean(x.(readingCodeNames{iReadingCode})(x.state==4)),addedTracksNA);
        intensityGroup={intensitiesInNAs, intensitiesInFCs, intensitiesInFAs};
        save(outFilePaths{procSet(kk),chanSet(kk), 2},'intensityGroup');
    end
    tracksAmpTotal = rmfield(addedTracksNA,{'xCoord','yCoord','startingFrameExtraExtra',...
    'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
    'endingFrame','amp','sigma'});
    toc
    %% Add report of mean intensity of the other channel per status
    disp('Saving...')
    try
        save(outFilePaths{procSet(kk),chanSet(kk), 1},'tracksAmpTotal');
    catch
        save(outFilePaths{procSet(kk),chanSet(kk), 1},'tracksAmpTotal','-v7.3');
    end
end
%% protrusion/retraction information - most of these are now done in analyzeAdhesionsMaturation
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.


disp('Done!')
end
%
