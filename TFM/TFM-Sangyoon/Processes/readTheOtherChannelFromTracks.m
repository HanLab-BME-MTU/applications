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
tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA');

% numChans = numel(p.ChannelIndex);
% s = load(adhAnalProc.outFilePaths_{1,p.ChannelIndex},'metaTrackData');
% metaTrackData = s.metaTrackData;
% fString = ['%0' num2str(floor(log10(metaTrackData.numTracks))+1) '.f'];
% numStr = @(trackNum) num2str(trackNum,fString);
% trackIndPath = @(trackNum) [metaTrackData.trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
% for ii=metaTrackData.numTracks:-1:1
%     curTrackObj = load(trackIndPath(ii),'curTrack');
%     tracksNA(ii,1) = curTrackObj.curTrack;
%     progressText((metaTrackData.numTracks-ii)/metaTrackData.numTracks,'Loading tracksNA') % Update text
% end
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
outputFile = cell(2, numel(MD.channels_));
for i = [p.ChannelIndex iChanSlave]
    [~, chanDirName, ~] = fileparts(MD.getChannelPaths{i});
    outFilename = [chanDirName '_Chan' num2str(i) '_tracksNA'];
    outputFile{1,i} = [p.OutputDirectory filesep outFilename '.mat'];

    outFilename = [chanDirName '_Chan' num2str(i) '_ampTotal2PerNAFCFA'];
    outputFile{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
end
theOtherChanReadProc.setOutFilePaths(outputFile);
mkClrDirWithBackup(p.OutputDirectory);
nFrames = MD.nFrames_;
%% Slave to Master registration using PIV
% p.doMultimodalRegistration = true;
% p.doFAregistration = false;
if p.doFAregistration
%     % Find the FA locations in the master channel at t=1
%     iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
%     FAPackage=MD.packages_{iFAPack}; iSDCProc=1;
%     SDCProc=FAPackage.processes_{iSDCProc};
% 
%     iFAsegProc = MD.getProcessIndex('FocalAdhesionSegmentationProcess');
%     FASegProc = MD.getProcess(iFAsegProc);
%     maskFAs = FASegProc.loadChannelOutput(p.ChannelIndex,1);
%     if ~isempty(SDCProc)
%         iBeadChan = 1; % might need to be updated based on asking TFMPackage..
%         s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
%         T = s.T;
%         
%         maskFAs = imtranslate(maskFAs,T(1,2:-1:1));
%         mainI = SDCProc.loadChannelOutput(p.ChannelIndex,1);
%         subI = SDCProc.loadChannelOutput(p.iChanSlave,1);
%     else
%         mainI = MD.getChannel(p.ChannelIndex).loadImage(1);
%         subI = MD.getChannel(p.iChanSlave).loadImage(1);
%     end
%     % Do the PIV from FA locations
%     % Track beads displacement in the xy coordinate system
%     xNA=arrayfun(@(x) x.xCoord(1),tracksNA);
%     yNA=arrayfun(@(x) x.yCoord(1),tracksNA);
% 
%     maskAdhesion = refineAdhesionSegmentation(maskFAs>0,mainI,xNA,yNA);
%     indivAdhs = regionprops(bwconncomp(maskAdhesion,4),mainI,'Centroid','Area','BoundingBox','MajorAxisLength','MeanIntensity');
% %     lengthAll = arrayfun(@(x) x.MajorAxisLength,indivAdhs);
%     indivSubAdhs = regionprops(bwconncomp(maskAdhesion,4),subI,'MeanIntensity');
%     orthoLengthAll = arrayfun(@(x) max(x.BoundingBox(3:4)),indivAdhs);
%     areaAll = arrayfun(@(x) x.Area,indivAdhs);
%     intensityAll = arrayfun(@(x) x.MeanIntensity,indivAdhs);
%     intensitySubAll = arrayfun(@(x) x.MeanIntensity,indivSubAdhs);
%     maxFlowSpeed = 5;
%     % Filter out segmentation near image boundary
%     dispMatX = arrayfun(@(x) x.Centroid(1),indivAdhs);
%     dispMatY = arrayfun(@(x) x.Centroid(2),indivAdhs);
%     bwstackImg = true(size(maskFAs));
%     bwstackImg = bwmorph(bwstackImg,'erode',maxFlowSpeed+10);
%     
%     [insideIdx] = maskVectors(dispMatX,dispMatY,bwstackImg); 
% 
%     % We go from obvious FAs that have high enough intensity and area but
%     % exclude too large FAs
%     idxStrongFAs = areaAll>20 & intensityAll>(mean(intensityAll)-std(intensityAll)) ...
%         & insideIdx & intensitySubAll>(mean(intensitySubAll)-1*std(intensitySubAll)); % ...
% %                    & orthoLengthAll<(mean(orthoLengthAll)+3*std(orthoLengthAll));
%     
%     numSFAs = sum(idxStrongFAs);
%     deformAll = zeros(numSFAs,4);
%     indexFAs = find(idxStrongFAs)';
%     addDist = 10;
% 
%     pivPar = [];      % variable for settings
%     pivData = [];     % variable for storing results
%     [pivPar, pivData] = pivParams(pivData,pivPar,'defaults');     
%     [pivData] = pivAnalyzeImagePair(mainI,subI,pivData,pivPar);
%     % show the FA segmentation
%     labelAdh = bwlabel(maskAdhesion,4);
%     % Show boundaries with only in indexFAs
%     labelAdh2 = labelAdh.*ismember(labelAdh,indexFAs);
%     bdryAdh = boundarymask(labelAdh2);
%     figure, imshow(imoverlay(uint8(mainI),bdryAdh,'cyan'), []), hold on
%     quiver(pivData.X,pivData.Y,pivData.U,pivData.V,0,'r')
%     
%     
%     % Decided go from each area with variable area
%     for ii=1:numSFAs
%         curAdh = indivAdhs(indexFAs(ii));
%         deformAll(ii,1:2) = curAdh.Centroid;
%         deformAll(ii,3)=interp2(pivData.X,pivData.Y,pivData.U,curAdh.Centroid(1),curAdh.Centroid(2));
%         deformAll(ii,4)=interp2(pivData.X,pivData.Y,pivData.V,curAdh.Centroid(1),curAdh.Centroid(2));
%     end
%     quiver(deformAll(:,1),deformAll(:,2),deformAll(:,3),deformAll(:,4),0,'y')
%     
%     % Make the transformation matrix (2d projective) out of deformAll
%     idxNoNaN = ~isnan(deformAll(:,3));
%     movingPoints = deformAll(idxNoNaN,1:2)+deformAll(idxNoNaN,3:4);
%     fixedPoints = deformAll(idxNoNaN,1:2);
%     curXform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity'); %'projective');%'affine'); %
%     invCurXform = curXform.invert;
%     disp('Projective transform has been found!')
% 
%     % Apply to the slave channel
%     ref_obj = imref2d(size(subI));
%     imageFileNames = MD.getImageFileNames;
% %     newSubI = imwarp(subI, invCurXform, 'OutputView', ref_obj);
%     
%     % Overwrite the side channel
%     if ~isempty(SDCProc)
%         % back up
%         backupFolder = [SDCProc.outFilePaths_{1,p.iChanSlave} ' Original']; 
%         if ~exist(backupFolder,'dir')
%             mkdir(backupFolder);
%         end
%         % Anonymous functions for reading input/output
%         outFile=@(chan,frame) [SDCProc.outFilePaths_{1,chan} filesep imageFileNames{chan}{frame}];
%         for ii=1:nFrames
%             subI = SDCProc.loadChannelOutput(p.iChanSlave,ii);
%             newSubI = imwarp(subI, curXform, 'OutputView', ref_obj);
%             copyfile(outFile(p.iChanSlave, ii),backupFolder,'f')
%             imwrite(uint16(newSubI), outFile(p.iChanSlave, ii));
%             
%             if ii==1
%                 mainI = SDCProc.loadChannelOutput(p.iChanMaster,ii);
%                 hFig = figure;
%                 hAx  = subplot(1,2,1);
%                 C = imfuse(mainI,subI,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%                 imshow(C,[],'Parent', hAx);
%                 
%                 title('Original result: red, ch 1, green, ch 2')
%                 hAx2  = subplot(1,2,2);
%                 C2 = imfuse(mainI,newSubI,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%                 imshow(C2,[],'Parent', hAx2);
%                 title('Aligned result: red, ch 1, green, adjusted ch 2')
%             end
%         end
%         disp(['Channel has been overwritten in ' SDCProc.outFilePaths_{1,p.iChanSlave}])
%     else
%         backupFolder = [MD.getChannel(p.iChanSlave).getPath ' Original']; 
%         if ~exist(backupFolder,'dir')
%             mkdir(backupFolder);
%         end
%         % Anonymous functions for reading input/output
%         outFile=@(chan,frame) [MD.getChannel(chan) filesep imageFileNames{chan}{frame}];
%         for ii=1:nFrames
%             subI = MD.getChannel(p.iChanSlave).loadImage(ii);
%             newSubI = imwarp(subI, invCurXform, 'OutputView', ref_obj);
%             copyfile(outFile(p.iChanSlave, ii),backupFolder)
%             imwrite(uint16(newSubI), outFile(p.iChanSlave, ii));
%         end
%         disp(['Channel has been overwritten in ' MD.getChannel(chan)])
%     end
% elseif p.doMultimodalRegistration
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

        iMainChan=p.ChannelIndex; iSlaveChan = p.iChanSlave;
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
%% Reading
iReadingCode=4;
for iCurChan = iChanSlave
    iReadingCode = iReadingCode+1;
    [h,w]=size(MD.channels_(iCurChan).loadImage(1));
    imgStack = zeros(h,w,nFrames);
    
    if ~isempty(SDCProc)
        for ii=1:nFrames
            imgStack(:,:,ii)=SDCProc.loadChannelOutput(iCurChan,ii);
        end
    else
        for ii=1:nFrames
            imgStack(:,:,ii)=MD.channels_(iCurChan).loadImage(ii); 
        end
    end
    %% Read force from imgStack
    % get the intensity
    disp('Reading the other channel...')
    tic
    if iReadingCode==5
        fieldsAll = fieldnames(tracksNA);
        unnecFields = setdiff(fieldsAll,{'xCoord','yCoord','startingFrameExtraExtra',...
            'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
            'endingFrame','amp','sigma'});
        essentialTracks = rmfield(tracksNA,unnecFields);
        addedTracksNA = readIntensityFromTracks(essentialTracks,imgStack,iReadingCode); % 5 means ampTotal2 from the other channel
    else
        addedTracksNA = readIntensityFromTracks(addedTracksNA,imgStack,iReadingCode); % 6 means ampTotal3 from the other channel
    end
    tracksAmpTotal = rmfield(addedTracksNA,{'xCoord','yCoord','startingFrameExtraExtra',...
    'startingFrameExtra','startingFrame','endingFrameExtraExtra','endingFrameExtra',...
    'endingFrame','amp','sigma'});
    toc
    %% Add report of mean intensity of the other channel per status
%     indNAs = arrayfun(@(x) sum(x.state(x.startingFrameExtra:x.endingFrameExtra)==2)...
%         /(x.endingFrameExtra-x.startingFrameExtra+1)>0.9, tracksNA);
%     indFCs = arrayfun(@(x) sum(x.state(x.startingFrameExtra:x.endingFrameExtra)==3)...
%         /(x.endingFrameExtra-x.startingFrameExtra+1)>0.9, tracksNA);
%     indFAs = arrayfun(@(x) sum(x.state(x.startingFrameExtra:x.endingFrameExtra)==4)...
%         /(x.endingFrameExtra-x.startingFrameExtra+1)>0.9, tracksNA);
    intensitiesInNAs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==2)),tracksNA,addedTracksNA);
    intensitiesInFCs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==3)),tracksNA,addedTracksNA);
    intensitiesInFAs = arrayfun(@(x,y) nanmean(y.ampTotal2(x.state==4)),tracksNA,addedTracksNA);
    intensityGroup={intensitiesInNAs, intensitiesInFCs, intensitiesInFAs};
%     save(outputFile{2,iChanSlave},'intensityGroup');

    amplitudeInNAs = arrayfun(@(x,y) nanmean(y.amp2(x.state==2)),tracksNA,addedTracksNA);
    amplitudeInFCs = arrayfun(@(x,y) nanmean(y.amp2(x.state==3)),tracksNA,addedTracksNA);
    amplitudeInFAs = arrayfun(@(x,y) nanmean(y.amp2(x.state==4)),tracksNA,addedTracksNA);
    amplitudeGroup={amplitudeInNAs, amplitudeInFCs, amplitudeInFAs};
    save(outputFile{2,iChanSlave},'intensityGroup','amplitudeGroup');
end
%% protrusion/retraction information - most of these are now done in analyzeAdhesionsMaturation
% time after protrusion onset (negative value if retraction, based
% on the next protrusion onset) in frame, based on tracksNA.distToEdge
% First I have to quantify when the protrusion and retraction onset take
% place.

disp('Saving...')
try
    save(outputFile{1,p.ChannelIndex},'tracksAmpTotal'); % the later channel has the most information.
catch
    save(outputFile{1,p.ChannelIndex},'tracksAmpTotal','-v7.3');
end
disp('Done!')
end
%
