function [tracksNA,idGroups] = assembleTracksFromFAPackage(MD)
% function [tracksNA] = assembleTracksFromFAPackage(MD) assemble tracksNA
% from MD that has run FA package. 
% Sangyoon Han July 2023

%% Input

try
    iProc = MD.getProcessIndex('InitialRiseTimeLagCalculationProcess',1,0);
    %If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(MD.processes_)+1;
        MD.addProcess(InitialRiseTimeLagCalculationProcess(MD));                                                                                                 
    end
    timeLagProc = MD.processes_{iProc};
catch
    iFAPack =  MD.getPackageIndex('FocalAdhesionPackage');
    FAPack = MD.getPackage(iFAPack);
    timeLagProc = FAPack.processes_{11};
end
p = parseProcessParams(timeLagProc); %,paramsIn);
iBeadChan = 1; % might need to be updated based on asking TFMPackage..
p.measureAmp2rates = true; % Should be replaced with a user input.

% ip.addParamValue('chanIntensity',@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse(MD,showAllTracks,plotEachTrack,varargin{:});
% iChanMaster=p.ChannelIndex;

%% Reading tracks from master channel
% Use the previous analys
%% Reading tracks from master channel
% Use the previous analysis folder structure
% It might be good to save which process the tracksNA was obtained.
disp('Reading tracksNA ...')
iFAPack =  MD.getPackageIndex('FocalAdhesionPackage');
FAPack = MD.getPackage(iFAPack);

tic
try
    iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
    adhAnalProc = MD.getProcess(iAdhProc);
catch
    adhAnalProc = FAPack.processes_{7};
end
% Doublecheck the p.ChannelIndex
p.ChannelIndex = adhAnalProc.funParams_.ChannelIndex;
tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA');
% numChans = numel(p.ChannelIndex);
%% Doublecheck if amp has similar range as amp2 and read it again at least for the extra part
% find one examplary track that has intermediate starting point
startingFrameExtraAll = arrayfun(@(x) x.startingFrameExtra,tracksNA);
startingFrameExtraExtraAll = arrayfun(@(x) x.startingFrameExtraExtra,tracksNA);
intermedTrackIDs = find(startingFrameExtraExtraAll<(startingFrameExtraAll-10));
if ~isempty(intermedTrackIDs)
    trackInspected = tracksNA(intermedTrackIDs(1));
    % See if amp is NaN before startingFrameExtra (after
    % startingFrameExtraExtra)
    if isnan(trackInspected.amp(trackInspected.startingFrameExtraExtra+1))
        % then we need to read this amp again
        disp('Reading image stack again to read extra time regime of amp...')
        tic
        imgStack = getAnyStacks(MD);
        toc
        disp('Reading from tracks'); tic
        tracksNA = readIntensityFromTracks(tracksNA,imgStack,1,'extraReadingOnly',true); toc;
        disp('Saving the tracks...'); tic
        numTracks=numel(tracksNA);
        fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
        numStr = @(trackNum) num2str(trackNum,fString);
        trackFolderPath = [adhAnalProc.funParams_.OutputDirectory filesep 'trackIndividual'];
        trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
        
        parfor k=1:numTracks
            curTrack=tracksNA(k);
            parsave(feval(trackIndPath,k),curTrack)
        end
        toc
    end
else
    disp('All tracks have the same or similar starting points')
end

%% Now we have to combine this with readings from step 9 and 10
iTheOtherProc=9; 
theOtherReadProc=FAPack.processes_{iTheOtherProc};

if ~isempty(theOtherReadProc)
    ampObj = load(theOtherReadProc.outFilePaths_{1,p.ChannelIndex},'tracksAmpTotal'); % the later channel has the most information.
    tracksAmpTotal = ampObj.tracksAmpTotal;
%     tracksAmpTotal = theOtherReadProc.loadChannelOutput(p.ChannelIndex);
    
    if isfield(tracksAmpTotal,'ampTotal2')
        [tracksNA(:).ampTotal2] = tracksAmpTotal.ampTotal2;
    end
    if isfield(tracksAmpTotal,'amp2')
        [tracksNA(:).amp2] = tracksAmpTotal.amp2;
    end
    if isfield(tracksAmpTotal,'ampTotal3')
        [tracksNA(:).ampTotal3] = tracksAmpTotal.ampTotal3;
    end
end

%% Combinging from Step 10
iForceRead=10;
forceReadProc=FAPack.processes_{iForceRead};
if ~isempty(forceReadProc)
    forceReadObj = load(forceReadProc.outFilePaths_{1,p.ChannelIndex},'tracksForceMag'); % the later channel has the most information.
    tracksForceMag = forceReadObj.tracksForceMag;
    idxTracksObj = load(forceReadProc.outFilePaths_{2,p.ChannelIndex},'idxTracks');
    if ~isfield(idxTracksObj,'idxTracks')
        idxTracksObj = load(forceReadProc.outFilePaths_{6,p.ChannelIndex},'idxTracks');
    end
    idxTracks = idxTracksObj.idxTracks;
    tracksNA = tracksNA(idxTracks);
    if isfield(tracksForceMag,'forceMag')
        [tracksNA(:).forceMag] = tracksForceMag.forceMag;
    end
end
toc
%% Reading classes
disp('Reading idsClassified ...')
try
    try
        iClaProc = MD.getProcessIndex('AdhesionClassificationProcess');
        classProc = MD.getProcess(iClaProc);
    catch
        classProc = FAPack.processes_{8};
    end
    % numChans = numel(p.ChannelIndex);
    idsClassified = load(classProc.outFilePaths_{4,p.ChannelIndex});
    idGroup1 = idsClassified.idGroup1;
    idGroup2 = idsClassified.idGroup2;
    idGroup3 = idsClassified.idGroup3;
    idGroup4 = idsClassified.idGroup4;
    idGroup5 = idsClassified.idGroup5;
    idGroup6 = idsClassified.idGroup6;
    idGroup7 = idsClassified.idGroup7;
    idGroup8 = idsClassified.idGroup8;
    idGroup9 = idsClassified.idGroup9;
catch
    disp('No Classified groups. Using no classification...')
    % Potentially we can use the dedactically chosen classes here (e.g.
    % shown in the code used for Tristan's movie analysis.
    idGroup1 = []; idGroup2 = []; idGroup3 = []; idGroup4 = [];
    idGroup5 = []; idGroup6 = []; idGroup7 = []; idGroup8 = []; idGroup9 = [];
end    
if ~isempty(forceReadProc)
    idGroup1 = idGroup1(idxTracks);
    idGroup2 = idGroup2(idxTracks);
    idGroup3 = idGroup3(idxTracks);
    idGroup4 = idGroup4(idxTracks);
    idGroup5 = idGroup5(idxTracks);
    idGroup6 = idGroup6(idxTracks);
    idGroup7 = idGroup7(idxTracks);
    idGroup8 = idGroup8(idxTracks);
    idGroup9 = idGroup9(idxTracks);
else
    disp('Traction reading was not done. No further filtering...')
end
idGroups = {idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};

disp('Track reading and integration are done.')
end

