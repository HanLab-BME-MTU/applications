function [ output_args ] = GCARunMovie(dataFolder,createNew,startAnal)
%Small function to run different stages of of the GCA segmentation
%was makeMovieDataMultChannels_NoBio until 20150321-

if ( nargin < 3 || isempty(startAnal))
    startAnal = 'all';
end
%
% % Retrieve current location
% % fullPath = which(mfilename);
% % path = fileparts(fullPath);
% % dataFolder=path;


if createNew == 1
    %% Channel creation
    channels = [ 'Channels' filesep 'C1_mCherry']; % misnomer
    
    %channels = ['images3'];
    % Create a channels object
    imgFolder = [dataFolder filesep channels];
    channel = Channel(imgFolder);
    
    
    % Set some channel properties
    % channel.fluorophore_='yfp';
    channel.emissionWavelength_=name2wavelength('GFP')*1e9;
    % channel.imageType_='TIRF';
    %% NOTE need to check if these are the same for Ludo's Experiments
    % NA = 1.4;
    % M = 60;
    % cameraPix= 6.45*1e-6;% in m
    % sigma = getGaussianPSFsigma(NA,M,cameraPix,name2wavelength('mCherry'));
    % channel.psfSigma_ = sigma; % used for the fitting.
    %
    %% MovieData creation
    %upOne = upDirectory(dataFolder,1);
    saveFolder = [upOne filesep 'ANALYSIS'];
    %saveFolder = [dataFolder filesep 'GrowthConeAnalyzer'];
    
    if ~isdir(saveFolder)
        mkdir(saveFolder)
    end
    
    % Constructor needs an array of channels and an output directory (for analysis)
    MD = MovieData(channel,saveFolder);
    
    
    % Set the path where to store the MovieData object.
    MD.setPath(saveFolder);
    MD.setFilename('movieData.mat');
    
    % Run sanityCheck on MovieData.
    % Check image size and number of frames are consistent.
    % Save the movie if successfull
    MD.sanityCheck;
    
    % Set some additional movie properties
    MD.numAperture_=1.4;
    MD.pixelSize_=215; % in nm after binning
    MD.timeInterval_=5;% in sec
    MD.camBitdepth_=16;
    
    %MD.notes_='Created for test purposes';
    
    % Save the movie
    MD.save;
    
else
    % load
    load([dataFolder  filesep 'movieData.mat'])
    saveFolder = MD.outputDirectory_;
end
%% Thresholding - I know this is weird but I run the thresholding here first
% so the Movie Data is already formated correctly (until I have some down
% time to fix it.
% in future this step will be skipped and I will archive my veil/stem steps
% directly in MD.
if (strcmpi(startAnal,'thresh1') || strcmpi(startAnal,'all'))
    %% Initial Threshold (Eventually this will be from neurite package)
    % for now run Hunter's gradient based to compare how performs
    %
    thresProc = ThresholdProcess(MD);
    %
    % % % Save the process in the movie object
    MD.addProcess(thresProc);
    
    % Create a segmentation package
    segPackage = SegmentationPackage(MD);
    
    %  Save the package in the movie object
    MD.addPackage(segPackage);
    %
    % % Associate the threshold process to the package
    MD.packages_{1}.setProcess(1,thresProc);
    
    params = MD.processes_{1}.funParams_;
    params.MethodIndx = 4;
    
    
    parseProcessParams(MD.processes_{end},params);
    %
    % % % Run the process
    MD.processes_{1}.run(4); % run the gradient based
    
    MD.save
end

%% Run Neurite Orientation Estimation
if (strcmpi(startAnal,'neuriteOrient') || strcmpi(startAnal,'all'));
    % Use the shade corrected images for the detection if you have them.
    GCAgetNeuriteOrientMovie(MD);
end

if (strcmpi(startAnal,'fixNeuriteOrient') || strcmpi(startAnal,'all'));
    GCAneuriteOrientConsistencyCheckMovie(MD);
end

%% % RunVeilStem
if strcmpi(startAnal,'veilStem') || strcmpi(startAnal,'all')
    GCAReconstructVeilStemMovie(MD);
end


%% do the switch ... my masks for gradient masks - quick fix until I can write my own process for veil/stem est.
if (strcmpi(startAnal,'BodySwitch') || strcmpi(startAnal,'all'));
    for iCh = 1:1
        
        % search for a refinement folder
        refineDataMat =  [ MD.outputDirectory_ filesep 'neurite_veilStem_refinements' filesep ...
            'analInfoWithScaleSave.mat'];
        if exist(refineDataMat,'file')== 2;
            load(refineDataMat);
        else
            load([MD.outputDirectory_ filesep  'neurite_body_masks' filesep ...
                'Neurite_Body_Masks_Channel_' num2str(iCh) filesep 'analInfoTestSave.mat']);
        end
        
        maskDir = [ saveFolder filesep 'masks' filesep 'masks_for_channel_' num2str(iCh)];
        
        % old mask folder
        %maskDirGrad = [saveFolder filesep 'masks' filesep 'masks_for_channel_' num2str(iCh) 'old'];
        
        
        %copyfile(maskDir,maskDirGrad);
        rmdir(maskDir,'s');
        mkdir(maskDir);
        
        makeNeuriteBodyMaskFolder(analInfo,maskDir,0,1,1,MD);
    end
    display('Finished Neurite Body Masks');
end % strcmpi 'BodySwitch'

%% Run Protrusion Software... make my own file to save the output in the folder
if (strcmpi(startAnal,'prot') || strcmpi(startAnal,'all') )
    
    
    
    prot = ProtrusionProcess(MD);
    MD.addProcess(prot);
    idxProtProc = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),MD.processes_));
    idxProtProc = idxProtProc(end); %
    protParams = MD.processes_{idxProtProc}.funParams_;
    protParams.ChannelIndex = 1; % use the DONOR CHANNEL
    protParams.SegProcessIndex = 1; % NEED to make my own neurite body est process
    parseProcessParams(MD.processes_{idxProtProc},protParams);
    %
    %
    %
    %
    MD.processes_{idxProtProc}.run(1);
    %
    % % associate it with the windowing package
    windPack = WindowingPackage(MD);
    MD.addPackage(windPack);
    % %
    MD.packages_{end}.setProcess(1,prot);
    %
    MD.save;
end
%% Run Reconstruction (organize input)
iCh = 1; %%% FIX MEEEE !!!

if strcmpi(startAnal,'reconstruct') || strcmpi(startAnal,'all')
    %idxShade = find(cellfun(@(x) sum(strcmpi(x.name_,'Shade Correction')),MD.processes_));
    
    projList{1} = MD.channels_(1).channelPath_;
    % projData.imDir = MD.channels_(1).channelPath_;
    idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),MD.processes_));
    restart =2; % zero if do not restart (currently numbers - need to fix input)
    saveDir = [MD.outputDirectory_ filesep 'filopodia_reconstruct'];
    saveDir = [saveDir filesep 'Filopodia_Reconstruct_Channel_' num2str(1)];
    if (exist(saveDir) && (restart ==1 || restart ==2))
        display('Restarting Reconstruction');
        load([saveDir filesep 'analInfoTestSave.mat']);
    else
        
        if ~isdir(saveDir)
            mkdir(saveDir)
        end
        % search for a refinement folder
        refineDataMat =  [ MD.outputDirectory_ filesep 'neurite_veilStem_refinements' filesep ...
            'analInfoWithScaleSave.mat'];
        if exist(refineDataMat,'file')== 2;
            load(refineDataMat);
            
        else
            load([MD.outputDirectory_ filesep  'neurite_body_masks' filesep  'Neurite_Body_Masks_Channel_' num2str(iCh) filesep 'analInfoTestSave.mat']);
        end
        % load([MD.outputDirectory_ filesep 'neurite_body_masks' filesep 'Neurite_BodyMask_Channel_' num2str(iCh) filesep 'params.mat']);
    end
    projData.imDir = projList{iCh}; % just use the same project list as for the
    % body estimation
    if iCh ==1
        % run orientation calcs
        protS =  load(MD.processes_{idxProt(end)}.outFilePaths_);
    else
        protS = []; % don't do orientation calcs
    end
    GCAReconstructFilopodia(projData.imDir,saveDir,1,0,restart,analInfo,35,protS);
    
end % strcmpi reconstruct
%% GCAfitFilopodia
if strcmpi(startAnal,'fitFilo') || strcmpi(startAnal,'all')
    GCAfitFilopodiaMovie(MD);
end


% %% run windowing...(note you can maybe do this later to save time is don't need the protrusion
% % cal right away.
% % Run windowing Note here figure out how to set the start window
% % automatically at the neurite entrance
if strcmpi(startAnal,'wind') || strcmpi(startAnal,'all')
    wind = WindowingProcess(MD);
    MD.addProcess(wind);
    windParams = MD.processes_{end}.funParams_;
    windParams.ChannelIdx = 1; % use the DONOR channel
    windParams.SegProcessIndex = 1;
    % start the windows at the neurite entrance point so more visually
    % appealing
    % load analInfo
    % check for the refinement directory
    veilStemRefine = [MD.outputDirectory_ filesep 'neurite_veilStem_refinements' filesep 'analInfoWithScaleSave.mat'];
    if exist(veilStemRefine,'file')~=0
        load(veilStemRefine)
    else
        
        analInfoDir = [MD.outputDirectory_ filesep 'neurite_body_masks' filesep 'Neurite_Body_Masks_Channel_1' ] ;
        
        load([analInfoDir filesep 'analInfoTestSave.mat']);
    end
    
    windParams.ParaSize = 3 ;
    SPIdx = analInfo(1).idxEnterNeurite;
    %
    [SPy,SPx] = ind2sub(MD.imSize_,SPIdx);
    % windParams.StartPoint = [SPy SPx];
    windParams.StartPoint = [SPx SPy]; % hunter's input is xy
    windParams.PerpSize = 3;
    windParams.MinSize = 500;
    %windParams.MethodName = 'ConstantNumber';
    %windParams.ReInit = 40;
    parseProcessParams(MD.processes_{end},windParams);
    %
    MD.processes_{end}.run(1);
    idxWind = find(cellfun(@(x) sum(strcmpi(x.name_,'Windowing')),MD.processes_));
    MD.packages_{end}.setProcess(2,wind);
    
    %
    MD.save;
    
    %% Run protrusion sampling process
    protSamp = ProtrusionSamplingProcess(MD);
    MD.addProcess(protSamp);
    MD.processes_{end}.run(1);
    MD.packages_{end}.setProcess(3,protSamp)
    
    MD.save;
end


