function [filoBranch] = GCAReconstructFilopodiaMovie(movieData,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the filoBranchInfo structure to.
%       If not input, the filoBranch will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "VI_filopodiaBranch_reconstruction"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the filopodia
%       reconstruct
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input). If not input, the
%       backbone information will be calculated from the channels (ie raw
%       images)
%
%  GCAReconstructFilopodia specific functions
%
% OUTPUT: (see main function GCAReconstructFilopodia- for details):
%       filoBranch Rx1 structure array where R = the number of frames
%       with added field .filoInfo a Rx1 structure where R is the number of
%       filopodia detected in that frame
%       .filoInfo has many subfields that store information regarding the
%       filopodias coordinates,orientation,grouping(for branching) etc.


%%
%% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VI_filopodiaBranch_reconstruction'];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1); % need separate ChannelIndex for veil and filo
ip.addParameter('ChannelIndexVeil',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

% Archiving the software version
ip.addParameter('getGITHashTag',false);

% Specific
% PARAMETERS
ip.addParameter('TSOverlaysReconstruct',false);
ip.addParameter('TSOverlaysRidgeCleaning',false); 

% Steerable Filter
ip.addParameter('FilterOrderFilo',4,@(x) ismember(x,[2,4]));
ip.addParameter('FiloScale',1.5);

% Estimate background to localize region of interest
ip.addParameter('filterBackEst',true); % flag to estimate high confidence
% background using the image intensity histogram : < mean + 2std
% is considered background
ip.addParameter('dilateLocalRegion',false); % flag to dilate further the
% local region of interest estimation (ie decrease the background
% estimation)
ip.addParameter('LRDilRad',10); % dilation radius of the structuring element
% applied to inital guess of the localized region of interest.


% Cleaning Response of Steerable Filter
ip.addParameter('multSTDNMSResponse',3);
ip.addParameter('minCCRidgeOutsideVeil',3);
ip.addParameter('filterBasedOnVeilStemAttachedDistr',true);

% Linking Parameters Embedded
%ip.addParameter('geoThreshEmbedded',0.9,@(x) isscalar(x));
ip.addParameter('geoThreshEmbedded',0.5,@(x) isscalar(x));
%ip.addParameter('maxRadiusLinkOutsideVeil',10);
ip.addParameter('maxRadiusLinkEmbedded',10);
ip.addParameter('curvBreakCandEmbed',0.05,@(x) isscalar(x));

% Linking Parameters Candidate Building
ip.addParameter('maxRadiusLink',5); %
ip.addParameter('geoThresh',0.9, @(x) isscalar(x));

% Linking Parameters Traditional Filopodia/Branch Reconstruct
ip.addParameter('maxRadiusConnectFiloBranch',15); % change default to 15
ip.addParameter('geoThreshFiloBranch',0.5);

% Option to detect Embedded actin bundles (for LifeAct and actin staining
% only)
ip.addParameter('detectEmbedded',true)

ip.addParameter('rotateVeilStemNormals',true);


ip.parse(varargin{:});
params = ip.Results;

%% Init:
nFrames = movieData.nFrames_-1; % to match the protrusion input length...
nChan = numel(params.ChannelIndex);
channels = params.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);

%% Loop for each channel
for iCh = 1:nChan
    
    display(['Reconstructing Filopodia Channel ' num2str(channels(iCh))]);
    
    %% Get Start and End Frames Based on Restart Choice
    
    outDir = [ip.Results.OutputDirectory filesep  'Channel_' num2str(channels(iCh))];
    
    if ~isdir(outDir);
        mkdir(outDir);
    end
    
    reconstructFile = [outDir filesep 'filoBranch.mat'];
    
    % If file exists
    if  exist(reconstructFile,'file')==2;
        load(reconstructFile) % load the file
        display('Loading Previously Run Filopodia Reconstructions ');
        if strcmpi(params.StartFrame,'auto')
            startFrame = numel(filoBranch)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = params.StartFrame; % use user input
            display(['Manual Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        startFrame = 1;
        display('No Filopodia Reconstruction Folder Found: Creating and Starting at Frame 1');
        
    end % exist(orientFile,'file') == 2
    
    if startFrame ~= 1
        load([outDir filesep 'params.mat']);
    end
    % startFrame = params.StartFrame;
    
    if strcmpi(params.EndFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = params.EndFrame;
        display(['Manual End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    %% Load veilStem from veil/stem folder
    
    load([movieData.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep ...
        'IV_veilStem_length' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep 'veilStem.mat']);
    
    % load( [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ...
    %    'StepsToReconstruct' filesep 'III_veilStem_reconstruction' ...
    %   filesep 'Channel_' num2str(channels(iCh)) filesep 'veilStem.mat']);
    
    %%
    % get the list of image filenames
    if params.ProcessIndex == 0
        imDir = movieData.channels_(channels(iCh)).channelPath_;
    else
        imDir = movieData.proceses_{params.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif','memo',imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
    
    %% Test if the protrusion process was run
    idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),movieData.processes_));
    if ~isempty(idxProt)
        %         load the three protrusion process associated cell structures,
        %         'normals','protrusion','smoothedEdge'
        load(movieData.processes_{idxProt(end)}.outFilePaths_);
        %         If for some reason there is more than one protrusion process
        %         associated with the data, tell the user that your are using the
        %         most recently run
        if length(idxProt) >1
            display(['There was more than one veil protrusion process associated with ' ...
                'this movie: local filopodia to veil calculations will be '...
                'will be performed using the most recently run protrusion process']);
        end % if length(protS)
    else % no veil protrusion vectors were found
        display(['The veil protrusion process was not found for this movie :'...
            'The filopodia reconstruction will continue, but no local filopodia orientation ' ...
            'calculations relative to the veil will be performed']);
        normals = [];
        smoothedEdge = [];
    end % ~isempty(idxProt)
    %
    %load([movieData.outputDirectory_ filesep 'protrusion' filesep 'protrusion_vectors.mat'] );
    
    %% Start Movie Loop %%%%
    for iFrame = startFrame:endFrame
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        veilStemMaskC = veilStem(iFrame).finalMask;
        
        % All pieces of veil/stem are stored: make sure only considering the largest
        % connected component
        veilStemMaskC = logical(getLargestCC(veilStemMaskC));
        
        protrusionC.normal = normals{iFrame};
        protrusionC.smoothedEdge = smoothedEdge{iFrame};
        if ip.Results.rotateVeilStemNormals
            EPLead = veilStem(iFrame).endPointLeadingProt;
            LPIndices = veilStem(iFrame).neuriteLongPathIndices;
            idxEnter = veilStem(iFrame).idxEnterNeurite;
        else
            EPLead = [];
            LPIndices = [];
            idxEnter = [];
        end
        
        [filoBranchC,TSFigs,TSFigsRecon,errorHist] =  GCAReconstructFilopodia(img,veilStemMaskC,protrusionC,EPLead,LPIndices,idxEnter,params);
        if errorHist == 1
            display(['Problem With Histogram Frame' num2str(iFrame)]);
        end
        %% Plot the results.
        if (ip.Results.TSOverlaysReconstruct || ip.Results.TSOverlaysRidgeCleaning)
            display('Saving Trouble Shoot Overlays')
            for iFig = 1:length(TSFigs)
                if ~isdir([outDir filesep TSFigs(iFig).group filesep num2str(iFig,'%02d') TSFigs(iFig).name]);
                    mkdir([outDir filesep TSFigs(iFig).group filesep num2str(iFig,'%02d') TSFigs(iFig).name]);
                end
            end
            type{1} = '.fig';
            type{2} = '.tif';
            
            
            if ~isempty(TSFigs)
                for iType = 1:numel(type)
                    arrayfun(@(x) saveas(TSFigs(x).h,...
                        [outDir filesep TSFigs(x).group filesep num2str(x,'%02d') TSFigs(x).name filesep num2str(iFrame,'%03d') type{iType}]),1:length(TSFigs));
                end
            end
        end
        
        if ip.Results.TSOverlaysReconstruct
            for iFig = 1:length(TSFigsRecon)
                cDir = [outDir filesep TSFigsRecon(iFig).group  filesep 'ReconIter' num2str(TSFigsRecon(iFig).ReconIt,'%02d') ...
                    filesep TSFigsRecon(iFig).name];
                if ~isdir(cDir);
                    mkdir(cDir);
                end
            end
            type{1} = '.fig';
            type{2} = '.tif';
            
            
            if ~isempty(TSFigsRecon)
                for iType = 1:numel(type)
                    
                    arrayfun(@(x) saveas(TSFigsRecon(x).h,...
                        [outDir filesep TSFigsRecon(x).group  filesep 'ReconIter' num2str(TSFigsRecon(x).ReconIt,'%02d') ...
                        filesep TSFigsRecon(x).name filesep num2str(iFrame,'%03d') type{iType}]),1:length(TSFigsRecon));
                    
                end
                arrayfun(@(x) saveas(TSFigsRecon(x).h,...
                    [outDir filesep TSFigsRecon(x).group  filesep 'ReconIter' num2str(TSFigsRecon(x).ReconIt,'%02d') ...
                    filesep TSFigsRecon(x).name filesep num2str(iFrame,'%03d') '.eps'],'psc2'),1:length(TSFigsRecon));
            end
            
            %if ~isdir([outDir filesep TSFigs(iFig).group filesep num2str(TSFigs(iFig).iter),
        end
        
        close all
        
        if ip.Results.getGITHashTag
            [~,hashTag] =  gcaArchiveGetGitHashTag;
        else
            hashTag = NaN;
        end
        filoBranchC.hashTag = hashTag;
        filoBranchC.timeStamp = clock;
        filoBranch(iFrame) = filoBranchC;
        p(iFrame) = params;
        save( [outDir filesep 'filoBranch.mat'], 'filoBranch','-v7.3');
        
        
        save([outDir filesep 'params.mat'],'p');
        display(['Finished Reconstructing Filopodia for Frame ' num2str(iFrame) ' for ' movieData.outputDirectory_])
    end % iFrame   
end   % for iCh
end % function 