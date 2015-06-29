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
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

% Specific
% PARAMETERS

% Steerable Filter 
ip.addParameter('FilterOrderFilo',4,@(x) ismember(x,[2,4]));
ip.addParameter('FiloScale',1.5);

% Cleaning Response of Steerable Filter 
ip.addParameter('multSTDNMSResponse',3);
ip.addParameter('minCCRidgeOutsideVeil',3);

% Linking Parameters Embedded 
ip.addParameter('geoThreshEmbedded',0.9,@(x) isscalar(x)); 
%ip.addParameter('maxRadiusLinkOutsideVeil',10);
ip.addParameter('maxRadiusLinkEmbedded',10);
ip.addParameter('curvBreakCandEmbed',0.05,@(x) isscalar(x)); 

% Linking Parameters Candidate Building
ip.addParameter('maxRadiusLink',5); 
ip.addParameter('geoThresh',0.9, @(x) isscalar(x));  

% Linking Parameters Traditional Filopodia/Branch Reconstruct
ip.addParameter('maxRadiusConnectFiloBranch',5); 



ip.parse(varargin{:});
p = ip.Results;
%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
channels = p.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);

%%
%% Loop for each channel
for iCh = 1:nChan
    
    display(['Reconstructing Filopodia Channel ' num2str(channels(iCh))]);
    
    %% Get Start and End Frames Based on Restart Choice
    
    
    
    saveDir = [movieData.outputDirectory_ filesep ...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VI_filopodiaBranch_reconstruction' ...
        filesep 'Channel_' num2str(channels(iCh))];
    
    if ~isdir(saveDir);
        mkdir(saveDir);
    end
    
    reconstructFile = [saveDir filesep 'filoBranch.mat'];
    
    % If file exists
    if  exist(reconstructFile,'file')==2;
        load(reconstructFile) % load the file
        display('Loading Previously Run Filopodia Reconstructions ');
        if strcmpi(p.startFrame,'auto')
            startFrame = numel(analInfo)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = p.restart.startFrame; % use user input
            display(['Manual Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        startFrame = 1;
        display('No Filopodia Reconstruction Folder Found: Creating and Starting at Frame 1');
        
    end % exist(orientFile,'file') == 2
    
    
    
    
    
    if strcmpi(p.EndFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = p.EndFrame;
        display(['Manual End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    %% Load veilStem from veil/stem folder
    
    load([movieData.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep ...
        'IV_veilStem_length' filesep 'Channel_' num2str(channels(iCh)) filesep 'veilStem.mat']);
    
   % load( [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ...
    %    'StepsToReconstruct' filesep 'III_veilStem_reconstruction' ...
     %   filesep 'Channel_' num2str(channels(iCh)) filesep 'veilStem.mat']);
    
    
     
    
    %%
    % get the list of image filenames
    if p.ProcessIndex == 0
        imDir = movieData.channels_(channels(iCh)).channelPath_;
    else
        imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
    
    %% Test if the protrusion process was run
    idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),movieData.processes_));
    if ~isempty(idxProt)
        % load the three protrusion process associated cell structures,
        % 'normals','protrusion','smoothedEdge'
        load(movieData.processes_{idxProt(end)}.outFilePaths_);
        % If for some reason there is more than one protrusion process
        % associated with the data, tell the user that your are using the
        % most recently run
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
    
    %% Start Movie Loop %%%%
    for iFrame = startFrame:endFrame
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        veilStemMaskC = veilStem(iFrame).finalMask;
        
%         if ~isempty(normals)
%             protrusionC.normal = normals{iFrame};
%             protrusionC.smoothedEdge = smoothedEdge{iFrame};
%             x{1} = protrusionC; % make a cell for parsing input. 
%             [filoBranch,TSFigs] = GCAReconstructFilopodiaWorkingInput(img,veilStemMaskC,x,p);
%         else
%             [filoBranch,TSFigs] = GCAReconstructFilopodiaWorkingInput(img,veilStemMaskC,p);
%         end
%         
        
        protrusionC.normal = normals{iFrame}; 
        protrusionC.smoothedEdge = smoothedEdge{iFrame}; 
        EPLead = veilStem(iFrame).endPointLeadingProt; 
        LPIndices = veilStem(iFrame).neuriteLongPathIndices;
    
        GCAReconstructFilopodiaWorkingInput(img,veilStemMaskC,protrusionC,EPLead,LPIndices,p); 
        
        
        
        
        
        
        % quick fix for the plots is to just make the frame number an input for not 20140812
        hashTag =  gcaArchiveGetGitHashTag;
        filoBranch.hashTag = hashTag; % make sure to add the hash tag first so the structure is similar (or initiate in begin)
        filoBranch(iFrame) = filoBranch;
        save( [outDir filesep 'filoBranch.mat'],'filoBranch');
    end % iFrame
    
end   % for iCh
end


