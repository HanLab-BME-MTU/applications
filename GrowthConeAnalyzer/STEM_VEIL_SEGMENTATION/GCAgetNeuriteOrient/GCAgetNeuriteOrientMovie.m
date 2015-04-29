function [ backboneInfo] = GCAgetNeuriteOrientMovie(movieData,paramsIn)
% GCAgetNeuriteOrientMovie.m :  movieWrapper function for GCAgetNeuriteOrient
% STEP I in Veil/Stem estimation of GCA package
%
% %%PERSONAL NOTE: was getNeuriteOrientationMovie.m until 20140529)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "neurite_orientation_estimations"
%       % PERSONAL NOTE : You might need
%       to truncate these names for the windows paths. %
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation
%       estimation.
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
%  See GCAgetNeuriteOrient specific input
%
%       ('BBScale -> Positive integer scalar or vector of scales to use for
%       the steerable filter estimation of the backbone. Note if a vector
%       is specified repsonses for all scales are calculated and the
%       scale with the largest steerable filter
%       response at each point is chosen for the final backbone estimation:
%       Default: 5-10 (in pixels))
%
%       (FilterOrderBB -> scalar 2 or 4: Default is 4 )
%
%       ('MaxRadiusLargeScaleLink  -> Positive scalar: Defines the Maximum 
%        Search Radius for Linking Large-Scale Ridges by the KD tree) 
%
%       ('threshNonMaxSupResponse -> percentile value of population below
%       which to threshold out Default: 25th percentile 
% 
%       (
%
%
%
%
%       (plots -> if true will make troubleshoot movie
%       folders with the larger outputDirectory for the given channel
%       'RidgeCandBeforeAfterClean'
%       'BeforeAndAfterConnect'
%       'CandSeeds'
%       'ScaleMaps'
%
%
% OUTPUT: (see main function GCAgetNeuriteOrient- for details):
%       backboneInfo Nx1 structure array where N = the number of frames
%       with fields
%       .backboneSeedMask   : the current candidate seed mask
%       .coordsEnterNeurite : the current coords of the entrance point of
%                             the neurite
%       .candidateSeeds     :
%       .ridgeScaleMap      :
% (PERSONAL NOTE: Maybe make these optional as a verbose option...)
%       .beforeConnect
%       .afterConnect
%       .maxNMSLarge
%       .cutoffResponse
%% Input

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    % Generic
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'GCSegmentation' filesep 'neurite_orientation_estimations'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
    paramsIn.restart.startFrame = 'auto';   % 'auto'
    paramsIn.restart.endFrame = 'auto'; % 'auto' or number, default auto.    
    paramsIn.plots = 1;
    
   % Specific
    paramsIn.BBScale = 5:10; % PERSONAL NOTE: not sure if this is optimal ...
    paramsIn.FilterOrderBB = 4;  
    paramsIn.MaxRadiusLargeScaleLink = 10; % Defines the Search Radius for Linking Linear Pieces by the KD tree
    paramsIn.ThreshNonMaxSuppResponse = 25; 
    
end

% FOR WHEN MAKE PROCESS
%Get the indices of any previous mask refinement processes from this function
% iProc = movieData.getProcessIndex('GetNeuriteOrientationProcess',1,0);
%
% %If the process doesn't exist, create it
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(GetNeuriteOrientationProcess(movieData,movieData.outputDirectory_));
% end
%
% %Parse input, store in parameter structure
% p = parseProcessParams(movieData.processes_{iProc},paramsIn);

p = paramsIn;

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);





%% Loop for each channel
for iCh = 1:nChan
    
    display(['Finding Neurite Orientation Channel ' num2str(iCh)]);
    %% Get Start and End Frames Based on Restart Choice
    
    
    
    orientFile = [movieData.outputDirectory_ filesep 'neurite_orientation_estimations' filesep 'Neurite_Backbone_Seed_Channel_1' filesep 'backboneInfo.mat'];
    
    
    % If file exists
    if  exist(orientFile,'file')==2;
        load(orientFile) % load the file
        display('Loading Previously Run Orientation Estimations');
        if strcmpi(p.restart.startFrame,'auto')
            startFrame = numel(backboneInfo)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = p.restart.startFrame; % use user input
            display(['Manual Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        
        startFrame = 1;
        
        display('No Neurite Orientation Folder Found: Creating and Starting at Frame 1');
        
        
    end % exist(orientFile,'file') == 2
    
    if strcmpi(p.restart.endFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Finding Neurite Orientations From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = p.restart.endFrame;
        display(['Manual End: Finding Neurite Orientations From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    
    %%
    
    % make final output dir where backboneInfo will be saved
    saveDir =  [p.OutputDirectory filesep 'Neurite_Backbone_Seed_Channel_' num2str(iCh)];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % get the list of image filenames
    if p.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
    else
        imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
    
    % initiate backbone structure
    % backbone = struct([]); % can't initiate structure here - the fields
    % don't match when you add to it- figure out how to better initiate (
    % might need all the fields)
    
    pInNeurite = p;
    pInNeurite.OutputDirectory = saveDir ;
    save([saveDir filesep 'paramsIn.mat'],'pInNeurite');
    
    %% Start Movie Loop %%%%
    for iFrame = startFrame:endFrame
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        [backboneFrame] = GCAgetNeuriteOrient(img,pInNeurite,iFrame); % quick fix for the plots is to just make the frame number an input for not 20140812
        hashTag =  gcaArchiveGetGitHashTag;
        backboneFrame.hashTag = hashTag; % make sure to add the hash tag first so the structure is similar (or initiate in begin)
        backboneInfo(iFrame) = backboneFrame;
        save( [saveDir filesep 'backboneInfo.mat'],'backboneInfo');
        save( [saveDir filesep 'params.mat'],'params.mat'); 
    end % iFrame
    
end   % for iCh
end


