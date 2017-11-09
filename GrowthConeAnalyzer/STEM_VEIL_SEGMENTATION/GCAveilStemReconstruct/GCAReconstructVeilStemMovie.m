function [ output_args ] = GCAReconstructVeilStemMovie( movieData,varargin)
%%% GCAveilStemReconstructMovie.m :  movieDataWrapper function for
%%% GCAveilStemReconstruct.m
% STEP III in Veil/Stem estimation of GCA Segmentation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData (REQUIRED)  - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
% Generic Fields:
%       ('Input/Output Fields Needed for Wrapper' -> Possible Values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "neurite_orientation"
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
%       ('Restart'  -> structure with fields = 'auto';   % 'auto'
%       paramsIn.restart.endFrame = 'auto'; % 'auto' or number, default auto.
%       paramsIn.plots = 1;
%
%  SPECIFIC INPUT
%  See GCAveilStemReconstruct.m for details.
%% %% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction'];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));
ip.addParameter('TSMovie',false,@(x) islogical(x)); 

ip.addParameter('LocalThresholdPatchSize',75,@(x) isscalar(x));
ip.addParameter('DiskSizeLarge',6,@(x) isscalar(x));
ip.addParameter('DiskSizeSmall',3,@(x) isscalar(x)); 

ip.parse(varargin{:});
p = ip.Results;

%% Initiate
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;

%% Start Channel Wrapper

for iCh = 1:nChan
    
    display(['Reconstructing the Veil/Stem for Channel ' num2str(iCh)]);
    

    
    %% Get Start and End Frames Based on Restart Choice
    
    % make final output dir where backboneInfo will be saved
    outDirC =  [ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh)];
    
    if ~isdir(outDirC)
        mkdir(outDirC);
    end
    
    if p.ProcessIndex == 0
        imgDir =  movieData.channels_(p.ChannelIndex).channelPath_; % currently mainly needed if you do the local thresholding/ otherwise just overlay
    else
        imgDir = movieData.proccesses_(p.ProcessIndex).outFilePaths_{p.ChannelIndex};
    end
    
    % collect images and initiate
    [listOfImages] = searchFiles('.tif',[],imgDir ,0);
    
    
    
%% Get Restart Information 
    veilStemFile = [outDirC filesep 'veilStem.mat'];
    
    % If veilStemInfo already exists load it
    if  exist(veilStemFile,'file')==2;
        load(veilStemFile) % load the file
        display('Loading Previously Run VeilStem Estimations');
        if strcmpi(ip.Results.StartFrame,'auto')
            startFrame = numel(veilStem)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = ip.Results.StartFrame; % use user input
            display(['Manual Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        
        startFrame = 1;
        
        display('No Veil/Stem Folder Found: Creating and Starting at Frame 1');
        
        
    end % exist(orientFile,'file') == 2
    
    % if restarting add saved parameters 
    if startFrame ~= 1 
        load([outDirC filesep 'params.mat']); 
    end 
    
    
    if strcmpi(ip.Results.EndFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = ip.Results.EndFrame;
        display(['Manual End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
%%    Load Backbone Information
    load([ movieData.outputDirectory_ filesep...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'II_neurite_orientation_refinements' filesep 'Channel_' num2str(iCh)...
        filesep 'backboneInfoFix.mat']);
  pBackbone =  load( [movieData.outputDirectory_ filesep ... 
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'I_neurite_orientation' filesep 'Channel_' num2str(iCh) ...
        filesep 'params.mat']); 
    pBackbone = pBackbone.p;
    BBScales = arrayfun(@(x) pBackbone(x).BBScale,1:length(pBackbone),'uniformoutput',0)'; 
%% Test for a manually chosen local patch size
    if exist([outDirC filesep 'LocalThreshPatchSizeTest' filesep ... 
          'manualPatchSizeSelect.mat'  ],'file')==2; 
      load([outDirC filesep 'LocalThreshPatchSizeTest' filesep ... 
          'manualPatchSizeSelect.mat']); 
      p.LocalThresholdPatchSize = patchSize;  
      display(['Manually Chosen Patch Size of ' num2str(p.LocalThresholdPatchSize) ' Will Be Used']); 
    end 
%% Main Function: 
    % (Note: Requires temporal information so input is per movie) 
    p.StartFrame = startFrame;
    p.EndFrame = endFrame;
    p.OutputDirectory = outDirC; 
    p.BBScales = BBScales; 
    
    if p.StartFrame == 1
        
    GCAReconstructVeilStem(listOfImages,backboneInfo,BBScales,p);
    else % for restart
        x{1} = veilStem;
        y{1} = paramsArchived; 
        GCAReconstructVeilStem(listOfImages,backboneInfo,BBScales,x,y,p); 
    end 
    
end % for iCh
end % function
