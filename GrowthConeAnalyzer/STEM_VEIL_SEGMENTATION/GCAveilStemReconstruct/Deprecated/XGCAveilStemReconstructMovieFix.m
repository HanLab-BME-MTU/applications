function [analInfo] = GCAveilStemReconstructMovie(movieData,paramsIn)
%%
% attempting to reconstruct the thin (ridge like ) and fat portions
% (amorphous) portions of the neurite using differential detection
% techniques.
% it assumes ridge-like features of the neurite body were already estimated
% using a large multi-scale steerable ridge detector
% (output of getNeuriteOrientationMovie.m stored in
% backboneInfo(iFrame).maxNMSLarge
%
% detector in checkFilopodia structures (ie thin ridge like structures) are
% (hopefully) eroded in this step and reconstructed with higher fidelity
% using this neurite body estimation as a template in reconstructFilopodiaMovie

% INPUT :
%
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analInfo structure to.
%       If not input, the analInfo will be saved in the same directory
%       as the movieData, in a sub-directory called "neurite_body_masks"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite_body_masks.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       subtracted to potentially be input). If not input, the
%       information will be calculated from the
%
%
%       ('MaskToErod' -> char )
%       'gradient' or 'local'
%
%       ('Disk Size' -> scalar)
%       Default = 6 for LifeAct: 4 for membrane markers
%       This parameter specifies the radius of the disk to be used for
%       erosion of thin objects from the input mask (ie the filopodia)
%       Larger values will erode out thicker structures. Note for the lifeAct
%       channels that have very strong filopodia signal very often these
%       gradients for crossing filopodia tend not to be well segmented so
%       disks use are typically slightly bigger than if using a membrane
%       marker where the filpodia often exhibit weak signal relative to
%       entire image and are those not segmented at all.













% MD (the movie Data info)
% process Idx (the index of the processed files you want to read in 0 is
% the raw images)
% chIdx = the channel index to process

% detection params.:
%                 .filterParams.filterOrder
%                 .filterParams.scales
%                 .localThresh = 1 if you want to use local thresholding 2
%                 if active  contours (need to change name)
%                 .saveViaBB = 1 if you want to save via backbone ( will in
%                 the future always use this option
%                 .localThresh.patchSize = patch size for thresholding (only if local thresholding)
%                 .BBScale scale of backbone estimation: should probably
%                 try to make a set of multiple finding the max response
% restart: flag in case crashed and need to start analysis from where left
% off, 1 if restart 0 if not
% analInfo: only required if restarting.

% OUTPUT: currently everything is set-up to read out from analInfo. It will
% add fields .bodyEst and .masks these can likely go into separate .mat
% structures  in separate folders eventually into masks for channel 1 so
% the protrusion retraction vectors can be directly fed into this.  right
% now though keep together
%

%% Check Input% NOTE TO SELF FIX INPUT

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'neurite_body_masks'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
    paramsIn.MaskToErod = 'local' ;
    paramsIn.DiskSizeForErod = 6;
    paramsIn.plots = 1;
    paramsIn.makeMovie = 0;
    % collect from 'LocalThreshPatchSize'
    patchSizeFilename = [movieData.outputDirectory_ filesep ...
        'neurite_body_masks' filesep 'LocalThreshPatchSizeTest' filesep 'manualPatchSizeSelect.mat'];
    if exist(patchSizeFilename,'file')==2
        
        load(patchSizeFilename);
        
        
        paramsIn.patchSize = str2double(patchSize);
    else % set patch size
        paramsIn.patchSize = 50; % Default
    end
    
    paramsIn.restart.startFrame = 'auto';   % 'auto'
    paramsIn.restart.endFrame = 'auto';
    
    
    % manual will require a startFrame
    % number if the user does not enter
    % one a warning will come up but it
    % will proceed with the first frame
end

p= paramsIn;

%% Initiate
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);

%% Start Reconstructing Body

for iCh = 1:nChan
    
    display(['Finding Veil/Stem For Channel ' num2str(iCh)]);
    
    if p.ProcessIndex == 0
        imgDir =  movieData.channels_(p.ChannelIndex).channelPath_; % currently mainly needed if you do the local thresholding/ otherwise just overlay
    else
        imgDir = movieData.proccesses_(p.ProcessIndex).outFilePaths_{p.ChannelIndex};
    end
    
    outDirChan = [p.OutputDirectory filesep 'Neurite_Body_Masks_Channel_' num2str(iCh)];
    if ~isdir(outDirChan)
        mkdir(outDirChan);
    end
    
    save([outDirChan filesep 'params.mat'],'p') ;
    
    % collect images and initiate
    [listOfImages] = searchFiles('.tif',[],imgDir ,0);
    
    
    % will have an input such that you can choose the segmentation process
    if strcmpi(p.MaskToErod, 'gradient')
        % load the previously run gradient based masks (note this is currently
        % set up for me where I know I ran the gradient based masks)
        %          idxThresh= find(cellfun(@(x) sum(strcmpi(x.name_,'Thresholding')),movieData.processes_));
        %          cDir= movieData.processes_{idxThresh}.outFilePaths_{1};
        %          upOne = upDirectory(cDir,1);
        cDir  = [movieData.outputDirectory_ filesep 'masks' filesep 'Gradient_masks_for_channel_' num2str(p.ChannelIndex)]; % will eventualy make this pick one of the segprocesses
        listOfMasks  = searchFiles('.tif',[],cDir,0);
    end
    
    %% Define first and last frame for restart
    bodyEstFile = [movieData.outputDirectory_ filesep ...
        'neurite_body_masks' filesep 'Neurite_Body_Masks_Channel_1' filesep 'analInfoTestSave.mat'];
    
    % If file exists
    if  exist(bodyEstFile,'file')==2;
        load(bodyEstFile) % load the file
        display('Loading Previously Run Orientation Estimations');
        if strcmpi(p.restart.startFrame,'auto')
            startFrame = numel(backboneInfo)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = p.restart.startFrame; % use user input
            display(['Manual Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        
        startFrame = 1;
        
        display('No Veil/Stem Folder Found: Creating and Starting at Frame 1');
        
        
    end % exist(orientFile,'file') == 2
    
    if strcmpi(p.restart.endFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = p.restart.endFrame;
        display(['Manual End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    
    % Load backbone info for now assume that the one you want is always in
    % this file
    load([movieData.outputDirectory_ filesep ...
        'neurite_orientation_estimation_fixes' filesep  'Neurite_Backbone_Seed_Channel_' num2str(iCh)...
        filesep 'backboneInfoFix.mat']);
    
    
    %% Start Loop Over Frames
    for iFrame = startFrame:nFrames
        % Get the necessary input
        img=double(imread([movieData.getChannelPaths{iCh} filesep  movieData.getImageFileNames{iCh}{iFrame}]));
        backboneInfoC= backboneInfo(iFrame);
        
        %% Input mask for erosion - INPUT you favorite mask making scheme here
        % or read in from the segmentation package files.
        if strcmpi(p.MaskToErod,'gradient')
            maskForErod = logical(imread([listOfMasks{iFrame,2} filesep listOfMasks{iFrame,1}]));
        else
            % we found albeit not perfect local thresholding works fairly
            % well in most cases
            [~,maskForErod] = gcaThresholdOtsu_local(img, p.patchSize,3);
        end
        
        %% Run the Veil/Stem Reconstruct Algorithm % shit I keep on forgetting you do need some temporal info here... 
        % can't have so no a wrapper unless you use the temporal info
        % later...
        [analInfoC,hTrouble] =  GCAveilStemReconstruct(img,maskForErod,backboneInfoC,p,iFrame);
  
        % if the trouble shoot plots were made save them 
        if ~isempty(hTrouble); 
              troubleDir  =  [ outDirChan filesep 'TroubleshootVeilStem']; 
              if ~isdir(troubleDir) 
                  mkdir(troubleDir)
              end 
              saveas(hTrouble,[troubleDir filesep num2str(iFrame,'%03d') '.png']); 
              saveas(hTrouble,[troubleDir filesep num2str(iFrame,'%03d') '.fig']); 
        end 
        
        % document git hash
        hashTag = gcaArchiveGetGitHashTag;
        analInfoC.hashTag = hashTag;
        analInfo(iFrame) = analInfoC;
        clear analInfoC 
        close(hTrouble); 
    end % for iFrame
    
   
    
end % for iCh
end % the END


