function [ backboneInfo] = GCAgetNeuriteOrientMovie(movieData,varargin)
% GCAgetNeuriteOrientMovie.m :  movieWrapper function for GCAgetNeuriteOrient
% STEP I in Veil/Stem estimation of GCA package
%
% %%PERSONAL NOTE: was getNeuriteOrientationMovie.m until 20140529)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData (REQUIRED)  - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%  
%
%
%   
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
%  See GCAgetNeuriteOrient.m for details. 
%
%       
%
%
%
%
%       
%
%
% OUTPUT: (see main function GCAgetNeuriteOrient- for details):
%      
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
 defaultOutDir = [movieData.outputDirectory_ filesep...
     'Segmentation' filesep 'neurite_orientation_estimations'];

 ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x)); 
 ip.addParameter('ChannelIndex',1); 
 ip.addParameter('ProcessIndex',0);
 ip.addParameter('StartFrame','auto'); 
 ip.addParameter('EndFrame','auto'); 

% Specific
 ip.addParameter('TSOverlays',true,@(x) islogical(x));
 ip.addParameter('BBScale',[5 6 7 8 9 10]);
 ip.addParameter('FilterOrderBB',4,@(x) ismember(x,[2,4]));
 ip.addParameter('MaxRadiusLargeScaleLink',10) ;
 ip.addParameter('ThreshNMSResponse',25); 
 ip.addParameter('MinCCRidgeBeforeConnect',3);
 ip.addParameter('MinCCRidgeAfterConnect',5);

ip.addParameter('MinCCEntranceRidgeFirstTry',10);
ip.addParameter('MaxDistBorderFirstTry',10);
 
 
 
ip.parse(varargin{:});

% % quick fix is to just set the defaults again here. 
% if isempty(paramsInS); 
%     
% ip.Results.paramInS.'TSOverlays',true,@(x) islogical(x));
% 
% ip.addParamValue('BBScale',[5 6 7 8 9 10]);
% ip.addParamValue('FilterOrderBB',4,@(x) ismember(x, 2,4));
% 
% ip.addParamValue('MaxRadiusLargeScaleLink',10) ;
% 
% ip.addParamValue('ThreshNMSResponse',25);
% ip.addParamValue('MinCCRidgeBeforeConnect',3);
% ip.addParamValue('MinCCRidgeAfterConnect',5);
% 
% ip.addParamValue('MinCCEntranceRidgeFirstTry',10);
% ip.addParamValue('MaxDistBorderFirstTry',10);


% %% OLD 
% if nargin < 2
%     % Generic
%     ip.Results.OutputDirectory = [movieData.outputDirectory_ filesep 'GCSegmentation' filesep 'neurite_orientation_estimations'];
%     ip.Results.ChannelIndex = 1;
%     ip.R.ProcessIndex = 0; % use raw images
%     paramsIn.startFrame = 'auto';   % 'auto'
%     paramsIn.endFrame = 'auto'; % 'auto' or number, default auto.    
%     paramsIn.plots = 1;
%     
%    % Specific
%     paramsIn.BBScale = 5:10; % PERSONAL NOTE: not sure if this is optimal ...
%     paramsIn.FilterOrderBB = 4;  
%     paramsIn.MaxRadiusLargeScaleLink = 10; % Defines the Search Radius for Linking Linear Pieces by the KD tree
%     paramsIn.ThreshNonMaxSuppResponse = 25; 
%     
% end
%%
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

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(ip.Results.ChannelIndex);

%% Loop for each channel
for iCh = 1:nChan
    
    display(['Finding Neurite Orientation Channel ' num2str(iCh)]);
    %% Get Start and End Frames Based on Restart Choice
    
    % make final output dir where backboneInfo will be saved
    saveDir =  [ip.Results.OutputDirectory filesep 'Neurite_Backbone_Seed_Channel_' num2str(iCh)];
    
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % If file exists
    if  exist(orientFile,'file')==2;
        load(orientFile) % load the file
        display('Loading Previously Run Orientation Estimations');
        if strcmpi(ip.Results.StartFrame,'auto')
            startFrame = numel(backboneInfo)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = ip.Results.StartFrame; % use user input
            display(['Manual Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if previous file doesn't exist... need to start at 1  
        startFrame = 1;
        display('No Neurite Orientation Folder Found: Creating and Starting at Frame 1');
    end % exist(orientFile,'file') == 2
    
    if strcmpi(ip.Results.EndFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Finding Neurite Orientations From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = ip.Results.EndFrame;
        display(['Manual End: Finding Neurite Orientations From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    
    %%
    
    
    
    % get the list of image filenames
    if ip.Results.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
    else
        imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
        
    %% Start Movie Loop %%%%
    for iFrame = startFrame:endFrame
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        % note here we just pass in ALL of the parameters even the
        % non-specific for now - if the unmatched option is set to true in 
        % the main function there is no problem with this. 
        [backboneFrame,TSFigs] = GCAgetNeuriteOrientFixInput(img,ip.Results); % quick fix for the plots is to just make the frame number an input for not 20140812
       
        % make the directories for the figures if required. 
        for iFig = 1:length(TSFigs)
           if ~isdir([saveDir filesep TSFigs(iFig).name]); 
               mkdir([saveDir filesep TSFigs(iFig).name]); 
           end  
        end 
            
            
        if ~isempty(TSFigs)
            arrayfun(@(x) saveas(x.h,...
                [saveDir filesep x.name filesep num2str(iFrame,'%03d') '.tif']),TSFigs);               
        end 
            
        
        hashTag =  gcaArchiveGetGitHashTag;
        backboneFrame.hashTag = hashTag; % make sure to add the hash tag first so the structure is similar (or initiate in begin)
        backboneInfo(iFrame) = backboneFrame;
        save( [saveDir filesep 'backboneInfo.mat'],'backboneInfo');
       
    end % iFrame
    
end   % for iCh

end % function


