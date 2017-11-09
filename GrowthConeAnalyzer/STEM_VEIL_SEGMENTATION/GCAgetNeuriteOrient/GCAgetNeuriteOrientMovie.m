function [ backboneInfo] = GCAgetNeuriteOrientMovie(movieData,varargin)
%% GCAgetNeuriteOrientMovie.m :  movieDataWrapper function for GCAgetNeuriteOrient
% STEP I in Veil/Stem estimation of GCA Segmenttation
%
% %%PERSONAL NOTE: was getNeuriteOrientationMovie.m until 20140529)%%
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
% for now check movieData separately. 
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
 defaultOutDir = [movieData.outputDirectory_ filesep...
     'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'I_neurite_orientation'];

 ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x)); 
 ip.addParameter('ChannelIndex',1); 
 ip.addParameter('ProcessIndex',0);
 ip.addParameter('StartFrame','auto'); 
 ip.addParameter('EndFrame','auto'); 

% Specific
 ip.addParameter('TSOverlays',true,@(x) islogical(x));
 ip.addParameter('screen2png',false); 
 
 ip.addParameter('BBScale',[5 6 7 8 9 10]);
 ip.addParameter('FilterOrderBB',4,@(x) ismember(x,[2,4]));
 
 ip.addParameter('MaxRadiusLargeScaleLink',10) ;
 ip.addParameter('MaxDistNoGeoTerm',3); 
 
 ip.addParameter('ThreshNMSResponse',25); 
 ip.addParameter('MinCCRidgeBeforeConnect',3);
 ip.addParameter('MinCCRidgeAfterConnect',5);

ip.addParameter('MinCCEntranceRidgeFirstTry',10);
ip.addParameter('MaxDistBorderFirstTry',10);
 
 
 
ip.parse(varargin{:});


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
params = ip.Results; 


%% Loop for each channel
for iCh = 1:nChan
    
    display(['Finding Neurite Orientation Channel ' num2str(iCh)]);
    %% Get Start and End Frames Based on Restart Choice
    
    % make final output dir where backboneInfo will be saved
    saveDir =  [ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh)];
    
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
   
    orientFile = [saveDir filesep 'backboneInfo.mat'];
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
        imDir = movieData.proceses_{ip.Results.ProcessIndex}.outfilePaths_;
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
        [backboneFrame,TSFigs] = GCAgetNeuriteOrient(img,ip.Results); % quick fix for the plots is to just make the frame number an input for not 20140812
       
        % make the directories for the figures if required. 
        for iFig = 1:length(TSFigs)
           if ~isdir([saveDir filesep TSFigs(iFig).name]); 
               mkdir([saveDir filesep TSFigs(iFig).name]); 
           end  
        end 
            type{1} = '.fig'; 
            type{2} = '.tif'; 
            
        if ~isempty(TSFigs)
            if ip.Results.screen2png
                
                arrayfun(@(x) helperScreen2png([saveDir filesep x.name filesep ...
                    num2str(iFrame,'%03d') '.png'],'figureHandle',x.h),TSFigs);
            else
                
                for iType = 1:numel(type)
                    
                    
                    arrayfun(@(x) saveas(x.h,...
                        [saveDir filesep x.name filesep num2str(iFrame,'%03d') type{iType}]),TSFigs);
                end
            end
        end 
            
        close all
         hashTag =  gcaArchiveGetGitHashTag;
         backboneFrame.hashTag = hashTag; % make sure to add the hash tag first so the structure is similar (or initiate in begin)
        backboneInfo(iFrame) = backboneFrame;
        save( [saveDir filesep 'backboneInfo.mat'],'backboneInfo');
        display(['Finished Finding Neurite Orientation for Frame ' num2str(iFrame)]);
        p(iFrame) = params; 
        save([saveDir filesep 'params.mat'],'p');  
    
    end % iFrame
    
end   % for iCh

end % function


