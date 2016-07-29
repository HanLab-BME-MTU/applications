function [ filoInfo] = GCAfitFilopodiaMovie(movieData,varargin)
%fitLinescansMovie: performs automated fitting of filopodia detections
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
%  SPECIFIC INPUT
%  See main function GCAFitFilopodia.m for details. 
%
% OUTPUT: (see main function GCAFitFilopodia- for details):
%       filoBranch is an Rx1 structure array where R = the number of frames
%       filoBranch holds the field .filoInfo a Rx1 structure where R is the number of
%       filopodia detected in that frame
%       .filoInfo has many subfields that store information regarding the
%       filopodia's coordinates,orientation,grouping(for branching) etc.
%       filoBranch is first created by GCAReconstructFilopodia. 
%       GCAFitFilopodia will add addtional fields to .filoInfo
%       corresponding to the the filopodia fitting. 
%% ----- Input ------ %%
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

%%Input check
ip = inputParser;

ip.CaseSensitive = false;



% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits_new'];


defaultInDir = [movieData.outputDirectory_ filesep ... 
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VI_filopodiaBranch_reconstruction']; 

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x)); 
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));
ip.addParameter('TSMovie',false,@(x) islogical(x)); 

ip.addParameter('InternalFiloOn',3,@(x) isscalar(x));
% ip.addParameter('NumPixForFitBack',10,@(x) isscalar(x));
ip.addParameter('ValuesForFit','Intensity',@(x) ischar(x)); % maybe remove 
ip.addParameter('PSFSigma',0.43,@(x) isnumeric(x)) ; %% NOTE CHANGE THIS TO BE READ IN FROM MD. 
ip.addParameter('fitLengthInPix',10,@(x) isscalar(x)); 

ip.parse(varargin{:});
params = ip.Results;
 %% Initiate
nFrames = movieData.nFrames_;
nChan = numel(params.ChannelIndex);
imSize = movieData.imSize_;
channels = params.ChannelIndex; 

%% Start Wrapper
for iCh = 1:nChan
    
    
    display(['Fitting Filopodia for Channel ' num2str(channels(iCh))]);
    

    
    %% Get Start and End Frames Based on Restart Choice
    
    
    
    % make final output dir where backboneInfo will be saved
    outDirC =  [ip.Results.OutputDirectory filesep 'Channel_' num2str(channels(iCh))];
    inDirC = [ip.Results.InputDirectory filesep 'Channel_' num2str(channels(iCh))]; 
    
    if ~isdir(outDirC)
        mkdir(outDirC);
    end
    
%     if params.ProcessIndex == 0
%         imgDir =  movieData.channels_(channels(iCh)).channelPath_; % currently mainly needed if you do the local thresholding/ otherwise just overlay
%     else
%         imgDir = movieData.proccesses_(params.ProcessIndex).outFilePaths_{params.ChannelIndex};
%     end
    
    % collect images and initiate
    %[listOfImages] = searchFiles('.tif',[],imgDir ,0);
    
   
%% Get Restart Information 
% Check for a fit file and load if exists 
    fitFile = [outDirC filesep 'filoBranch.mat'];

    
    % If fitFile already exists load it
    if  exist(fitFile,'file')==2;
        load(fitFile) % load the file
        display('Loading Previously Run Fits');
        if strcmpi(ip.Results.StartFrame,'auto')
             startFrame = find(arrayfun(@(x) ~isfield(filoBranch(x).filoInfo, 'Ext_exitFlag')...
                ,1:length(filoBranch)),1,'first')-1;
          
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting FiloBranch Fitting at Frame ' num2str(startFrame)]);
        else
            startFrame = ip.Results.StartFrame; % use user input
            display(['Manual Start: Starting FiloBranch Fitting at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist- force a start at 1 
        
        startFrame = 1;
        
        display('No Filopodia Fit Folder Found: Creating and Starting at Frame 1');
        reconstructFile = [inDirC filesep 'filoBranch.mat']; 
        load(reconstructFile); 
        display('Loading Reconstruction File'); 
        
    end % exist fitFile
    
    % if restarting add saved parameters 
     if startFrame ~= 1 
         load([outDirC filesep 'params.mat']); 
     end 
    
    
    if strcmpi(ip.Results.EndFrame,'auto');
        endFrame = nFrames-1;
        display(['Auto End: Fitting Filopodia/Branches From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = ip.Results.EndFrame;
        display(['Manual End: Fitting Filopodia/Branches From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
   
    %% Start Loop Over Movie
    % GET FRAME INFORMATION - this function wraps per frame
    for iFrame = startFrame:endFrame
        % get the filoInfo for the current frame
        filoInfo = filoBranch(iFrame).filoInfo;
        imgPath = [movieData.getChannelPaths{(channels(iCh))} filesep movieData.getImageFileNames{(channels(iCh))}{iFrame}];
        img = double(imread(imgPath));
        % make a specific output directory for the plotting for each frame
%         pSpecific = p;
%         pSpecific.sigma = movieData.channels_.psfSigma_;
%         if isempty(pSpecific.sigma)
%             display(['Using sigma 0.43']);
%             pSpecific.sigma = 0.43;
%         end
%         if pSpecific.SavePlots == 1
        if ip.Results.TSOverlays == 1
             params.OutputDirectory = [outDirC filesep 'Linescans' filesep 'Frame ' num2str(iFrame,'%03d')];
             
             
             mkClrDir(params.OutputDirectory)
        end
        
        [filoInfo] = GCAfitFilopodiaParamFreeInt(filoInfo,img,params) ;
        
        
        
        
        
        % rewrite the filoInfo with the extra filo Info fields.
        filoBranch(iFrame).filoInfo = filoInfo;
        display(['Finished Fitting Filopodia for  Channel ' num2str((channels(iCh))) 'Frame ' num2str(iFrame)]);
        filoBranch(iFrame).reconstructInfo.createTimeFiloFit = clock;
        %hashTag = gcaArchiveGetGitHashTag;
        %filoBranch(iFrame).reconstructInfo.hashTagFiloFit = hashTag;
        p(iFrame) = params; 
        save([outDirC filesep 'filoBranch.mat'],'filoBranch','-v7.3')
        save([outDirC filesep 'params.mat'],'p'); 
        
    end % for iFrame
end % for iCh
