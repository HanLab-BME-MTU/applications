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
%  GCAgetNeuriteOrient specific functions
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
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'neurite_orientation_estimations'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
    paramsIn.startFrame = 1;   %default = 1 
    paramsIn.startChoice =  'manual' ; % auto, manual % auto will look for analInfo and redo the last frame 
                                     % manual will require a startFrame
                                     % number if the user does not enter
                                     % one a warning will come up but it
                                     % will proceed with the first frame
    paramsIn.plots = 1;
    
    
    % Specific
    paramsIn.BBScale = 5:10; % PERSONAL NOTE: not sure if this is optimal ...
    paramsIn.FilterOrderBB = 4;
     
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
 orientFile = [movieData.outputDirectory_ filesep 'neurite_orientation_estimations' filesep 'Neurite_Backbone_Seed_Channel_1' filesep 'backboneInfo.mat'];
          
 if  exist(orientFile,'file')==2;
     load(orientFile)
     display('Loading Previously Run Orieintation Estimations');
     if strcmpi(p.startChoice,'auto')
         startFrame = numel(backboneInfo)-1;
         display(['Auto Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
     else
         startFrame = p.startFrame; % use user input
         if startFrame == 1 
             mkClrDir(p.OutputDirectory)
         end 

         display(['Manual Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
     end
 else
       mkClrDir(p.OutputDirectory)
     display('Auto Start: Creating Neurite Backbone Orientation Estimation Starting at Frame 1')
         
     startFrame = 1;
 end


%%




for iCh = 1:nChan
    display(['Finding Neurite Orientation Channel ' num2str(iCh)]); 
   
    % make final output dir where backboneInfo will be saved 
    saveDir =  [p.OutputDirectory filesep 'Neurite_Backbone_Seed_Channel_' num2str(iCh)]; 
    mkClrDir(saveDir)
    
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
    
    % if restart == 1
    %    % startFrame = numel(analInfo)-1;
    %
    % else
    %     startFrame = 1;
    % end
    
    %

    % initiate backbone structure 
   % backbone = struct([]); % can't initiate structure here - the fields
   % don't match when you add to it- figure out how to better initiate (
   % might need all the fields) 
    
    
    %%%% Start Frame Loop %%%%
    for iFrame = startFrame:nFrames
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        [backboneFrame] = GCAgetNeuriteOrient(img,p,iFrame); % quick fix for the plots is to just make the frame number an input for not 20140812
        backboneInfo(iFrame) = backboneFrame; 
         
        
    end % iFrame
    save( [paramsIn.OutputDirectory filesep 'backboneInfo.mat'],'backboneInfo'); 
    % choose correct orientation
    % collect all the indexes dilate this mask...
    % get CC
    
    
    
    
    %     idxRidge = find(cleanedRidge==1);
    %    [yRidge,xRidge] = ind2sub(size(imgCrop),idxRidge);
    %    idxBound = find(boundaryMask == 1)   ;
    %    [yBound,xBound] = ind2sub(size(imCrop),idxBound);
    %     [pixelListBoundary] =  KDTreeBallQuery( [xRidge,yRidge],[xBound,yBound], 3);
    
end   % for iCh
end
% function [coords] = getEndpoints(pixIdx,size)
% 
% endpoints = zeros(size);
% endpoints(pixIdx)=1;
% sumKernel = [1 1 1];
% % find endpoints of the floating candidates to attach (note in the
% % future might want to prune so that the closest end to the
% % body is the only one to be re-attatched: this will avoid double connections)
% endpoints = double((endpoints.* (conv2(sumKernel, sumKernel', padarrayXT(endpoints, [1 1]), 'valid')-1))==1);
% [ye,xe] = find(endpoints~=0);
% coords(:,1) = xe;
% coords(:,2) = ye;
% end
