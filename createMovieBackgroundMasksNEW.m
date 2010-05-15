function movieData = createMovieBackgroundMasksNEW(movieData,paramsIn)
%CREATEMOVIEBACKGROUNDMASKS creates background masks by growing the foreground masks
%                                               
% movieData = createMovieBackgroundMasks(movieData);                                              
%
% movieData = createMovieBackgroundMasks(movieData,paramsIn)
% 
% This function uses the (already created) image masks to generate
% background masks. This is accomplished by "growing" the masks (dilation).
% 
% 
% Input:
%   
%   movieData - Structure describing movie created with setupMovieData.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       masks to. Masks for different channels will be saved as
%       sub-directories of this directory.
%       If not input, the masks will be saved to the same directory as the
%       movieData, in a sub-directory called "BackgroundMasks"
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       The integer index of the channels to create background masks from.
%       These channels must have already been segmented. 
%       Optional. If not input, background masks are created for all
%       channels which have foreground masks.
%
%       ('GrowthRadius' - positive integer scalar)
%       The radius (in pixels) to grow the foreground masks to
%       produce the background masks.
%       Optional. Default is 20 pixels.
% 
%       'BatchMode' - If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - the updated movieData structure with the background mask
%   creation logged in it.
%
%
% Additionally, the masks are written to the movie's analysis directory in
% a sub-folder called "backgroundMasks"
%
%
% Hunter Elliott, 11/2009
%
%% ----- Parameters ----- %%

pString = 'bkgrnd_'; %Prefix for saving masks to file

%% ------------ Input ----------- %%

if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end
if nargin < 2
    paramsIn = [];
end

%Make sure the move has been segmented
iSegProc = find(cellfun(@(x)(isa(x,'SegmentationProcess') && ~ isa(x,'BackgroundMaskProcess')),movieData.processes_),1);   
%WHAT IF THERE ARE MULTIPLE SEGMENTATION PROCESSES!!!! TEMP TEMP
if isempty(iSegProc) 
    error('Must create foreground masks before creating background masks!')
else
   %Check which channels have foreground masks 
   hasMasks = cellfun(@(x)(~isempty(x)),movieData.processes_{iSegProc}.maskPaths_);
end


%Get the indices of any previous background mask processes from this function                                                                              
iProc = find(cellfun(@(x)(isa(x,'BackgroundMaskProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(BackgroundMaskProcess(movieData));                                                                                                 
end

nChan = length(movieData.channelPath_);


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


%----Param Check-----%

if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

%Make sure all the selected channels have foreground masks.
if any(~hasMasks(p.ChannelIndex))
    warning('Cannot create background masks for some channels, because they do not have foreground masks! \n Please segment these channels before creating background masks!')
    p.ChannelIndex = p.ChannelIndex(hasMasks);
end



if p.GrowthRadius ~= round(p.GrowthRadius) || p.GrowthRadius < 1
    error('Input variable p.GrowthRadius must be a positive integer!')
end
   
nChanBack = numel(p.ChannelIndex);

%Set up the mask directories as sub-directories of the output directory
for j = 1:nChanBack;
    
    %Create string for current directory
    currDir = [p.OutputDirectory filesep 'BackgroundMasks_channel_' num2str(p.ChannelIndex(j))];    
    %Save this in the process object
    movieData.processes_{iProc}.setMaskPath(p.ChannelIndex(j),currDir);
   
    %Check/create directory
    if ~exist(currDir,'dir')
       mkdir(currDir)
    end        
end

%% ------------ Init --------------%%


growDisk = strel('disk',p.GrowthRadius);


%% ------- Background mask creation -------------%%

disp('Starting background mask creation...')


maskFileNames = movieData.processes_{iSegProc}.getMaskFileNames(p.ChannelIndex);

nMasks = movieData.nFrames_;
nMaskTot = nMasks * nChanBack;

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, thresholding channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end    

for iChan = 1:nChan
           
    disp(['Creating background masks for channel ' num2str(p.ChannelIndex(iChan)) '...']);

    currBkgrndMaskDir = movieData.processes_{iProc}.maskPaths_{p.ChannelIndex(iChan)};        
    currMaskDir = movieData.processes_{iSegProc}.maskPaths_{p.ChannelIndex(iChan)};        
    
    if ~p.BatchMode        
        waitbar((iChan-1)*nMasks / nMaskTot,wtBar,['Please wait, thresholding channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end        

    
    for iMask = 1:nMasks
        
        %Load the current foreground mask
        currMask = imread([currMaskDir filesep maskFileNames{iChan}{iMask}]);
        
        %Grow and invert this mask to create the background mask
        backgroundMask = ~imdilate(currMask,growDisk);
        
        %Write it to file        
        imwrite(backgroundMask,[currBkgrndMaskDir filesep ...
            pString maskFileNames{iChan}{iMask}]);

        if ~p.BatchMode && mod(iMask,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iMask + (iChan-1)*nMasks) / nMaskTot,wtBar)
        end                
        
    end
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

%% ------ Log processing in moviedata and save ---- %%

movieData.processes_{iProc}.setSuccess(true);
movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk


disp('Finished creating background masks!')




