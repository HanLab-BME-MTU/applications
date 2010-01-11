function movieData = createMovieBackgroundMasks(movieData,varargin)
%                                               
% movieData = createMovieBackgroundMasks(movieData);                                              
%
% movieData = createMovieBackgroundMasks(movieData,'OptionName',optionValue,...)
% 
% This function uses the (already created) image masks to generate
% background masks. This is accomplished by "growing" the masks (dilation).
% 
% 
% Input:
%   
%   movieData - Structure describing movie created with setupMovieData.m
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%       ('OptionName' -> possible values)
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       The integer index of the channels to create background masks from.
%       These channels must have already been segmented. 
%       Optional. If not input, background masks are created for all
%       channels which have foreground masks.
%
%       ('GrowthRadius' - positive integer scalar)
%       growRadius - The radius (in pixels) to grow the foreground masks to
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

movieData = setupMovieData(movieData,'backgroundMasks');

[iChannels,batchMode,growRadius] = parseInput(varargin);


%----Defaults-----%

if isempty(batchMode)
    batchMode =false;
end
if isempty(iChannels);
    iChannels = find(cellfun(@(x)(~isempty(x) | ischar(x)),movieData.masks.channelDirectory));%The ischar check covers the empty string case which adds flexibility
end
if isempty(growRadius)
    growRadius = 20;
elseif growRadius ~= round(growRadius) || growRadius < 1
    error('Input variable growRadius must be a positive integer!')
end
    

%% ------------ Init --------------%%

%Make sure the move has been segmented
if ~checkMovieMasks(movieData,iChannels)
    error('Must create foreground masks before creating background masks!')
end

nChan = length(iChannels);

bkMaskDir = movieData.backgroundMasks.directory;

growDisk = strel('disk',growRadius);

%% ------- Background mask creation -------------%%

disp('Starting background mask creation...')


maskFileNames = getMovieMaskFileNames(movieData,iChannels);


for iChan = 1:nChan
    
    
    nMasks = movieData.nImages(iChannels(iChan));        
    
    %Check/create directory for this channel
    currBkgrndMaskDir = [bkMaskDir filesep ...
        movieData.masks.channelDirectory{iChannels(iChan)}];    
    
    disp(['Creating background masks for channel ' movieData.masks.channelDirectory{iChannels(iChan)}]);
    
    
    if ~exist(currBkgrndMaskDir,'dir')
        mkdir(currBkgrndMaskDir);
    end    
    
    for iMask = 1:nMasks
        
        %Load the current foreground mask
        currMask = imread(maskFileNames{iChan}{iMask});
        
        %Grow this mask to create the background mask
        backgroundMask = ~imdilate(currMask,growDisk);
        
        %Write it to file
        iLastSep = max(regexp(maskFileNames{iChan}{iMask},filesep));
        imwrite(backgroundMask,[currBkgrndMaskDir filesep ...
            pString maskFileNames{iChan}{iMask}(iLastSep+1:end)]);
                

    end
end

disp('Finished!')

%% ------ Log processing in moviedata and save ---- %%

movieData.backgroundMasks.channelDirectory(iChannels) = ...
                         movieData.masks.channelDirectory(iChannels); 
movieData.backgroundMasks.status = 1;
movieData.backgroundMasks.dateTime = datestr(now);
movieData.backgroundMasks.iFrom(iChannels) = iChannels; %From and to channel indices are the same
movieData.backgroundMasks.parameters.growthRadius = growRadius;

updateMovieData(movieData);




function [iChannels,batchMode,growRadius] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
growRadius = [];


if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
       
       case 'ChannelIndex'
           
           iChannels = argArray{i+1};
              
       case 'BatchMode'
           batchMode = argArray{i+1};

       case 'GrowthRadius'
           growRadius = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end

