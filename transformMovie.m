function movieData = transformMovie(movieData,varargin)

% movieData = transformMovie(movieData)
% 
% movieData = transformMovie(movieData,'OptionName',optionValue)
% 
% This function performs a spatial transformation on the selected channels
% of the input movie and writes the transformed images to a new channel in
% the movie. 
% This can be done using a pre-created transform that has been saved to
% file, or a new one can be calculated from alignment images (e.g. of a
% grid micrometer or multi-spectral beads)
% 
% 
% Input:
% 
%   movieData - movie information structure, as created by setupMovieData.m
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer scalar or vector)
%       The integer indices of the channel(s) to perform transformation
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels. (Unless in
%       batch mode, in which case this MUST be input)
%       If multiple channels are specified, the same transformation will be
%       applied to all channels.
%       
%
%       ('TransformFileName' -> Character string)
%       The FULL filename of the .mat file containing the transform to
%       apply to the images. The transform should be of the format used by
%       imtransform.m 
%       If not input, the user will be asked to locate a file containing a
%       transform, UNLESS batchmode is enabled, in which case an error will
%       be generated.
% 
%
%       ('TransformMasks' -> True/False)
%       If true, the masks for a given channel will also be transformed,
%       and saved to a mask directory for the output channel(s). If true,
%       the specified channels MUST have masks. Default is true.
%
%
%       ('MaskChannels' -> Positive integer scalar or vector)
%       The indices of the channel(s) to transform masks from if
%       TransformMasks is true. If not input, the masks will be used from
%       the channels which are transformed.
%
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed. Default is false.
%
% 
% Output:
%   
%   The transformed images are written to a new channel in the movie's
%   image directory.
%
%   The location of these images and all parameters used in the
%   transformation are stored in the movieData structure and it is saved.
%
% Hunter Elliott
% 11/2009

%% ------ Parameters ----%%

pString = 'x_'; %The string to prepend before the transformed image directory & channel name

%% ------- Input ------- %%


movieData = setupMovieData(movieData,'transformation');

[batchMode,iChannels,xFormFile,doMasks,iMaskChan] = parseInput(varargin);

%----Defaults----%

if isempty(batchMode)    
    batchMode = false;
end
if isempty(doMasks)
    doMasks = true;
end
if isempty(iChannels)   
    if ~batchMode
        iChannels = selectMovieChannels(movieData,1,'Select the channel(s) to spatially transform:');
    else
        error('In batch mode, you must specify the channels to transform!')
    end
end
if isempty(iMaskChan)
    iMaskChan = iChannels;
end
if doMasks && ~checkMovieMasks(movieData,iMaskChan)
    error('Cannot transform masks because the specified channels do not have valid masks! Disable the "TransformMasks" option, or create masks!!!')
end
if isempty(xFormFile)
    if ~batchMode
        [xFormFile xDir] = uigetfile('*.mat','Select a file containing the transformation:');
        if xFormFile == 0
            error('Must specify a transformation file to continue!')
        end
        xFormFile = [xDir xFormFile];
    else
       error('In batch mode, the transformation file name must be specified!') 
    end
end



%% ------- Init ------ %%

nChan = length(iChannels);


disp('Loading transformation...')
    

%Load the transform and check it

xForm = load(xFormFile);
tstField = fieldnames(xForm);

if length(tstField) > 1
    error('The specified .mat file should only contain one variable, which is the transform structure!')
else
   xForm = xForm.(tstField{1});     
end

if ~istransform(xForm)
    error('The specified .mat file does not contain a valid transformation! The transformation should use the same format used by imtransform.m!!!')
end


%Check / create the channel for the transformed images
iTransformed = zeros(1,nChan);
xDir = cell(1,nChan); %The directory names for the transformed images
for i = 1:nChan
    %Check if this channel already exists
    tmp = find(strcmp([pString movieData.channelDirectory{iChannels(i)}],movieData.channelDirectory),1);                    
    if ~isempty(tmp)
        iTransformed(i) = tmp;
    %If not, create a channel for it
    else
        iTransformed(i) = length(movieData.channelDirectory)+1;
        movieData.channelDirectory{iTransformed(i)} = [pString movieData.channelDirectory{iChannels(i)}]; %Name the correct channel after the old one
        movieData.nImages(iTransformed(i)) = movieData.nImages(iChannels(i)); %The corrected images should have the same number
        movieData.imSize(:,iTransformed(i)) = movieData.imSize(:,iChannels(i)); %... and size
    end
        
    xDir{i} = [movieData.imageDirectory filesep movieData.channelDirectory{iTransformed(i)}];
    
    if ~exist(xDir{i},'dir')
       mkdir(xDir{i});
    end
end


imNames = getMovieImageFileNames(movieData,iChannels);
if doMasks
    
    %Get the original mask file names    
    maskNames = getMovieMaskFileNames(movieData,iMaskChan);
    
    %Check/create the directories for the transformed masks
    xMaskDir = cell(1,nChan);
    for i = 1:nChan
        
        xMaskDir{i} = [movieData.masks.directory filesep movieData.channelDirectory{iTransformed(i)}];
    
        movieData.masks.channelDirectory{iTransformed} = [movieData.channelDirectory{iTransformed(i)}];
        
        if ~exist(xMaskDir{i},'dir')
            mkdir(xMaskDir{i})
        end                
    end    
end

%% ------- Spatial Transformation ------ %%
%Transform all images in requested channels and write them to a new
%channel.


disp('Transforming all images....')


for iChan = 1:nChan            
    
    nImages = movieData.nImages(iChannels(iChan));
    
    m = movieData.imSize(1,iChannels(iChan)); %Get the image size
    n = movieData.imSize(2,iChannels(iChan));
    
    disp(['Transforming channel ' movieData.channelDirectory{iChannels(iChan)}])
    disp(['Writing transformed images to channel ' movieData.channelDirectory{iTransformed}])
    
    if doMasks        
        disp(['Also transforming masks from ' movieData.channelDirectory{iMaskChan}])
    end
    
    for iImage = 1:nImages
        
        currIm = imread(imNames{iChan}{iImage});                                
        currIm = imtransform(currIm,xForm,'XData',[1 n],'YData',[1 m],'FillValues',1);
        iLastSep = max(regexp(imNames{iChan}{iImage},filesep));
        imwrite(currIm,[xDir{iChan} filesep pString ... 
            imNames{iChan}{iImage}(iLastSep+1:end)]);
        
        if doMasks            
            currMask = imread(maskNames{iChan}{iImage});            
            currMask = imtransform(currMask,xForm,'XData',[1 n],'YData',[1 m],'FillValues',0);                    
            
            iLastSep = max(regexp(maskNames{iChan}{iImage},filesep));
            imwrite(currMask,[xMaskDir{iChan} filesep pString ... 
            maskNames{iChan}{iImage}(iLastSep+1:end)]);
        
            
        end
    
    end    
end





%% ------ Output and Finalization ----- %%


movieData.transformation.dateTime = datestr(now);
movieData.transformation.status = 1;
movieData.transformation.iFrom(iTransformed) = iChannels;
movieData.transformation.fileName = xFormFile;
movieData.transformation.transformedMasks(iTransformed) = doMasks;
if doMasks
    movieData.transformation.iTransformedMasksFrom(iTransformed) = iMaskChan;
end

updateMovieData(movieData);


disp('Finished!')


function [batchMode,iChannels,xFormFile,doMasks,iMaskChan] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
xFormFile = [];
doMasks = [];
iMaskChan = [];

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
              
       case 'BatchMode'
           batchMode = argArray{i+1};
           
       case 'ChannelIndex'
           iChannels = argArray{i+1};
           
       case 'TransformFileName'
           xFormFile = argArray{i+1};
           
       case 'TransformMasks'
           doMasks = argArray{i+1};           
           
       case 'MaskChannels'
           iMaskChan = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end