function movieData = backgroundSubtractMovie(movieData,varargin)
%BACKGROUNDSUBTRACTMOVIE corrects background by subtracting the average unmasked value
%
% movieData = backgroundSubtractMovie(movieData)
%
% movieData = backgroundSubtractMovie(movieData,'OptionName',optionValue)
%
% This function performs background subtraction on the movie described by
% the input movieData. This is accomplished by averaging the intensity in
% the areas covered by background masks, as created using
% createMovieBackgroundMasks.m
%
%
% Input:
%
%   movieData - The structure describing the movie, as created using
%   setupMovieData.m
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) 
%       The integer indices of the channel(s) to perform background subtraction
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels
%
%       ('BackgroundMaskIndex' -> Positive integer scalar or vector)
%       The integer index of the channels to use the background masks from.
%       If not input, the masks are used from the same channels you are
%       correcting.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
% Hunter Elliott, 11/2009
%

%%  --------- Parameters ------- %%

pString = 'bs_'; %The string to prepend before the background-subtracted image directory & channel name
saveName = 'background_subtraction_values'; %File name for saving subtracted values

%% ----------- Input ------------ %%

%Validate & initialize the movieData
movieData = setupMovieData(movieData,'backgroundSubtraction');

[batchMode,iChannels,iBackMasks] = parseInput(varargin);

%--- Default Values----%

if isempty(batchMode)
    batchMode = false;
end

if ~checkMovieBackgroundMasks(movieData);
    error('Must create background masks before running background subtraction!')
end
if isempty(iChannels)
    if batchMode
        error('In batch mode, you must specify the channels to background subtract!')
    else        
        iChannels = selectMovieChannels(movieData,1,'Select the channel(s) to background subtract:');
    end
end
if isempty(iBackMasks)
    iBackMasks = iChannels;
end
if ~checkMovieBackgroundMasks(movieData,iBackMasks)
    error('Invalid background masks for specified channels! Please create background masks prior to subtraction!')
end
nChan = length(iChannels);
if length(iBackMasks) ~= nChan
    error('Must specify a background mask channel for each corrected channel!')
end


%% ------- Init ----- %%


disp('Starting background subtraction...')


%Go through the channels and check output directories/channels
iCorrected = zeros(1,nChan); %The index for the output channels
corrDir = cell(1,nChan); %The directory names for the background-subtracted images
for i = 1:nChan
    %Check if there is a channel with this name
    tmp = find(strcmp([pString movieData.channelDirectory{iChannels(i)}],movieData.channelDirectory),1);                    
    if ~isempty(tmp)
        iCorrected(i) = tmp;
    %If not, create a channel for it
    else
        iCorrected(i) = length(movieData.channelDirectory)+1;
        movieData.channelDirectory{iCorrected(i)} = [pString movieData.channelDirectory{iChannels(i)}]; %Name the correct channel after the old one
        movieData.nImages(iCorrected(i)) = movieData.nImages(iChannels(i)); %The corrected images should have the same number
        movieData.imSize(:,iCorrected(i)) = movieData.imSize(:,iChannels(i)); %... and size
    end
        
    corrDir{i} = [movieData.imageDirectory filesep movieData.channelDirectory{iCorrected(i)}];
    
    if ~exist(corrDir{i},'dir')
       mkdir(corrDir{i});
    end
end

imNames = getMovieImageFileNames(movieData,iChannels);



%Get the background mask names
bakDir = cell(1,nChan);
bakNames = cell(1,nChan);
for i = 1:nChan
   bakDir{i} = [movieData.backgroundMasks.directory filesep ...
       movieData.backgroundMasks.channelDirectory{iBackMasks(i)}];
   bakNames{i} = imDir([bakDir{i}]);      
        
end



%% ---- Background Subtraction ---- %%
%Go through each image and subtract the average value behind the background
%mask from the image.


backgroundValues = cell(1,nChan);

for iChan = 1:nChan

    disp(['Correcting channel "' movieData.channelDirectory{iChannels(iChan)} ...
        '" using background masks from channel "' movieData.channelDirectory{iBackMasks(iChan)} '"']);            
    disp(['Resulting images will be stored in channel ' movieData.channelDirectory{iCorrected(iChan)}])
    
    
    nImages = movieData.nImages(iChannels(iChan));
    backgroundValues{iChan} = zeros(1,nImages);
    for iImage = 1:nImages
    
        currIm = imread(imNames{iChan}{iImage});
        ogClass = class(currIm); %Determine class so it is not altered later
        
        
        currBackMask = imread([bakDir{iChan} filesep bakNames{iChan}(iImage).name]);
        
        %Get average background intensity 
        backgroundValues{iChan}(iImage) = mean(double(currIm(currBackMask(:))));        
        currIm = double(currIm) - backgroundValues{iChan}(iImage);
        currIm(currIm < 0) = 0; %Clip negative values to zero
        currIm = cast(currIm,ogClass); %Return to original data type

        %Write to the new directory
        iLastSep = max(regexp(imNames{iChan}{iImage},filesep));
        imwrite(currIm,[corrDir{iChan} filesep pString ... 
            imNames{iChan}{iImage}(iLastSep+1:end)]);

    end
end



%% ----- Output / Finalization ----- %%

disp('Finishing, saving results...')

for i = 1:nChan
   subtractedValues = backgroundValues{i}; %#ok<NASGU>
   save([movieData.backgroundSubtraction.directory filesep saveName '_' ...
       movieData.channelDirectory{iCorrected(i)}],'subtractedValues');              
end



movieData.backgroundSubtraction.dateTime = datestr(now);
movieData.backgroundSubtraction.status = 1;
movieData.backgroundSubtraction.iBackgroundMasks(iCorrected) = iBackMasks; %Indices of channels background masks were from
movieData.backgroundSubtraction.iFrom(iCorrected) = iChannels; %Indices of the images that were background subtracted
movieData.backgroundSubtraction.fileNamePrefix = saveName;
movieData.backgroundSubtraction.channelNamePrefix = pString;

updateMovieData(movieData);


function [batchMode,iChannels,iBackMasks] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
iBackMasks = [];

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
           
       case 'BackgroundMaskIndex'
           iBackMasks = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end