function movieData = bleedthroughCorrectMovie(movieData,varargin)
%BLEEDTHROUGHCORRECTMOVIE corrects for bleedthrough of other fluorophores in the input movie.
% 
% movieData = bleedthroughCorrectMovie(movieData)
% 
% movieData = bleedthroughCorrectMovie(movieData,'OptionName',optionValue,...)
% 
%
% This function corrects for bleedthrough of other fluorophores into a
% particular channel in the input movie. This is done using bleedthrough
% coefficients calculated from a "bleedthrough movie" where only one
% fluorophore is present, using processBleedthroughMovie.m. These
% coefficients are then used, in combination with images of the fluorophore
% channels which are bleeding into the image to be corrected, to remove the
% effects of bleedthrough.
% 
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
%       The integer index of the channel to perform bleedthrough correction
%       on. This index corresponds to the channel directory's location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels. (Unless in
%       batch mode, in which case this MUST be input)
%
%       ('BleedImageChannels'-> Positive integer scalar or vector)
%       The integer indices of the channels whose fluorophore is bleeding
%       through into the channel to be corrected. You will need to have a
%       bleedthrough coefficient for each of these channels.
%       If not input, user will be asked to select, unless batch mode is
%       enabled, in which case an error will be generated.
%
%       ('BleedCoefficients' -> Positive scalar or vector)
%       The bleedthrough coefficients for each bleed channel specified by
%       BleedImageChannels, in the same order. This should be the average
%       coefficient produced by calculateMovieBleedthrough.m
%       MUST BE INPUT!
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
%
%   
% Output:
%
%   The bleedthrough-corrected images are saved as a new channel in the movie -
%   They are written as .tif files to the image directory of the input
%   movie.
% 
%   The location of the output images, and all parameters, are stored in
%   the movieData structure.
%
% Hunter Elliott
% 11/2009

%% ------ Parameters ------- %%

pString = 'btc_'; %The string to prepend before the bleedthrough-corrected image directory & channel name

%% ----------- Input ------------ %%

%Validate the input movieData
movieData = setupMovieData(movieData,'bleedthroughCorrection');

[batchMode,iChannels,iBleedIm,bleedCoef] = parseInput(varargin);



% --- Defaults ---- %

if isempty(batchMode)
    batchMode = false;
end

%Ask the user for the image channel to correct if not input
if isempty(iChannels)  
    if batchMode
        error('In batch mode, you must specify the channel to perform bleedthrough correction on!')
    else
        %As the user to select a channel.
        iChannels = selectMovieChannels(movieData,0,'Select the channel to bleedthrough correct:');        
    end    
end

%Ask the user for the channels which bleed into the images to be corrected, if not input
if isempty(iBleedIm)
    if batchMode
        error('In batch mode, you must specify the channels containing the bleedthrough correction images!')
    else
        %As the user to select channel(s)        
        iBleedIm = selectMovieChannels(movieData,1,...
            ['Select the bleedthrough channels for image channel ' movieData.channelDirectory{iChannels(1)} ':']);        
    end
end


if isempty(bleedCoef)
    error('Bleedthrough coefficients MUST be input!')   
elseif length(iBleedIm) ~= length(bleedCoef)
    error('You must specify the same number of bleedthrough coefficients and bleedthrough images!')
end
nBleed = length(bleedCoef);
    


%% ------------ Init ---------- %%


disp('Starting bleedthrough correction...')

%Get image file names for needed channels
bleedImNames = getMovieImageFileNames(movieData,iBleedIm);
imNames = getMovieImageFileNames(movieData,iChannels(1));

nChan = length(iChannels);

%Go through the channels and check output directories/channels
iCorrected = zeros(1,nChan); %The index for the output channels
corrDir = cell(1,nChan); %The directory names for the bleedthrough-corrected images

%Check if this channel has been bleedthrough-corrected before
tmp = find(strcmp([pString movieData.channelDirectory{iChannels(1)}],movieData.channelDirectory),1);                    
if ~isempty(tmp)
    iCorrected(1) = tmp;
%If not, create a channel for it
else
    iCorrected(1) = length(movieData.channelDirectory)+1;
    movieData.channelDirectory{iCorrected(1)} = [pString movieData.channelDirectory{iChannels(1)}]; %Name the correct channel after the old one
    movieData.nImages(iCorrected(1)) = movieData.nImages(iChannels(1)); %The corrected images should have the same number
    movieData.imSize(:,iCorrected(1)) = movieData.imSize(:,iChannels(1)); %... and size
end

corrDir{1} = [movieData.imageDirectory filesep movieData.channelDirectory{iCorrected(1)}];

if ~exist(corrDir{1},'dir')
   mkdir(corrDir{1});
end

%this is to make compatible with input-output format of functions which
%process multiple channels
if length(iChannels) > 1
    iCorrected(2:end) = iChannels(2:end);
end

%% -------------- Apply bleedthrough correction ------------%%
%Applies the bleedthrough correction from above to each selected channel

disp('Applying bleedthrough correction to images...')

%Go through each image and apply the appropriate bleedthrough correction

disp(['Correcting channel "' movieData.channelDirectory{iChannels(1)} ...
    '" using images from "' movieData.channelDirectory{iBleedIm(1)} ...
    '" with coefficient ' num2str(bleedCoef(1)) ' and channel "' ...
    movieData.channelDirectory{iBleedIm(2)} '" with coefficient ' num2str(bleedCoef(2))]);     
disp(['Resulting images will be stored in channel ' movieData.channelDirectory{iCorrected(1)}])

for iImage = 1:movieData.nImages(iChannels(1))

    %Load the image to be corrected
    currIm = imread(imNames{1}{iImage});
    %check the bit-depth of the image
    ogClass = class(currIm);
    currIm = double(currIm);

    for iBleed = 1:nBleed

        %Load the bleed image
        currBleedIm = double(imread(bleedImNames{iBleed}{iImage}));

        %Subtract the bleedthrough from this channel
        currIm = currIm - (currBleedIm * bleedCoef(iBleed));
        
        %Remove negative values (these usually occur in the background)
        currIm(currIm < 0) = 0;

    end                


    %Cast to original class
    currIm = cast(currIm,ogClass);

    %Write it to disk
    iLastSep = max(regexp(imNames{1}{iImage},filesep));%find last file seperator
    imwrite(currIm,[corrDir{1} filesep pString ... 
        imNames{1}{iImage}(iLastSep+1:end)]);



end

%% ------------- Output ------- %%

disp('Saving results...')


%Log the bleedthrough correction in the movieData structure and save it
movieData.bleedthroughCorrection.status = 1;
movieData.bleedthroughCorrection.dateTime = datestr(now);
movieData.bleedthroughCorrection.iBleed(:,iCorrected) = iBleedIm; %Indexes of the bleedthrough channels used to correct
movieData.bleedthroughCorrection.bleedCoef = bleedCoef; %Coefficients used in correction
movieData.bleedthroughCorrection.iFrom(iCorrected) = iChannels; %Indexes of the channels that were corrected/input
movieData.bleedthroughCorrection.channelNamePrefix = pString;

updateMovieData(movieData);

disp('Finished!')


end

function [batchMode,iChannels,iBleedIm,bleedCoef] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
iBleedIm = [];
bleedCoef = [];

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
           
       case 'BleedImageChannels'
           iBleedIm = argArray{i+1};
           
       case 'BleedCoefficients'
           bleedCoef = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end

end