function movieData = darkCurrentCorrectMovie(movieData,varargin)

% movieData = darkCurrentCorrectMovie(movieData)
% 
% movieData = darkCurrentCorrectMovie(movieData,'OptionName',optionValue,...)
% 
% This function performs dark-current correction on the input movie. This
% is accomplished by subtracting a "dark-current image" from each image in
% the channels to be corrected. The dark current image is an image taken
% with no light incident on the camera. If multiple dark-current images are
% taken, they will be averaged together. This is highly recommended.
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
%       The integer indices of the channel(s) to perform dark-current
%       correction on. This index corresponds to the channel directories
%       location in the cell array movieData.channelDirectory. If not
%       input, the user will be asked to select from the available
%       channels. (Unless in batch mode, in which case this MUST be input)
%
%       ('DarkImageChannels'-> Positive integer scalar or vector)
%       The integer indices of the channels which contain the dark-current
%       images. This must be the same size as ChannelIndex!
%       If not input, user will be asked to select, unless batch mode is
%       enabled, in which case an error will be generated.
%
%       ('MedianFilter' - True/False)
%       If true, the final (averaged) dark correction image will be median
%       filtered with a 3x3 neighborhood.
%       Optional. Default is true.
%
%       ('GaussFilterSigma' -> Positive scalar, >= 1.0)
%       This specifies the sigma (in pixels) of the gaussian filter to
%       apply to the final (averaged) dark correction image. If less than
%       one, no gaussian filtering is performed.
%       Optional. Default is no filtering.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
%
%   
% Output:
%
%   The dark-corrected images are saved as a new channel in the movie -
%   They are written as .tif files to the image directory of the input
%   movie.
% 
%   The location of the output images, and all parameters, are stored in
%   the movieData structure.
%
% Hunter Elliott
% 11/2009



%% ------ Parameters ------- %%

pString = 'dcc_'; %The string to prepend before the dark-corrected image directory & channel name
saveName = 'dark_current_correction'; %File name for saving processed/avged dark current images. Actual file name will have channel name appended.

%% ----------- Input ------------ %%

%Validate the input movieData
movieData = setupMovieData(movieData,'darkCurrentCorrection');

[batchMode,iChannels,iDarkIm,medFilt,sigGFilt] = parseInput(varargin);



% --- Defaults ---- %

if isempty(batchMode)
    batchMode = false;
end
if isempty(medFilt)
    medFilt = true;
end
if isempty(sigGFilt)
    sigGFilt = 0;
end

%Ask the user for the image channels to correct if not input
if isempty(iChannels)  
    if batchMode
        error('In batch mode, you must specify the channels to perform dark current correction on!')
    else
        %As the user to select a channel.
        iChannels = selectMovieChannels(movieData,1,'Select the channels to dark current correct:');        
    end    
end

nChan = length(iChannels); % Number of channels to dark current correct

%Ask the user for the channels with the dark current correction images, if not input
if isempty(iDarkIm)
    if batchMode
        error('In batch mode, you must specify the channels containing the dark current correction images!')
    else
        %As the user to select channel(s)
        for i = 1:nChan
            iDarkIm(i) = selectMovieChannels(movieData,0,...
                ['Select the dark current correction channel for image channel ' movieData.channelDirectory{iChannels(i)} ':']);
        end
    end
end

if length(iDarkIm) ~= nChan
    error('You must specify a dark current correction image channel for each channel you are correcting!')
end

%Check that the dark current-correction images are the same size as the images to
%correct
if ~all(all(movieData.imSize(:,iDarkIm) == movieData.imSize(:,iChannels)))
    error('Dark current correction images must be the same size as the images to be corrected!')
end


%% ------------ Init ---------- %%


%Create gaussian filter, if needed
if sigGFilt >= 1
    gFilt = fspecial('gaussian',roundOddOrEven(sigGFilt*6,'odd',Inf),sigGFilt);        
end


disp('Starting dark current correction...')

%Get image file names for needed channels
darkImNames = getMovieImageFileNames(movieData,iDarkIm);
imNames = getMovieImageFileNames(movieData,iChannels);

%Go through the channels and check output directories/channels
iCorrected = zeros(1,nChan); %The index for the output channels
corrDir = cell(1,nChan); %The directory names for the dark current-corrected images
for i = 1:nChan
    %Check if this channel has been dark current-corrected before
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



%% ----------- Get and Process Dark Current Correction Images ------------- %%
%Loads, averages, and filters the dark current correction images

disp('Loading and processing correction image(s)...')

%Go through each requested channel and process the dark current correction
darkIm = cell(1,nChan);
for iChan = 1:nChan
    
    % ---- Average the dark images --- %
    nImages = movieData.nImages(iDarkIm(iChan));
    for iImage = 1:nImages    
        
        currIm = imread(darkImNames{iChan}{iImage});
        
        if iImage == 1
           darkIm{iChan} = zeros(size(currIm));
        end
        
        %Average the images together
        darkIm{iChan} = darkIm{iChan} + double(currIm) ./ nImages;                
               
        
    end
    
    %---Filter the averaged dark image---%
    
    %Median filter
    if medFilt
        %Add a border to prevent distortion        
        darkIm{iChan} = medfilt2(darkIm{iChan},'symmetric'); %Uses default 3x3 neighborhood        
    end
    
    %Gaussian filter
    if sigGFilt >= 1
        darkIm{iChan} = imfilter(darkIm{iChan},gFilt,'replicate');                    
    end
        
    
end



%% -------------- Apply dark Correction ------------%%
%Applies the dark correction from above to each selected channel

disp('Applying dark current correction to images...')

%Go through each image and apply the appropriate dark current correction
for iChan = 1:nChan
    
    disp(['Correcting channel "' movieData.channelDirectory{iChannels(iChan)} ...
        '" using images from "' movieData.channelDirectory{iDarkIm(iChan)} '"']);     
    disp(['Resulting images will be stored in channel ' movieData.channelDirectory{iCorrected(iChan)}])
    
    for iImage = 1:movieData.nImages(iChannels(iChan))
    
        %Load the image to be corrected
        currIm = imread(imNames{iChan}{iImage});
        
        ogClass = class(currIm);
    
        %Correct it
        currIm = double(currIm) - darkIm{iChan};
        %Cast to original class
        currIm = cast(currIm,ogClass);
        
        %Write it to disk
        iLastSep = max(regexp(imNames{iChan}{iImage},filesep));%find the last file seperator
        imwrite(currIm,[corrDir{iChan} filesep pString ... 
            imNames{iChan}{iImage}(iLastSep+1:end)]);
        
         
    end
end


%% ------------- Output ------- %%

disp('Saving results...')


%Save the averaged/filtered dark current images
for i = 1:nChan
    processedDarkImage = darkIm{i}; %#ok<NASGU> %Get this element of save array because the save function sucks.
    save([movieData.darkCurrentCorrection.directory filesep saveName '_' movieData.channelDirectory{iDarkIm(i)} '.mat'],'processedDarkImage');
end

%Log the dark current correction in the movieData structure and save it
movieData.darkCurrentCorrection.parameters.medianFilter = medFilt;
movieData.darkCurrentCorrection.parameters.sigmaGaussFilter = sigGFilt;
movieData.darkCurrentCorrection.status = 1;
movieData.darkCurrentCorrection.dateTime = datestr(now);
movieData.darkCurrentCorrection.iDarkImages(iCorrected) = iDarkIm; %Indexes of the dark current images used to correct
movieData.darkCurrentCorrection.iFrom(iCorrected) = iChannels; %Indexes of the channels that were used to create the dark current corrected images
movieData.darkCurrentCorrection.fileNamePrefix = saveName;
movieData.darkCurrentCorrection.channelNamePrefix = pString;

updateMovieData(movieData);

disp('Finished!')


end

function [batchMode,iChannels,iDarkIm,medFilt,sigGFilt] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
iDarkIm = [];
medFilt = [];
sigGFilt = [];

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
           
       case 'DarkImageChannels'
           iDarkIm = argArray{i+1};

       case 'MedianFilter'
           medFilt = argArray{i+1};
           
       case 'GaussFilterSigma'
           sigGFilt = argArray{i+1};
                      
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end

end






















