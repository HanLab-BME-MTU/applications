function movieData = shadeCorrectMovie(movieData,varargin)

% 
% movieData = shadeCorrectMovie(movieData)
% 
% movieData = shadeCorrectMovie(movieData,'OptionName',optionValue,...)
% 
%
% This function corrects the input movie for uneven illumination using
% "shade correction" images - Images taken of a blank area of a coverslip.
% If multiple shade correction images were taken, they are averaged. Also
% they can be spatially and/or median filtered to reduce noise and "hot
% pixels"
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
%       The integer indices of the channel(s) to perform shade correction
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels. (Unless in
%       batch mode, in which case this MUST be input)
%
%       ('ShadeImageChannels'-> Positive integer scalar or vector)
%       The integer indices of the channels which contain the shade
%       correction images. This must be the same size as ChannelIndex!
%       If not input, user will be asked to select, unless batch mode is
%       enabled, in which case an error will be generated.
%
%       ('MedianFilter' - True/False)
%       If true, the final (averaged) shade correction image will be median
%       filtered with a 3x3 neighborhood.
%       Optional. Default is true.
%
%       ('GaussFilterSigma' -> Positive scalar, >= 1.0)
%       This specifies the sigma (in pixels) of the gaussian filter to
%       apply to the final (averaged) shade correction image. If less than
%       one, no gaussian filtering is performed.
%       Optional. Default is no filtering.
%
%       ('Normalize' -> 0,1,2)
%       If set to 1, they will be divided by their own mean, resulting in
%       their mean being equal to one.        
%       If set to 0, the processed/averaged shade images will be divided by
%       their combined mean intensity. This results in their means being
%       nearer to one, minimizing rounding error in the final
%       shade-corrected image, while preserving their relative intensities.
%       If set to 0, no normalization will be done.
%       Default is 1.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
%
%   
% Output:
%
%   The shade-corrected images are saved as a new channel in the movie -
%   They are written as .tif files to the image directory of the input
%   movie.
% 
%   The location of the output images, and all parameters, are stored in
%   the movieData structure.
%
% Hunter Elliott
% 11/2009

%% ------ Parameters ------- %%

pString = 'sc_'; %The string to prepend before the shade-corrected image directory & channel name
saveName = 'shade_correction'; %File name for saving processed/avged shade images. Actual file name will have channel name appended.

%% ----------- Input ------------ %%

%Validate the input movieData
movieData = setupMovieData(movieData,'shadeCorrection');

[batchMode,iChannels,iShadeIm,medFilt,sigGFilt,doNorm] = parseInput(varargin);



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
if isempty(doNorm)
    doNorm = 1;
end

%Ask the user for the image channels to correct if not input
if isempty(iChannels)  
    if batchMode
        error('In batch mode, you must specify the channels to perform shade correction on!')
    else
        %As the user to select a channel.
        iChannels = selectMovieChannels(movieData,1,'Select the channels to shade correct:');        
    end    
end

nChan = length(iChannels); % Number of channels to shade correct

%Ask the user for the channels with the shade correction images, if not input
if isempty(iShadeIm)
    if batchMode
        error('In batch mode, you must specify the channels containing the shade correction images!')
    else
        %As the user to select channel(s)
        for i = 1:nChan
            iShadeIm(i) = selectMovieChannels(movieData,0,...
                ['Select the shade correction channel for image channel ' movieData.channelDirectory{iChannels(i)} ':']);
        end
    end
end

if length(iShadeIm) ~= nChan
    error('You must specify a shade correction image channel for each channel you are correcting!')
end

%Check that the shade-correction images are the same size as the images to
%correct
if ~all(all(movieData.imSize(:,iShadeIm) == movieData.imSize(:,iChannels)))
    error('Shade correction images must be the same size as the images to be corrected!')
end


%% ------------ Init ---------- %%


%Create gaussian filter, if needed
if sigGFilt >= 1
    gFilt = fspecial('gaussian',roundOddOrEven(sigGFilt*6,'odd',Inf),sigGFilt);        
end


disp('Starting shade correction...')

%Get image file names for needed channels
shadeImNames = getMovieImageFileNames(movieData,iShadeIm);
imNames = getMovieImageFileNames(movieData,iChannels);

%Go through the channels and check output directories/channels
iCorrected = zeros(1,nChan); %The index for the output channels
corrDir = cell(1,nChan); %The directory names for the shade-corrected images
for i = 1:nChan
    %Check if this channel has been shade-corrected before
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



%% ----------- Get and Process Shade Correction Images ------------- %%
%Loads, averages, filters  and normalizes the shade correction images

disp('Loading and processing correction image(s)...')

%Go through each requested channel and process the shade correction
shadeIm = cell(1,nChan);
for iChan = 1:nChan
    
    % ---- Average the shade images --- %
    nImages = movieData.nImages(iShadeIm(iChan));
    for iImage = 1:nImages    
        
        currIm = imread(shadeImNames{iChan}{iImage});
        
        if iImage == 1
           shadeIm{iChan} = zeros(size(currIm));
        end
        
        %Average the images together
        shadeIm{iChan} = shadeIm{iChan} + double(currIm) ./ nImages;                
               
        
    end
    
    %---Filter the averaged shade image---%
    
    %Median filter
    if medFilt                
        shadeIm{iChan} = medfilt2(shadeIm{iChan},'symmetric');
    end
    
    %Gaussian filter
    if sigGFilt >= 1
        shadeIm{iChan} = imfilter(shadeIm{iChan},gFilt,'replicate');                    
    end
        
    
end

%---Normalize the Shade Images---%
if doNorm == 1
    %Divide each image by its own mean
    shadeIm = cellfun(@(x)(x ./ mean(x(:))),shadeIm,'UniformOutput',false);            
elseif doNorm == 2
    %Calculate the combined mean of all the shade images
    combMean = mean(cellfun(@(x)(mean(x(:))),shadeIm));
    %Divide each image by this mean
    shadeIm = cellfun(@(x)(x ./ combMean),shadeIm,'UniformOutput',false);
end


%% -------------- Apply shade Correction ------------%%
%Applies the shade correction from above to each selected channel

disp('Applying shade correction to images...')

%Go through each image and apply the appropriate shade correction
for iChan = 1:nChan
    
    disp(['Correcting channel "' movieData.channelDirectory{iChannels(iChan)} ...
        '" using images from "' movieData.channelDirectory{iShadeIm(iChan)} '"']);     
    disp(['Resulting images will be stored in channel ' movieData.channelDirectory{iCorrected(iChan)}])
    
    for iImage = 1:movieData.nImages(iChannels(iChan))
    
        %Load the image to be corrected
        currIm = imread(imNames{iChan}{iImage});
        
        ogClass = class(currIm);
    
        %Correct it
        currIm = double(currIm) ./ shadeIm{iChan};
        %Cast to original class
        currIm = cast(currIm,ogClass);
        
        %Write it to disk
        iLastSep = max(regexp(imNames{iChan}{iImage},filesep));
        imwrite(currIm,[corrDir{iChan} filesep pString ... 
            imNames{iChan}{iImage}(iLastSep+1:end)]);
        
        
 
    end
end


%% ------------- Output ------- %%

disp('Saving results...')


%Save the averaged/filtered shade images
for i = 1:nChan
    processedShadeImage = shadeIm{i}; %#ok<NASGU> %Get this element of save array because the save function sucks.
    save([movieData.shadeCorrection.directory filesep saveName '_' movieData.channelDirectory{iShadeIm(i)} '.mat'],'processedShadeImage');
end

%Log the shade correction in the movieData structure and save it
movieData.shadeCorrection.parameters.medianFilter = medFilt;
movieData.shadeCorrection.parameters.sigmaGaussFilter = sigGFilt;
movieData.shadeCorrection.parameters.normalize = doNorm;
movieData.shadeCorrection.status = 1;
movieData.shadeCorrection.dateTime = datestr(now);
movieData.shadeCorrection.iShadeImages(iCorrected) = iShadeIm; %Indexes of the shade images used to correct
movieData.shadeCorrection.iFrom(iCorrected) = iChannels; %Indexes of the channels that were used to create the shade corrected images
movieData.shadeCorrection.fileNamePrefix = saveName;
movieData.shadeCorrection.channelNamePrefix = pString;

updateMovieData(movieData);

disp('Finished!')


end

function [batchMode,iChannels,iShadeIm,medFilt,sigGFilt,doNorm] = parseInput(argArray)


%Init output
batchMode = [];
iChannels = [];
iShadeIm = [];
medFilt = [];
sigGFilt = [];
doNorm = [];

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
           
       case 'ShadeImageChannels'
           iShadeIm = argArray{i+1};

       case 'MedianFilter'
           medFilt = argArray{i+1};
           
       case 'GaussFilterSigma'
           sigGFilt = argArray{i+1};
           
       case 'Normalize'
           doNorm = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end

end