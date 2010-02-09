function movieData = ratioMovie(movieData,varargin)
%RATIOMOVIE creates a new ratio channel by dividing one movie channel by another
% 
% movieData = ratioMovie(movieData)
% 
% movieData = ratioMovie(movieData,'OptionName',optionValue)
%
% This function divides the images in one channel by those in another and
% stores the resulting image in a new channel. Rounding error is minimized
% by having the resulting ratios fill the entire range of the input image
% bit-depth.
% 
% 
% Input:
% 
%   movieData - The structure describing the movie, as created using
%   setupMovieData.m
% 
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer vector, 1x2)
%       The index of the channels to use as the numerator (1st element) and
%       denominator (2nd element). If not input, the user will be asked. 
%
%       ('ApplyMasks' -> True/False)
%       If true, pixels which are not masked in EITHER channel will be set
%       to zero in the ratio image. Requires that both channels have masks.
%       Default is true.
%
%       ('CreateMasks' -> True/False)
%       If true, the intersection of the masks from both channels will be
%       used to create a new mask for the ratio channel and these masks
%       will be written to a mask directory for this channel.
%       Default is true.
%
%       ('MaskChannelIndex' -> Positive integer vector, 1x2) The index of
%       the channel of the masks to use for the numerator (1st element) and
%       denominator (2nd element), if masks are applied or created.
%       Default is to use the same channel as the images.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
% Output:
%
%   The ratio images will be written to a  new channel, named after the
%   channels used to make the ratios. All parameters will be stored in the
%   movieData structure.
%
% Hunter Elliott, 11/2009
%

%%  --------- Parameters ------- %%

pString = 'ratio_'; %The string to prepend before the ratio image directory & channel name

%% ---------- Input -------------%%

movieData = setupMovieData(movieData,'ratioing');

[batchMode,iNum,iDenom,applyMasks,makeMasks,iNumMask,iDenomMask] = parseInput(varargin);

%----Defaults----%

if isempty(batchMode)
    batchMode = false;
end
if isempty(applyMasks)
    applyMasks = true;
end
if isempty(makeMasks)
    makeMasks = true;
end
if isempty(iNum)
    if ~batchMode
        iNum = selectMovieChannels(movieData,0,'Please select the numerator channel:');
    else
        error('In batch mode, the numerator channel must be specified!')
    end
end
if isempty(iDenom)
    if ~batchMode
        iDenom = selectMovieChannels(movieData,0,'Please select the denominator channel:');
    else
        error('In batch mode, the denominator channel must be specified!')
    end
end
if makeMasks || applyMasks
    if isempty(iNumMask)
        iNumMask = iNum;
    end
    if isempty(iDenomMask)
        iDenomMask = iDenom;
    end
end


%% ------------- Init ------------ %%

disp('Initializing ratioing...')

if movieData.nImages(iNum) ~= movieData.nImages(iDenom)
    error('Numerator and denominator channels must have the same number of images!')
end
if movieData.imSize(:,iNum) ~= movieData.imSize(:,iDenom)
    error('Numerator and denominator channel images must be the same size!')
end
if (applyMasks || makeMasks) && ~checkMovieMasks(movieData,[iNumMask iDenomMask]);
    error('If masks are to be applied or created, both channels must have valid masks! Check masks or disable the CreateMasks and ApplyMasks options!')
end

nImages = movieData.nImages(iNum);

numImNames = getMovieImageFileNames(movieData,iNum);
denomImNames = getMovieImageFileNames(movieData,iDenom);

if applyMasks || makeMasks
    numMaskNames = getMovieMaskFileNames(movieData,iNumMask);
    denomMaskNames = getMovieMaskFileNames(movieData,iDenomMask);
end

%Check/Create the new channel

ratioChanName = [pString movieData.channelDirectory{iNum} '_to_' ...
                         movieData.channelDirectory{iDenom}];
                                          
iRatio = find(strcmp(ratioChanName,movieData.channelDirectory),1);

if isempty(iRatio)
    iRatio = length(movieData.channelDirectory) + 1;    
    movieData.channelDirectory{iRatio} = ratioChanName;
    movieData.nImages(iRatio) = nImages;
    movieData.imSize(:,iRatio) = movieData.imSize(:,iNum);
end

ratioDir = [movieData.imageDirectory filesep movieData.channelDirectory{iRatio}];

if ~exist(ratioDir,'dir')
    mkdir(ratioDir)
end


%Check/create channel for masks if requested
if makeMasks
    
    movieData.masks.channelDirectory{iRatio} = ratioChanName;
    maskDir = [movieData.masks.directory filesep movieData.masks.channelDirectory{iRatio}];
    
    if ~exist(maskDir,'dir')
        mkdir(maskDir)
    end    
end

%% ------ Pre-Ratio -----%
%Goes through all ratios and determines max&min values for writing to file
%by checking the maximum and minimum ratio values and the image bit-depth

%TEMP - This is unnecessary if there is enough memory to hold all images...
%Check total memory first???

disp('Pre-ratioing (calculating scale factor)...')

maxRatios = zeros(1,nImages);
minRatios = zeros(1,nImages);

for iImage = 1:nImages
    
    currNum = imread(numImNames{1}{iImage});    
   
    currDenom = imread(denomImNames{1}{iImage});
    
    if iImage == 1 %Get and check the image bit-depths
        ogClass = class(currNum);
        if ~strcmp(ogClass,class(currDenom))
            error('Numerator and denominator images must have the same bit-depth!')
        end        
    end
    
    currRatio = double(currNum) ./ double(currDenom);
    
    if applyMasks || makeMasks  %If masks are to be applied, don't include the masked values      
        currNumMask = imread(numMaskNames{1}{iImage});
        currDenomMask = imread(denomMaskNames{1}{iImage});
                
        if applyMasks
            currRatio(~currNumMask(:) | ~currDenomMask(:)) = 0;        
        end
         
    end    
    
    %Remove any infinities from division-by-zero
    currRatio(~isfinite(currRatio(:))) = 0;
    
    %Get max and min values for scalefactor calculation
    maxRatios(iImage) = max(currRatio(currRatio(:) ~= 0)); %Ignore zeros, in case masks have been applied/infs removed
    minRatios(iImage) = min(currRatio(currRatio(:) ~= 0));
    
    
end

%Calculate scale factor to minimize rounding error
scaleFactor = double(intmax(ogClass)) / (max(maxRatios)-min(minRatios));



%% ------ Ratio -----%%
% Ratios the channels and writes the resulting ratio image to a new channel

disp(['Creating ratio images by dividing channel ' movieData.channelDirectory{iNum} ...
      ' by channel ' movieData.channelDirectory{iDenom}])
  
disp(['Resulting images will be written to channel ' movieData.channelDirectory{iRatio}])
  
nDig = floor(log10(nImages))+1;
fString = ['%0' num2str(nDig) '.0f']; %Get formatting string for image names
  
for iImage = 1:nImages
    
     
    currNum = imread(numImNames{1}{iImage});
   
    currDenom = imread(denomImNames{1}{iImage});
    
    currRatio = double(currNum) ./ double(currDenom);
    
    if applyMasks  %If masks are to be applied, don't include the masked values      
        currNumMask = imread(numMaskNames{1}{iImage});
        currDenomMask = imread(denomMaskNames{1}{iImage});
        
        combMask = currNumMask & currDenomMask; %Get the intersection of the masks
        
        currRatio(~combMask) = 0; %Remove un-masked pixels        
        
        if makeMasks            
            imwrite(combMask,[maskDir filesep 'mask_' pString num2str(iImage,fString) '.tif'])
        end
    end    
    
    %Remove any infinities from division-by-zero (this shouldn't happen if
    %the masks are perfect, but let's be realistic here .... )
    currRatio(~isfinite(currRatio(:))) = 0;
  
    
    currRatio = cast(currRatio .* scaleFactor,ogClass);
    
    imwrite(currRatio,[ratioDir filesep pString num2str(iImage,fString) '.tif'])

    
end
   

%% ------ Output / Finalization ---- %%


movieData.ratioing.status = 1;
movieData.ratioing.dateTime = datestr(now);
movieData.ratioing.iFrom(1,iRatio) = iNum;
movieData.ratioing.iFrom(2,iRatio) = iDenom;
movieData.ratioing.scaleFactor(iRatio) = scaleFactor;
movieData.ratioing.applyMasks(iRatio) = applyMasks;
movieData.ratioing.createMasks(iRatio) = makeMasks; 
if makeMasks || applyMasks
    movieData.ratioing.iDenomMasksFrom(iRatio) = iDenomMask;
    movieData.ratioing.iNumMasksFrom(iRatio) = iNumMask;
end

updateMovieData(movieData);

disp('Finished!')













function [batchMode,iNum,iDenom,applyMasks,makeMasks,iNumMask,iDenomMask] = parseInput(argArray)


%Init output
batchMode = [];
iNum = [];
iDenom = [];
applyMasks = [];
makeMasks = [];
iNumMask = [];
iDenomMask = [];

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
           
       case 'Numerator'
           iNum = argArray{i+1};
           
       case 'Denominator'
           iDenom = argArray{i+1};
           
       case 'ApplyMasks'
           applyMasks = argArray{i+1};
           
       case 'CreateMasks'
           makeMasks = argArray{i+1};
           
       case 'ChannelIndex'
           iNum = argArray{i+1}(1);
           iDenom = argArray{i+1}(2);
           
       case 'MaskChannelIndex'
           iNumMask = argArray{i+1}(1);
           iDenomMask = argArray{i+1}(2);
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end