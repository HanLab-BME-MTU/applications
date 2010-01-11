function movieData = photobleachCorrectMovieRatios(movieData,varargin)

% 
% movieData = photobleachCorrectFRETmovie(movieData)
% 
% movieData = photobleachCorrectFRETmovie(movieData,'OptionName',optionValue)
%
% This function applies a photo-bleach correction to the ratio images in
% the input movie. Both the type of correction and the channel(s) to
% correct can be selected. 
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
%
%   Possible Option Names:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The
%       integer indices of the RATIO channel(s) to perform photobleach
%       correction on. This index corresponds to the channel directories
%       location in the cell array movieData.channelDirectory. If not
%       input, the user will be asked to select from the available
%       channels. (Unless in batch mode, in which case this MUST be input)
%
%       ('CorrectionType' -> character string) Character string describing
%       the photo-bleach correction method to use. The options are:
%
%           "RatioOfAverages" The average intensity of the images used to
%           make the ratios is calculated in each Frame. Then these averages
%           are ratioed, and a double-exponential is fit to the resulting
%           timeseries. 
%
%           "AverageOfRatios" The average masked value of the ratio images is
%           calculated in each frame, and then a double-exponential is fit to
%           this timeseries.
%
%           "RatioOfTotals" The total intensity of the images used to make
%           the ratios is calculated in each frame, and a double-eponential
%           is fit to the ratio of these values.
%           NOTE: In all cases, the correction is applied to the ratio images.
%         
%       ('BatchMode' -> True/False) If true, all graphical outputs and user
%       interaction is suppressed. 
% 
% 
% Output:       
% 
%   The corrected images will be written to a new image channel in the
%   movie, and all parameters will be stored in the movieData structure.
% 
% Hunter Elliott, 11/2009
% 
%%  --------- Parameters ------- %%

pString = 'pbc_'; %The string to prepend before the corrected image directory & channel name
fitFileName = 'photobleach_correction.mat'; %File name for saving fit results to

%% ------ Input ------ %%

movieData = setupMovieData(movieData,'photobleachCorrection');

[batchMode,pbCorrectType,iChannel] = parseInput(varargin);

%---Defaults--%

if isempty(batchMode)
    batchMode = false;
end
if isempty(pbCorrectType)
    pbCorrectType = 'RatioOfAverages';
end
if isempty(iChannel)
    if ~batchMode
        iChannel = selectMovieChannels(movieData,1,'Select the ratio channel to photobleach correct:');
    else
        error('In batch mode, you must specify the channel to photobleach correct!')
    end
end



%% -------- Init -------- %%

disp('Starting photobleach correction...')

%Make sure this is a ratio channel 
if ~checkMovieProcedure(movieData,'ratioing')
    error('This function only photobleach corrects ratio images! Create a ratio channel using ratioMovie.m before running!')
end
%Get the num/denom index
iNum = movieData.ratioing.iFrom(1,iChannel);
iDenom = movieData.ratioing.iFrom(2,iChannel);
if iNum < 1 || iDenom < 1
    error('Invalid ratio channel! The specified channel must have been created using ratioMovie.m!')
end

%Calculate the movie intensity vs. time. This is re-run every in case
%changes have been made (its fast anyways).
chanUsed = [iChannel iNum iDenom]; %List of all channels used in this function
nChanUsed = length(chanUsed);
movieData = calculateMovieIntensityVsTime(movieData,'ChannelIndex',chanUsed,'BatchMode',1);

%Make sure there are enough frames
if movieData.nImages(iChannel) <= 4
    error('Specified channel must have AT LEAST 5 timepoints for photobleach correction!!!')    
end
nImages = movieData.nImages(iChannel);

%Load the intensity statistics for each of the required channels
allMeanI = cell(1,nChanUsed);
allTotalI = cell(1,nChanUsed);
for i = 1:nChanUsed
    
    clear('totalIntensity','meanIntensity') %Erase the variables from the previous channel
    
    load([movieData.intensityVsTime.directory filesep  ...
        movieData.intensityVsTime.fileNamePrefix movieData.channelDirectory{chanUsed(i)}])
    
    if ~(exist('totalIntensity','var') && exist('meanIntensity','var'))
        error(['Problem with intensity vs time file for channel ' ... 
            movieData.channelDirectory{chanUsed(i)} ... 
            '! should contain variables named totalIntensity and meanIntensity!'])
    end
        
    
    allMeanI{i} = meanIntensity;
    allTotalI{i} = totalIntensity;           
    
end



%Check/create channel for resulting images
chanName = [pString movieData.channelDirectory{iChannel}];

iCorrected = find(strcmp(chanName,movieData.channelDirectory),1);

if isempty(iCorrected)
    iCorrected = length(movieData.channelDirectory)+1;
    movieData.channelDirectory{iCorrected} = chanName;
    movieData.nImages(iCorrected) = movieData.nImages(iChannel);
    movieData.imSize(:,iCorrected) = movieData.imSize(iChannel);
end

corrDir = [movieData.imageDirectory filesep movieData.channelDirectory{iCorrected}];

if ~exist(corrDir,'dir')
    mkdir(corrDir)
end


%Get ratio image file names
imageFileNames = getMovieImageFileNames(movieData,iChannel);


%% ----- Calculate Photobleach correction ----- %%

disp('Calculating fit...')


fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     %Double-exponential function for fitting
timePoints = (0:1:nImages-1) * movieData.timeInterval_s;     %time data
bInit = [1 0 1 0]; %Initial guess for fit parameters.

 
switch pbCorrectType


    case 'RatioOfAverages'

        fitData = allMeanI{chanUsed == iNum} ./ allMeanI{chanUsed == iDenom};

    case 'AverageRatios' 

         fitData = allMeanI{chanUsed == iChannel};
         
    case 'RatioOfTotals'
        
        fitData = allTotalI{chanUsed == iNum} ./ allTotalI{chanUsed == iDenom};

    otherwise

        error(['Invalid photobleach correction method!! "' pbCorrectType '" is not a recognized method!'])
end

%If there was a stimulation event, remove the points after it.
if isMovieStimulated(movieData)
    storeFitData = fitData; %save the removed points for later plotting
    fitData(movieData.stimulation.iFrame:end) = NaN;
end

%Fit function to ratio timeseries
fitOptions = statset('Robust','off','MaxIter',500,'Display','off');
[bFit,resFit,jacFit,covFit,mseFit] = nlinfit(timePoints(:),fitData(:),fitFun,bInit,fitOptions); %#ok<NASGU>

%Evaluate the fitted function at each timepoint
fitValues = fitFun(bFit,timePoints);


%% ----- Apply photobleach correction to ratio images ----- %%


disp(['Applying photobleach correction method ' pbCorrectType ' to ratio channel ' movieData.channelDirectory{iChannel}])
disp(['Writing corrected images to channel ' movieData.channelDirectory{iCorrected}])

%Disable convert-to-integer warning
warning('off','MATLAB:intConvertNonIntVal');

%Go through all the images and correct them
for iImage = 1:nImages
   
    %Load the image
    currIm = imread(imageFileNames{1}{iImage});
    
    %Correct the image
    currIm = currIm ./ fitValues(iImage);
    
    %Write it back to file.
    iLastSep = max(regexp(imageFileNames{1}{iImage},filesep));
    imwrite(currIm,[corrDir filesep pString ...
        imageFileNames{1}{iImage}(iLastSep+1:end)]);    
    
end



%% ------- Make and Save Figure ------- %%


disp('Making figures...')

if batchMode
    fitFig = figure('Visible','off');
else
    fitFig = figure;
end

hold on
title('Photobleach Correction Fit')
xlabel('Time, seconds')
ylabel(pbCorrectType)
plot(timePoints,fitData)
plot(timePoints,fitValues,'k')
legend(pbCorrectType,'Fit')

if isMovieStimulated(movieData)
    plot(timePoints,storeFitData,':')
end

hgsave(fitFig,[movieData.photobleachCorrection.directory filesep 'photobleach correction fit.fig']);

if batchMode
    close(fitFig)
end


%% ----- Output/Finalization ---- %%

save([movieData.photobleachCorrection.directory filesep fitFileName],'fitData','fitValues',...
    'timePoints','covFit','mseFit','resFit','jacFit',...
    'fitFun','bFit');


movieData.photobleachCorrection.status = 1;
movieData.photobleachCorrection.iFrom(iCorrected) = iChannel;
movieData.photobleachCorrection.dateTime = datestr(now);
movieData.photobleachCorrection.method = pbCorrectType;
movieData.photobleachCorrection.fileName = fitFileName; 

updateMovieData(movieData);


disp('Finished!')

function [batchMode,pbCorrectType,iChannel] = parseInput(argArray)


%Init output
batchMode = [];
pbCorrectType = [];
iChannel = [];

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
           
       case 'CorrectionType'
           pbCorrectType = argArray{i+1};
           
       case 'ChannelIndex'
           iChannel = argArray{i+1};
 
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end
