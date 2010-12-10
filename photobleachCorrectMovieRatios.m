function movieData = photobleachCorrectMovieRatios(movieData,paramsIn)
%PHOTOBLEACHCORRECTMOVIERATIOS applies a photobleach correction to ratio images for the input movie
% 
% movieData = photobleachCorrectMovieRatios(movieData)
% 
% movieData = photobleachCorrectMovieRatios(movieData,paramsIn)
%
% This function applies a photo-bleach correction to the ratio images in
% the input movie. Both the type of correction and the channel(s) to
% correct can be selected. 
%
% Input:
% 
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the corrected images to.
%       Corrected images for different channels will be saved as
%       sub-directories of this directory. If not input, the corrected
%       images will be saved to the same directory as the movieData, in a
%       sub-directory called "bleedthrough_corrected_images"
%
%       ('ChannelIndex'-> Positive integer scalar) The integer index of the
%       NUMERATOR of the ratio channel to perform photbleach correction on.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, the user will be asked to select
%       from the movie's channels.
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
%
% Output:
%
%   movieData - the updated movieData object with the correction
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The corrected images are written to the directory specified by the
%   parameter OuptuDirectory, with each channel in a separate
%   sub-directory. They will be stored as double-precision .mat files.
%
% 
% Hunter Elliott
% 11/2009
% Revamped 6/2010
%

%%  --------- Parameters ------- %%

pString = 'pbc_'; %The string to prepend before the corrected image directory & channel name
dName = 'photobleach_corrected_images_for_channel_'; %String for naming the directories for each corrected channel
fitFileName = 'photobleach_correction.mat'; %File name for saving fit results to
figName = 'photobleach correction fit.fig'; %Name for saving figure to file

%% ------ Input ------ %%


%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

nImages = movieData.nFrames_;

%Make sure there are enough frames
if nImages <= 4
    error('Input movie must have AT LEAST 5 timepoints for photobleach correction!!!')    
end

if nargin < 2
    paramsIn = [];
end

%Make sure the movie has been ratioed
iRProc = find(cellfun(@(x)(isa(x,'RatioProcess')),movieData.processes_),1);                          

if isempty(iRProc)
    error('The input movie has not been ratioed! Please perform ratioing prior to photobleach correction!')
end

%Get the indices of any previous photbleach correction processes from this
%function
iProc = find(cellfun(@(x)(isa(x,'PhotobleachCorrectionProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(PhotobleachCorrectionProcess(movieData,movieData.outputDirectory_));                                                                                                 
end


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

nChan = numel(movieData.channels_);

if isempty(p.ChannelIndex)
    if ~p.BatchMode
        p.ChannelIndex = selectMovieChannels(movieData,1,'Select the ratio channel to photobleach correct:');
    else
        error('In batch mode, you must specify the channel to photobleach correct!')
    end
end

if length(p.ChannelIndex) ~=1
    error('You can only photobleach-correct one ratio channel at a time!')
end

if p.ChannelIndex > nChan || p.ChannelIndex < 1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel number specified! Check ChannelIndex input!!')
end


%% -------- Init -------- %%

disp('Starting photobleach correction...')

%Get input directories/image names

iNum = movieData.processes_{iRProc}.funParams_.ChannelIndex(1);
iDenom = movieData.processes_{iRProc}.funParams_.ChannelIndex(2);

numDir = movieData.processes_{iRProc}.inImagePaths_{iNum};
numFileNames = movieData.processes_{iRProc}.getInImageFileNames(iNum);

denomDir = movieData.processes_{iRProc}. ...
    inImagePaths_{movieData.processes_{iRProc}.funParams_.ChannelIndex(2)};
denomFileNames = movieData.processes_{iRProc}.getInImageFileNames(iDenom);

ratDir = movieData.processes_{iRProc}. ...
    outImagePaths_{movieData.processes_{iRProc}.funParams_.ChannelIndex(1)};
ratioFileNames = movieData.processes_{iRProc}.getOutImageFileNames(iNum);

%Set-up output directory
outDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex)];

%Set up output directory
mkClrDir(outDir);
movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex,outDir);


nImTot = nImages*2;

%% ------- Calculate Intensity Vs. Time ----- %%
%Calculate the movie intensity vs. time in the needed channels

meanNum = zeros(1,nImages);
meanDenom = zeros(1,nImages);
meanRat = zeros(1,nImages);
totalNum = zeros(1,nImages);
totalDenom = zeros(1,nImages);

disp('Calculating average intensities...')

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, calculating intensity statistics...... ');        
end        



for iImage = 1:nImages
    
    currNum = double(imread([numDir filesep numFileNames{1}{iImage}]));
    meanNum(iImage) = mean(currNum(:));
    totalNum(iImage) = sum(currNum(:));
    currDenom = double(imread([denomDir filesep denomFileNames{1}{iImage}]));
    meanDenom(iImage) = mean(currDenom(:));
    totalDenom(iImage) = sum(currDenom(:));
    currRat = load([ratDir filesep ratioFileNames{1}{iImage}]);
    currRat = currRat.currRatio;
    meanRat(iImage) = mean(currRat(currRat(:) > 0));        

    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iImage / nImTot,wtBar)
    end                        
    
    
end




%% ----- Calculate Photobleach correction ----- %%

disp('Calculating fit...')


fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     %Double-exponential function for fitting
%Check if time was defined in moviedata
if ~isempty(movieData.timeInterval_)
    timePoints = (0:1:nImages-1) * movieData.timeInterval_;     %time data
else
    timePoints = (0:1:nImages-1);
end
bInit = [1 0 1 0]; %Initial guess for fit parameters.

 
switch p.CorrectionType


    case 'RatioOfAverages'

        fitData = meanNum ./ meanDenom;

    case 'AverageOfRatios' 

         fitData = meanRat;
         
    case 'RatioOfTotals'
        
        fitData = totalNum ./ totalDenom;

    otherwise

        error(['Invalid photobleach correction method!! "' p.CorrectionType '" is not a recognized method!'])
end


%Fit function to ratio timeseries
fitOptions = statset('Robust','on','MaxIter',500,'Display','off');
[bFit,resFit,jacFit,covFit,mseFit] = nlinfit(timePoints(:),fitData(:),fitFun,bInit,fitOptions);
%Get confidence intervals of fit and fit values
[fitValues,deltaFit] = nlpredci(fitFun,timePoints(:),bFit,resFit,'covar',covFit,'mse',mseFit);

%Check the fit jacobian
[dummy,R] = qr(jacFit,0); %#ok<ASGLU>
if ~p.BatchMode && condest(R) > 1/(eps(class(bFit)))^(1/2)        
    warndlg('WARNING: The photobleach correction fit is not very good. Please use extreme caution in interpreting ratio changes over time in the photobleach corrected ratios!')
end


%% ----- Apply photobleach correction to ratio images ----- %%


disp(['Applying photobleach correction method ' p.CorrectionType ' to ratio channel ' ratDir ])
disp(['Writing corrected images to channel ' outDir])

%Disable convert-to-integer warning
warning('off','MATLAB:intConvertNonIntVal');

if ~p.BatchMode        
    waitbar(nImages / nImTot,wtBar,'Please wait, applying photobleach correction ...');
end        


%Go through all the images and correct them
for iImage = 1:nImages
   
    %Load the image
    currRat = load([ratDir filesep ratioFileNames{1}{iImage}]);
    currRat = currRat.currRatio;
    
    %Correct the image. We multiply by the average of the first ratio to
    %prevent normalization.
    currRat = currRat ./ fitValues(iImage) .* meanRat(1); %#ok<NASGU>
    
    %Write it back to file.    
    save([outDir filesep pString ratioFileNames{1}{iImage}],'currRat');
    
    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar((iImage +nImages)/ nImTot,wtBar)
    end                        

end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------- Make and Save Figure ------- %%


disp('Making figures...')

if p.BatchMode
    fitFig = figure('Visible','off');
else
    fitFig = figure;
end

hold on
title('Photobleach Correction Fit')
if ~isempty(movieData.timeInterval_)
    xlabel('Time, seconds')
else
    xlabel('Frame Number')
end
ylabel(p.CorrectionType)
plot(timePoints,fitData)
plot(timePoints,fitValues,'r')
plot(timePoints,fitValues+deltaFit,'--r')
legend(p.CorrectionType,'Fit','Fit 95% C.I.')
plot(timePoints,fitValues-deltaFit,'--r')

hgsave(fitFig,[p.OutputDirectory filesep figName]);
%Log this file name in the parameter structure
p.figName = figName;
movieData.processes_{iProc}.setPara(p);


if ishandle(fitFig) %make sure user hasn't closed it.
    close(fitFig)
end


%% ----- Output/Finalization ---- %%

save([p.OutputDirectory filesep fitFileName],'fitData','fitValues',...
    'timePoints','covFit','mseFit','resFit','jacFit',...
    'fitFun','bFit');

%Log the correction in the movieData object and save it

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData;


disp('Finished!')
