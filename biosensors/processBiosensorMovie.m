function runInfo = processBiosensorMovie(movieData,numerator,denominator,varargin)
%PROCESSBIOSENSORMOVIE runs via command line all necessary processing steps on the input biosensor movie 
% 
% runInfo = processBiosensorMovie(movieData,numerator,denominator)
% runInfo = processBiosensorMovie(...,'OptionName1',optionValue1,'OptionName2,optionValue2
% runInfo = processBiosensorMovie(...,optionStruc); 
%
% This function runs all necessary processing steps on the input biosensor
% movie, including segmentation, all corrections and ratioing. 
% Right now this function is very simple and crude, so don't expect too
% much. I'll probably improve it later.... probably.
% 
% **WARNING** This function has the following known, damning flaws: 
%       -Does not currently support bleedthrough correction
%       -Does not create/use the BiosensorsPackage object
%       -Has tons of hard-coded parameters
%       -Does not handle process dependencies or errors in any way - all steps are
%       run every time this is called, and if one step errors the whole
%       processing stops.
%       -Assumes a transformation is required.
% 
% **NOTE:*** This is the non-user-friendly, command-line-only function for
% processing - for a more user-friendly function which is point-and-click
% GUI based, please call movieSelectorGUI and then select the Biosensors
% Package option. Also, there's a manual for the software the inside the
% "Doc" folder, so read that shit or risk doing bad science and annoying
% people.
% 
% Input:
% 
%   movieData - The MovieData object describing the movie to process.
% 
%   numerator - A string (character array) identifying the channel which
%   will be the numerator of the ratio. This string should match the name
%   of the folder containing the images for this channel.
% 
%   denominator - A string (character array) identifying the channel which
%   will be the numerator of the ratio. This string should match the name
%   of the folder containing the images for this channel.
%   
%   'OptionName',optionValue - Name and value for additional inputs. May
%   also be input as a structure. * indicates required parameter
%
%       ('OptionName'->optionValue)
% 
%       ('ShadeImageDirectories'->Cell array of character arrays) This cell
%       array specifies the folder(s) containing the shade-correction
%       images for each channel. It must have a directory for each channel
%       to be corrected, though multiple channels may use the same
%       corrections. If not specified, the user will be asked.
%
%       ('DarkImageDirectories' -> Cell array of character strings)
%       This cell array contains directories for dark-current images which
%       should be used to correct the raw images. Must contain one valid
%       directory for each channel specified by ChannelIndex.
%       Optional. If not input, the user will be asked to select a folder
%       for each channel.
%
%       ('TransformFilePaths' -> Cell array of Character strings) A cell
%       array specifying The FULL path and filename of the .mat file
%       containing the transform to apply to the images in each channel.
%       Should contain one element for each channel to be transformed. The
%       transform should be of the format used by imtransform.m. If not
%       input, the user will be asked to locate a file containing a
%       transform for each channel, UNLESS batchmode is enabled, in which
%       case an error will be generated.
%
%       *('TransformChannel'->'numerator' or 'denominator') Specifies which
%       channel to apply the transformation to. Required input.
% 
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
%
% Output:
% 
%   runInfo - A cell array the same size as the number of processing steps,
%   with some info regarding each step including the parameters which were
%   passed to each processing function
%
% Hunter Elliott
% 8/2011
%

%% --------------------- Input ---------------------- %%

%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.KeepUnmatched = true; %Keep extra parameters for passing to processing functions
ip.addRequired('movieData',@(x)(isa(x,'MovieData')));
ip.addRequired('numerator',@(x)(ischar(x) && ~isempty(x)));
ip.addRequired('denominator',@(x)(ischar(x) && ~isempty(x)));

ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.addParamValue('ShadeImageDirectories',{},@iscell);
ip.addParamValue('DarkImageDirectories',{},@iscell);
ip.addParamValue('TransformFilePaths',{},@iscell);
ip.addParamValue('TransformChannel',[],@ischar);

%ip.addParamValue('ForceRun',0,(@(x)(numel(x) == 1 || numel(x) == nSteps)));

ip.parse(movieData,numerator,denominator,varargin{:});
%procP = ip.Unmatched;%Unrecognized parameters will be passed to processing functions

p = ip.Results;


% ------------- Basic Checking --------------- %

nChan = numel(movieData.channels_);

if nChan < 2
    error('The input movie must have at least 2 channels!')
end



%% ---------------------- Init -------------------- %%

%Convert the channel name into a channel index
chanPaths = arrayfun(@(x)(x.channelPath_),movieData.channels_,'UniformOutput',false);%Get channel paths
chanPaths = cellfun(@(x)(x(max(regexp(x,filesep))+1:end)),chanPaths,'UniformOutput',false);%Retain only the image folder itself
numInd = find(strcmp(numerator,chanPaths));%Index of numerator channel
denInd = find(strcmp(denominator,chanPaths));%Index of denominator channel

%Make sure exactly one index was found for both channels
if isempty(numInd)
    error('PROCBIOSENSMOV:numChanMatch',['No channel present in movie which matches numerator string "' numerator '"'])
elseif numel(numInd) > 1
    error('PROCBIOSENSMOV:numChanMatch',['More than one channel matches the numerator string "' numerator '"'])
end
if isempty(denInd)
    error(['No channel present in movie which matches denominator string "' denominator '"'])
elseif numel(denInd) > 1
    error('PROCBIOSENSMOV:denChanMatch',['More than one channel matches the denominator string "' denominator '"'])
end

if ~isempty(p.TransformChannel)
    
    if strcmpi(p.TransformChannel,'numerator')
        p.TransformChannel = numInd;
    elseif strcmpi(p.TransformChannel,'denominator')
        p.TransformChannel = denInd;
    else
        error('PROCBIOSENSMOV:numChanMatch','The input TransformChannel must be either "numerator" or "denominator"!')
    end    
else
    error('PROCBIOSENSMOV:numChanMatch','Required option "TransformChannel" was not input!')
end
    


%Structure with parameters which are general to every processing step
genP.ChannelIndex = [numInd denInd];
genP.BatchMode = p.BatchMode;


iStep = 1;

disp('Starting biosensor movie processing...')

%% --------------------- Processing ---------------- %%
% Go through each step and run it. Eventually this should be made more
% general - converted to a loop rather than a series of commands....


%Set up parameters and run the thresholding process
sp = genP;%Get the general parameters
sp.MaxJump = 1.1; %Use threshold-jump suppression
movieData = thresholdMovie(movieData,sp);
threshProc = movieData.getProcessIndex('ThresholdingProcess',1,0);%Get the thresholding process index for later use
runInfo{iStep} = sp;
iStep = iStep + 1;

%Set up parameters and run the background mask generation
sp = genP;%Get the general parameters
sp.SegProcessIndex = threshProc;
movieData = createMovieBackgroundMasks(movieData,sp);
runInfo{iStep} = sp;
iStep = iStep + 1;

%Set up parameters and run the mask refinement process
sp = genP;%Get the general parameters
sp.SegProcessIndex = threshProc;
movieData = refineMovieMasks(movieData,sp);
refineProc = movieData.getProcessIndex('MaskRefinementProcess',1,0);%Get the refinement process index for later use

%Run dark-current correction if directories were input
if ~isempty(p.DarkImageDirectories)
    sp = genP;
    sp.DarkImageDirectories = p.DarkImageDirectories;
    movieData = darkCurrentCorrectMovie(movieData,sp);
    runInfo{iStep} = sp;
    iStep = iStep + 1;            
else
    disp('Not running dark-current correction - no dark image directory specified.')
end
   
%Set up and run shade correction
sp = genP;
sp.ShadeImageDirectories = p.ShadeImageDirectories;
movieData = shadeCorrectMovie(movieData,sp);
runInfo{iStep} = sp;
iStep = iStep + 1;            

%background subtraction
sp = genP;
movieData = backgroundSubtractMovie(movieData,sp);

%transformation
sp = genP;
sp.ChannelIndex = p.TransformChannel;
sp.TransformFilePaths = p.TransformFilePaths;
sp.SegProcessIndex = refineProc;
movieData = transformMovie(movieData,sp);
maskXFProc = movieData.getProcessIndex('MaskTransformationProcess',1,0);%Get the mask transform process index for later use
runInfo{iStep} = sp;
iStep = iStep + 1;            

%ratio
sp = genP;
sp.SegProcessIndex = [refineProc maskXFProc];
movieData = ratioMovie(movieData,sp);
runInfo{iStep} = sp;
iStep = iStep + 1;            

%Photobleach correction
sp = genP;
sp.ChannelIndex = numInd;
movieData = photobleachCorrectMovieRatios(movieData,sp);
runInfo{iStep} = sp;

disp('Finished biosensor movie processing...')

















